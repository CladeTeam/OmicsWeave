"""CENOMultiome training data with a plain PyTorch DataLoader.

Illustrative: paths and the manifest are placeholders, and the data are
taken to be converted already (``ow convert``), one set per cell:

- reference: bgzip FASTA ``hg38.fa.gz`` with its ``.fai`` / ``.gzi``
- bp tracks (bigWig): ATAC, DNase, CAGE +/-, RNA polyA +/-, RNA total +/-
- coarse tracks (bigWig): 11 histone marks and CTCF, pooled to 128 bp
- Hi-C: ``.hic`` (from ``.mcool``), read at 2048 bp bins
- expression: one CellPack matrix (``.cpk``), a row per cell, read with
  ``pycellpack``

A sample is one window of one cell: the window's bases (shared by every
cell, so read once per window and batch), each assay's signal with a valid
mask (an assay a cell lacks is all zeros, valid 0), the Hi-C tile, the
cell's expression, and which observed assays are inputs and which are
held out as targets. Windows are fixed anchors shifted at random, reverse
complemented half the time with the +/- strands swapped.

Readers open once per worker, on its first batch, and answer every query
after; the batch sampler runs in the main process and groups the cells of
one window, which ``__getitems__`` serves together.
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import torch
from pycellpack import Store
from torch.utils.data import DataLoader, Dataset, Sampler

from omicsweave.formats import BgzipFasta, BigWig, Hic
from omicsweave.query import ContactReader, SequenceReader, TrackReader
from omicsweave.types import ContactWindow, Region

BP_ASSAYS = (
    "atac",
    "dnase",
    "cage_plus",
    "cage_minus",
    "rna_polya_plus",
    "rna_polya_minus",
    "rna_total_plus",
    "rna_total_minus",
)
COARSE_ASSAYS = (
    "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac",
    "H3K27me3", "H3K36me3", "H3K79me2", "H4K20me1", "H2AFZ", "CTCF",
)  # fmt: skip
# Each bp assay's strand partner (itself when unstranded): reverse complement
# swaps them.
STRAND_SWAP = [0, 1, 3, 2, 5, 4, 7, 6]

WINDOW = 131_072
COARSE = 128  # bp per coarse-track bin
HIC_BIN = 2048  # bp per Hi-C bin: a 64 x 64 tile per window
SHIFT = 1024  # largest anchor shift, either way


@dataclass(frozen=True, slots=True)
class Cell:
    """One cell's files; an assay it lacks is ``None``."""

    name: str
    row: int  # its row in the expression matrix
    bp: Mapping[str, BigWig | None]
    coarse: Mapping[str, BigWig | None]
    hic: Hic | None


def load_cells(manifest: str) -> list[Cell]:
    """Cells from a JSON manifest: ``[{"name", "row", "bp": {assay: path},
    "coarse": {assay: path}, "hic": path}, ...]``."""
    with open(manifest) as f:
        entries = json.load(f)
    return [
        Cell(
            e["name"],
            e["row"],
            {a: BigWig(e["bp"][a]) if a in e["bp"] else None for a in BP_ASSAYS},
            {a: BigWig(e["coarse"][a]) if a in e["coarse"] else None for a in COARSE_ASSAYS},
            Hic(e["hic"]) if e.get("hic") else None,
        )
        for e in entries
    ]


@dataclass(frozen=True, slots=True)
class Key:
    """What one sample is: a window of one cell, and its augmentation."""

    region: Region
    cell: int
    reverse: bool
    targets: tuple[int, ...]  # bp assays held out of the input


class AnchorBatches(Sampler[list[Key]]):
    """Batches of ``cells_per_window`` cells of one shifted anchor window;
    each sample hides a random share of its observed bp assays as targets."""

    def __init__(
        self,
        anchors: Sequence[Region],
        cells: Sequence[Cell],
        *,
        cells_per_window: int = 8,
        target_share: float = 0.3,
        seed: int = 0,
    ):
        self.anchors, self.cells = list(anchors), list(cells)
        self.cells_per_window, self.target_share, self.seed = cells_per_window, target_share, seed
        self.epoch = 0

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def __iter__(self) -> Iterator[list[Key]]:
        rng = np.random.default_rng(self.seed + self.epoch)
        for i in rng.permutation(len(self.anchors)):
            anchor = self.anchors[i]
            shift = int(rng.integers(-SHIFT, SHIFT + 1))
            region = Region(anchor.chrom, anchor.start + shift, anchor.end + shift)
            reverse = bool(rng.random() < 0.5)
            batch = []
            for c in rng.choice(len(self.cells), self.cells_per_window, replace=False):
                observed = [a for a, f in enumerate(self.cells[c].bp.values()) if f is not None]
                n = int(round(self.target_share * len(observed)))
                targets = tuple(sorted(rng.choice(observed, n, replace=False).tolist()))
                batch.append(Key(region, int(c), reverse, targets))
            yield batch

    def __len__(self) -> int:
        return len(self.anchors)


_COMPLEMENT = np.arange(256, dtype=np.uint8)
for a, b in (b"AT", b"TA", b"CG", b"GC", b"at", b"ta", b"cg", b"gc"):
    _COMPLEMENT[a] = b


def _flip(a: np.ndarray) -> np.ndarray:
    """``a`` reversed along its position axis (the last)."""
    return np.ascontiguousarray(a[..., ::-1])


@dataclass
class _Readers:
    sequence: SequenceReader
    bp: list[list[TrackReader | None]]  # [cell][assay]
    coarse: list[list[TrackReader | None]]
    hic: list[ContactReader | None]
    scale: np.ndarray  # [cell, assay]: genome-wide mean of each bp track
    expression: Store


def _track(file: BigWig | None) -> TrackReader | None:
    return None if file is None else TrackReader(file)


def _mean(reader: TrackReader | None) -> float:
    s = None if reader is None else reader.summary
    return 1.0 if s is None or s.bases_covered == 0 else s.sum / s.bases_covered


class CenoWindows(Dataset):
    """Samples keyed by :class:`Key`; served a batch at a time."""

    def __init__(self, reference: BgzipFasta, cells: Sequence[Cell], expression: str):
        self.reference, self.cells, self.expression = reference, list(cells), expression
        self._readers: _Readers | None = None

    def _open(self) -> _Readers:
        # Opened in the worker on its first batch; every later query reuses
        # them. 54 cells x 20 tracks stay open per worker.
        if self._readers is None:
            bp = [[_track(f) for f in c.bp.values()] for c in self.cells]
            self._readers = _Readers(
                sequence=SequenceReader(self.reference),
                bp=bp,
                coarse=[[_track(f) for f in c.coarse.values()] for c in self.cells],
                hic=[None if c.hic is None else ContactReader(c.hic) for c in self.cells],
                scale=np.array([[_mean(r) for r in row] for row in bp], np.float32),
                expression=Store(self.expression).astype(np.float32),
            )
        return self._readers

    @staticmethod
    def _signal(
        readers: Sequence[TrackReader | None], region: Region
    ) -> tuple[np.ndarray, np.ndarray]:
        values = np.zeros((len(readers), len(region)), np.float32)
        valid = np.zeros((len(readers), len(region)), np.bool_)
        for row, ok, reader in zip(values, valid, readers, strict=True):
            if reader is not None:
                reader.signal(region, out=row)
                ok[:] = ~np.isnan(row)
        np.nan_to_num(values, copy=False)
        return values, valid

    def _sample(self, key: Key, bases: np.ndarray, expression: np.ndarray) -> dict[str, np.ndarray]:
        r = self._readers
        assert r is not None
        bp, bp_valid = self._signal(r.bp[key.cell], key.region)
        bp /= r.scale[key.cell][:, None]
        coarse, coarse_valid = self._signal(r.coarse[key.cell], key.region)
        coarse = coarse.reshape(len(COARSE_ASSAYS), -1, COARSE).mean(-1)
        coarse_valid = coarse_valid.reshape(len(COARSE_ASSAYS), -1, COARSE).all(-1)

        # The tile is read on the 2048 bp grid nearest the shifted window, so
        # it always holds WINDOW / HIC_BIN bins.
        start = key.region.start // HIC_BIN * HIC_BIN
        grid = Region(key.region.chrom, start, start + WINDOW)
        window = ContactWindow(grid, grid, HIC_BIN)
        hic = r.hic[key.cell]
        contacts = np.zeros(window.shape, np.float32)
        if hic is not None:
            hic.matrix(window, out=contacts)  # raw counts, in place
            np.log1p(contacts, out=contacts)

        observed = np.array([f is not None for f in self.cells[key.cell].bp.values()])
        target = np.zeros(len(BP_ASSAYS), np.bool_)
        target[list(key.targets)] = True
        inputs = bp * (observed & ~target)[:, None]

        if key.reverse:
            bases = _COMPLEMENT[bases[::-1]]
            inputs, bp, bp_valid = (_flip(a[STRAND_SWAP]) for a in (inputs, bp, bp_valid))
            coarse, coarse_valid = _flip(coarse), _flip(coarse_valid)
            contacts = np.ascontiguousarray(contacts[::-1, ::-1])
            observed, target = observed[STRAND_SWAP], target[STRAND_SWAP]
        return {
            "sequence": bases,  # [L] uint8 ASCII tokens
            "inputs": inputs,  # [8, L] observed bp assays, targets zeroed
            "targets": bp,  # [8, L]
            "valid": bp_valid,  # [8, L]
            "observed": observed,  # [8]
            "target": target,  # [8]
            "coarse": coarse,  # [12, L/128]
            "coarse_valid": coarse_valid,
            "hic": contacts,  # [64, 64] log1p counts
            "has_hic": np.array(hic is not None),
            "expression": expression,  # [genes]
        }

    def __getitems__(self, keys: list[Key]) -> list[dict[str, np.ndarray]]:
        r = self._open()
        # One read per distinct window: the batch's cells share its bases.
        bases = {k.region: r.sequence.fetch(k.region).bases for k in keys}
        # One CellPack decode for the batch's cells.
        expression = r.expression[[self.cells[k.cell].row for k in keys]].read_dense()
        return [self._sample(k, bases[k.region], e) for k, e in zip(keys, expression, strict=True)]

    def __getitem__(self, key: Key) -> dict[str, np.ndarray]:
        return self.__getitems__([key])[0]


def train_step(batch: dict[str, torch.Tensor]) -> None:
    """The model's forward and backward pass (not shown)."""


def anchors_of(chroms: Mapping[str, int], keep: Sequence[str]) -> list[Region]:
    """Non-overlapping anchors over the chromosomes of one split, clear of
    the ends by one shift."""
    return [
        Region(c, s, s + WINDOW)
        for c in keep
        for s in range(SHIFT, chroms[c] - WINDOW - SHIFT, WINDOW)
    ]


def main() -> None:
    reference = BgzipFasta("/data/hg38/hg38.fa.gz")
    cells = load_cells("/data/ceno/cells.json")  # 54 cells
    with SequenceReader(reference) as fa:
        chroms = dict(fa.chroms)
    train = [f"chr{i}" for i in range(1, 23) if i not in (8, 9)]

    batches = AnchorBatches(anchors_of(chroms, train), cells, cells_per_window=8, seed=7)
    loader = DataLoader(
        CenoWindows(reference, cells, "/data/ceno/expression.cpk"),
        batch_sampler=batches,  # one window x 8 cells per batch
        num_workers=16,
        prefetch_factor=2,  # batches each worker keeps ready
        persistent_workers=True,  # readers stay open across epochs
        pin_memory=True,
    )
    for epoch in range(2):
        batches.set_epoch(epoch)  # main process: workers just receive keys
        for batch in loader:
            # sequence [8, L], inputs / targets / valid [8, 8, L],
            # coarse [8, 12, L / 128], hic [8, 64, 64], expression [8, genes]
            train_step(batch)


if __name__ == "__main__":
    main()
