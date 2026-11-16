"""NonB (G4 / R-loop) training data with a plain PyTorch DataLoader.

Illustrative: paths and manifests are placeholders. Every file is read by
region in a format built for it; text formats are converted once
beforehand:

- reference: bgzip FASTA ``hg38.fa.gz`` with its ``.fai`` / ``.gzi``
  (``ow convert hg38.fa hg38.fa.gz``)
- signal: one bigWig per track (CUT&Tag G4 / R-loop, ...), up to 146;
  bedGraph or wig is converted
  (``ow convert G4_K562_rep1.bedGraph G4_K562_rep1.bw --chrom-sizes hg38.fa.gz.fai``)
- peaks: one bigBed per track, plus their union ``union_peaks.bb``; BED is
  not read by region and is converted
  (``ow convert G4_K562.bed G4_K562.bb --chrom-sizes hg38.fa.gz.fai``)

Two pipelines:

- tracks: 131 kb windows tiling each chromosome, a fresh random overlap
  each epoch (``seed + epoch``); bp signal of every track, normalized by
  its mean over covered bases and averaged into ``resolution`` bins.
  Tracks without a positive, finite header mean keep their original scale.
- peaks: 5 kb one-hot windows for binary classification; positives centred
  on peaks, negatives drawn where no peak lies; per-track 0/1 bin coverage.

Readers open once per worker, on its first sample, and answer every query
after. The sampler runs in the main process and hands out regions, so a new
epoch's windows reach persistent workers. The track pipeline gives each
worker one block cache shared by its readers, so overlapping windows reuse
the reference and bigWig blocks they share.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset, Sampler

from omicsweave.cache import BlockCache
from omicsweave.formats import BgzipFasta, BigBed, BigWig
from omicsweave.query import IntervalReader, SequenceReader, TrackReader
from omicsweave.types import Region

# Chromosome splits: whole chromosomes, so neighbouring windows never leak.
SPLITS = {
    "train": [f"chr{i}" for i in range(1, 23) if i not in (8, 9)] + ["chrX"],
    "val": ["chr8"],
    "test": ["chr9"],
}

# ACGT -> 0..3, anything else (N, IUPAC) -> 4, which one-hot encodes as zeros.
_CODE = np.full(256, 4, np.uint8)
for i, base in enumerate(b"ACGT"):
    _CODE[base] = _CODE[base + 32] = i
_ONE_HOT = np.vstack([np.eye(4, dtype=np.float32), np.zeros(4, np.float32)])


def one_hot(bases: np.ndarray) -> np.ndarray:
    """``(L, 4)`` float32 one-hot of ASCII bases; N is all zeros."""
    return _ONE_HOT[_CODE[bases]]


def bin_coverage(start: np.ndarray, end: np.ndarray, region: Region, resolution: int) -> np.ndarray:
    """0/1 per ``resolution``-base bin of ``region``: covered by any interval."""
    nbins = -(-len(region) // resolution)
    # Columns are uint32: widen before subtracting.
    start = np.clip(start.astype(np.int64), region.start, region.end) - region.start
    end = np.clip(end.astype(np.int64), region.start, region.end) - region.start
    nonempty = start < end
    first, last = start[nonempty] // resolution, -(-end[nonempty] // resolution)
    edges = np.zeros(nbins + 1, np.int32)
    np.add.at(edges, first, 1)
    np.add.at(edges, last, -1)
    return (np.cumsum(edges[:-1]) > 0).astype(np.float32)


# ---------------------------------------------------------------- tracks


class TiledWindows(Sampler[Region]):
    """Windows tiling ``chroms``; each epoch starts the tiling at a random
    offset and jitters every window, so neighbours overlap differently."""

    def __init__(
        self, chroms: Mapping[str, int], length: int, *, seed: int = 0, shuffle: bool = True
    ):
        self.chroms, self.length, self.seed, self.shuffle = dict(chroms), length, seed, shuffle
        self.epoch = 0

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def windows(self) -> list[Region]:
        rng = np.random.default_rng(self.seed + self.epoch)
        out = []
        for chrom, size in self.chroms.items():
            if size < self.length:
                continue
            offset = int(rng.integers(min(self.length, size - self.length + 1)))
            starts = np.arange(offset, size - self.length + 1, self.length)
            jitter = rng.integers(-self.length // 4, self.length // 4 + 1, len(starts))
            starts = np.clip(starts + jitter, 0, size - self.length)
            out += [Region(chrom, int(s), int(s) + self.length) for s in starts]
        if self.shuffle:
            rng.shuffle(out)
        return out

    def __iter__(self) -> Iterator[Region]:
        return iter(self.windows())

    def __len__(self) -> int:
        return len(self.windows())


class TrackWindows(Dataset):
    """One window: its bases and the binned, normalized signal of every track.
    Each worker keeps ``cache_bytes`` of decoded blocks.

    ``sequence``: uint8 ASCII ``[L]`` (for a tokenizer); ``signal``: float32
    ``[n_tracks, L / resolution]``. Reverse complement reads the window on
    the minus strand: strand-free tracks only reverse.
    """

    def __init__(
        self,
        reference: BgzipFasta,
        tracks: Sequence[BigWig],
        *,
        resolution: int = 128,
        reverse_complement: bool = True,
        cache_bytes: int = 256 << 20,
    ):
        self.reference, self.tracks = reference, list(tracks)
        self.resolution, self.rc = resolution, reverse_complement
        self.cache_bytes = cache_bytes
        self._readers: tuple[SequenceReader, list[TrackReader], np.ndarray] | None = None

    def _open(self) -> tuple[SequenceReader, list[TrackReader], np.ndarray]:
        # Opened in the worker on its first sample; readers are kept for every
        # later query. Keep dataset construction free of open readers: the
        # shared cache belongs to this process and cannot cross fork/pickle.
        if self._readers is None:
            cache = BlockCache(max_bytes=self.cache_bytes)
            sequence = SequenceReader(self.reference, cache=cache)
            signals = [TrackReader(track, cache=cache) for track in self.tracks]
            # A header summary is optional; empty/all-zero tracks need no scaling.
            means = np.ones(len(signals), np.float32)
            for i, reader in enumerate(signals):
                summary = reader.summary
                if summary is not None and summary.bases_covered:
                    mean = summary.sum / summary.bases_covered
                    if np.isfinite(mean) and mean > 0:
                        means[i] = mean
            self._readers = sequence, signals, means
        return self._readers

    def __getitem__(self, region: Region) -> dict[str, np.ndarray]:
        sequence, signals, means = self._open()
        if self.rc and np.random.random() < 0.5:
            # The minus strand: bases come reverse complemented, bins reversed.
            region = region.with_strand("-")
        bases = sequence.fetch(region).bases
        rows = []
        for reader, mean in zip(signals, means, strict=True):
            bins = reader.bins(region, self.resolution, "sum")
            # edges are local to the oriented query, including a partial last bin.
            rows.append(np.nan_to_num(bins["sum"]) / (np.diff(bins.edges) * mean))
        signal = (
            np.stack(rows).astype(np.float32)
            if rows
            else np.empty((0, -(-len(region) // self.resolution)), np.float32)
        )
        return {"sequence": bases, "signal": signal}


# ----------------------------------------------------------------- peaks


class PeakWindows(Sampler[tuple[Region, int]]):
    """Positives centred on union peaks, and as many negatives drawn uniformly
    from the chromosomes where the union has no peak. Yields (region, label)."""

    def __init__(self, union: BigBed, chroms: Mapping[str, int], length: int, *, seed: int = 0):
        self.union = union
        self.chroms = {c: n for c, n in chroms.items() if n >= length}
        self.length, self.seed = length, seed
        self.epoch = 0

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def __iter__(self) -> Iterator[tuple[Region, int]]:
        rng = np.random.default_rng(self.seed + self.epoch)
        half = self.length // 2
        items = []
        with IntervalReader(self.union) as peaks:
            for chrom, size in self.chroms.items():
                hits = peaks.fetch((chrom, 0, size))
                for centre in (hits.start.astype(np.int64) + hits.end) // 2:
                    s = int(np.clip(centre - half, 0, size - self.length))
                    items.append((Region(chrom, s, s + self.length), 1))
            positives = len(items)
            names, sizes = list(self.chroms), np.array(list(self.chroms.values()))
            while len(items) < 2 * positives:
                chrom = names[rng.choice(len(names), p=sizes / sizes.sum())]
                s = int(rng.integers(0, self.chroms[chrom] - self.length + 1))
                region = Region(chrom, s, s + self.length)
                if len(peaks.fetch(region)) == 0:
                    items.append((region, 0))
        rng.shuffle(items)
        return iter(items)


class PeakDataset(Dataset):
    """One window: ``sequence`` one-hot ``[L, 4]``, its ``label``, and
    ``coverage`` ``[n_tracks, L / resolution]``: which bins each track's
    peaks cover."""

    def __init__(self, reference: BgzipFasta, peaks: Sequence[BigBed], *, resolution: int = 100):
        self.reference, self.peaks, self.resolution = reference, list(peaks), resolution
        self._readers: tuple[SequenceReader, list[IntervalReader]] | None = None

    def _open(self) -> tuple[SequenceReader, list[IntervalReader]]:
        if self._readers is None:
            self._readers = (
                SequenceReader(self.reference),
                [IntervalReader(p) for p in self.peaks],
            )
        return self._readers

    def __getitem__(self, item: tuple[Region, int]) -> dict[str, np.ndarray]:
        region, label = item
        sequence, peak_sets = self._open()
        rows = [
            bin_coverage(h.start, h.end, region, self.resolution)
            for h in (p.fetch(region) for p in peak_sets)
        ]
        coverage = (
            np.stack(rows)
            if rows
            else np.empty((0, -(-len(region) // self.resolution)), np.float32)
        )
        return {
            "sequence": one_hot(sequence.fetch(region).bases),
            "coverage": coverage,
            "label": np.array(label, np.float32),
        }


# ------------------------------------------------------------------ main


def train_step(*tensors: torch.Tensor) -> None:
    """The model's forward and backward pass (not shown)."""


def main() -> None:
    reference = BgzipFasta("/data/hg38/hg38.fa.gz")
    tracks = [
        BigWig(f"/data/nonb/tracks/{name}.bw") for name in ("G4_K562_rep1", "Rloop_HEK293_rep1")
    ]
    peaks = [BigBed(f"/data/nonb/peaks/{name}.bb") for name in ("G4_K562", "Rloop_HEK293")]
    with SequenceReader(reference) as fa:
        chroms = {c: n for c, n in fa.chroms.items() if c in SPLITS["train"]}

    windows = TiledWindows(chroms, 131_072, seed=42)
    tracks_loader = DataLoader(
        TrackWindows(reference, tracks, resolution=128),
        sampler=windows,
        batch_size=8,
        num_workers=16,
        prefetch_factor=4,  # batches each worker keeps ready
        persistent_workers=True,  # readers stay open across epochs
        multiprocessing_context="spawn",  # native readers/caches are created in each worker
        pin_memory=True,
    )
    for epoch in range(3):
        windows.set_epoch(epoch)  # main process: the workers just receive regions
        for batch in tracks_loader:
            # sequence [B, L] uint8, signal [B, T, L / 128] float32
            train_step(batch["sequence"], batch["signal"])

    candidates = PeakWindows(BigBed("/data/nonb/peaks/union_peaks.bb"), chroms, 5_000, seed=42)
    peaks_loader = DataLoader(
        PeakDataset(reference, peaks, resolution=100),
        sampler=candidates,
        batch_size=256,
        num_workers=8,
        prefetch_factor=128,
        persistent_workers=True,
        multiprocessing_context="spawn",
    )
    for epoch in range(3):
        candidates.set_epoch(epoch)
        for batch in peaks_loader:
            # sequence [B, 5000, 4], label [B], coverage [B, T, 50]
            train_step(batch["sequence"], batch["label"], batch["coverage"])


if __name__ == "__main__":
    main()
