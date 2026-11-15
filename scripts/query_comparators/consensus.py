"""Real Python/CLI consensus comparisons; imported by the central runner.

Window FASTA creation happens in prepare/build_cases, outside all factories
and timing. OmicsWeave factories load the whole BCF store and open the full
reference. pysam and CLI have no persistent consensus reader: their factories
are no-ops and each real query opens its inputs again. Report open and query
costs separately, and disclose the precomputed window inputs.
"""

from __future__ import annotations

import hashlib
import os
import random
import struct
import subprocess
from pathlib import Path

WIDTHS = (128, 4096, 131072)
SEED = 0xC0A5E85
CHROMS = ("chr1", "chr2")
SELECTOR_FLAGS = {
    "all": ("-s", "-"),
    "hap1": ("-s", "S", "-H", "1"),
    "hap2": ("-s", "S", "-H", "2"),
}


def _window_path(fixtures: Path, q: tuple[str, int, int]) -> Path:
    key = f"{q[0]}:{q[1]}-{q[2]}".encode("ascii")
    return fixtures / "consensus-windows" / (hashlib.sha256(key).hexdigest()[:24] + ".fa")


def prepare(fixtures: Path) -> dict[int, tuple[tuple[str, int, int], ...]]:
    """Prepare twelve deterministic windows from existing indexed fixtures.

    Required files: ref.fa/ref.fa.fai and variants.bcf/variants.bcf.csi.
    Reference contigs chr1 and chr2 must each contain 2**22 uppercase ACGTN
    bases. BCF records must be sorted biallelic A>C SNVs at REF=A positions,
    with sample S, GT=0|1, approximately one record per 101 reference bases.
    The central fixture generator owns those source files; this helper only
    writes window FASTAs, verifies their reference bytes, and returns regions.
    """
    import pysam

    fixtures = Path(fixtures)
    for name in ("ref.fa", "ref.fa.fai", "variants.bcf", "variants.bcf.csi"):
        if not (fixtures / name).is_file():
            raise FileNotFoundError(fixtures / name)
    windows = fixtures / "consensus-windows"
    windows.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    regions: dict[int, tuple[tuple[str, int, int], ...]] = {}
    with pysam.FastaFile(os.fspath(fixtures / "ref.fa")) as reference:
        for chrom in CHROMS:
            if reference.get_reference_length(chrom) != 2**22:
                raise ValueError(f"unexpected reference length: {chrom}")
        for width in WIDTHS:
            queries = []
            for i in range(4):
                chrom = CHROMS[i % len(CHROMS)]
                last = reference.get_reference_length(chrom) - width
                start = rng.randrange(last + 1)
                # N is a valid reference base; prefer a called first base so
                # all backends also exercise an ordinary window boundary.
                for _ in range(32):
                    if reference.fetch(chrom, start, start + 1) != "N":
                        break
                    start = (start + 1) % (last + 1)
                q = (chrom, start, start + width)
                bases = reference.fetch(*q).encode("ascii")
                if len(bases) != width or set(bases) - set(b"ACGTN"):
                    raise ValueError(f"invalid reference bytes for {q!r}")
                header = f">{chrom}:{start + 1}-{start + width}\n".encode("ascii")
                body = b"\n".join(bases[j : j + 60] for j in range(0, len(bases), 60)) + b"\n"
                expected = header + body
                path = _window_path(fixtures, q)
                if not path.exists() or path.read_bytes() != expected:
                    path.write_bytes(expected)
                queries.append(q)
            regions[width] = tuple(queries)
    return regions


def _fasta_sequence(raw, q) -> bytes:
    data = raw.encode("ascii") if isinstance(raw, str) else bytes(raw)
    lines = data.splitlines()
    if sum(line.startswith(b">") for line in lines) != 1:
        raise ValueError(f"expected exactly one consensus sequence for {q!r}")
    bases = b"".join(line for line in lines if line and not line.startswith(b">"))
    if len(bases) != q[2] - q[1] or set(bases) - set(b"ACGTN"):
        raise ValueError(f"invalid SNP-only consensus output for {q!r}")
    return bases


def _pack(sequences) -> bytes:
    # Prefix the count and every sequence length to preserve batch boundaries.
    values = tuple(sequences)
    return struct.pack("<Q", len(values)) + b"".join(
        struct.pack("<Q", len(value)) + value for value in values
    )


def _normalize_omics_single(raw, q) -> bytes:
    if raw.chain is not None:
        raise ValueError("chain must be disabled for this comparison")
    return raw.bases.tobytes()


def _normalize_omics_batch(raw, q) -> bytes:
    if len(raw.alleles) != 3 or len(raw.offsets) != 4:
        raise ValueError("expected three selectors in one region")
    if any(chain is not None for chain in raw.chains):
        raise ValueError("chain must be disabled for this comparison")
    offsets = tuple(int(offset) for offset in raw.offsets)
    if offsets[0] != 0 or offsets[-1] != len(raw.bases):
        raise ValueError("invalid packed consensus offsets")
    return _pack(raw.bases[a:b].tobytes() for a, b in zip(offsets, offsets[1:], strict=False))


def _normalize_external_batch(raw, q) -> bytes:
    if len(raw) != 3:
        raise ValueError("expected three real consensus calls")
    return _pack(_fasta_sequence(item, q) for item in raw)


def build_cases(ow, fixtures: Path, tools: dict) -> list[dict]:
    import pysam.bcftools

    fixtures = Path(fixtures)
    regions = prepare(fixtures)
    bcf = os.fspath(fixtures / "variants.bcf")
    reference = os.fspath(fixtures / "ref.fa")
    bcftools = os.fspath(tools["bcftools"])
    consensus_dispatch = pysam.bcftools.consensus
    arguments = {
        (q, selected): (*flags, "-f", os.fspath(_window_path(fixtures, q)), bcf)
        for queries in regions.values()
        for q in queries
        for selected, flags in SELECTOR_FLAGS.items()
    }

    def omics_factory(mode):
        def factory():
            reader = ow.query.ConsensusReader(bcf, reference)
            selectors = (
                ow.types.Haplotype("S", "1"),
                ow.types.Haplotype("S", "2"),
                ow.types.AllAlt(),
            )
            if mode == "batch":

                def query(q):
                    return reader.fetch_batch(q, selectors)

                normalize = _normalize_omics_batch
            elif mode == "hap2":

                def query(q):
                    return reader.fetch(q, selectors[1])

                normalize = _normalize_omics_single
            else:

                def query(q):
                    return reader.fetch(q)

                normalize = _normalize_omics_single
            return {"query": query, "normalize": normalize, "close": reader.close}

        return factory

    def external_factory(kind, mode):
        def factory():
            modes = ("hap1", "hap2", "all") if mode == "batch" else (mode,)

            def one(q, selected):
                args = arguments[q, selected]
                if kind == "python":
                    return consensus_dispatch(*args, catch_stdout=True, raw=True)
                return subprocess.run(
                    [bcftools, "consensus", *args],
                    capture_output=True,
                    timeout=30,
                    check=True,
                ).stdout

            def query(q):
                if mode == "batch":
                    return tuple(one(q, selected) for selected in modes)
                return one(q, modes[0])

            return {
                "query": query,
                "normalize": _normalize_external_batch if mode == "batch" else _fasta_sequence,
                "close": lambda: None,
            }

        return factory

    cases = []
    for width in WIDTHS:
        for mode, interface in (("all", "fetch"), ("hap2", "fetch"), ("batch", "fetch_batch")):
            backends = [
                {
                    "name": "omicsweave",
                    "kind": "omics",
                    "factory": omics_factory(mode),
                    "native_batch": mode == "batch",
                    "note": "Factory loads the full BCF store and opens the full reference; chain=False.",
                },
                {
                    "name": "pysam.bcftools.consensus",
                    "kind": "python",
                    "factory": external_factory("python", mode),
                    "native_batch": False,
                    "note": "No persistent reader: no-op factory; every real consensus call opens inputs and returns raw FASTA bytes.",
                },
                {
                    "name": "bcftools consensus",
                    "kind": "cli",
                    "factory": external_factory("cli", mode),
                    "native_batch": False,
                    "note": "No-op factory; every query includes a real subprocess, file opening, and captured FASTA stdout.",
                },
            ]
            cases.append(
                {
                    "id": f"consensus.{interface}.{mode}.{width}",
                    "interface": f"ConsensusReader.{interface}",
                    "format": "BCF+FASTA",
                    "workload": f"SNV-only, width={width}, selector={mode}, chain=False",
                    "queries": regions[width],
                    "backends": backends,
                    "note": (
                        "Four distinct fixed-seed windows per width, prepared outside timing; "
                        "window FASTA inputs benefit bcftools and must be disclosed. "
                        "AllAlt uses fetch's default; Hap2 uses S/2. "
                        "Batch is one region with [S/1,S/2,AllAlt], "
                        "compared with three successive real external consensus calls. "
                        "No chain filesystem output, fixture preparation, FASTA normalization, "
                        "or checksum calculation belongs in raw query timing."
                    ),
                }
            )
    return cases
