"""Real Python/CLI query adapters; the coordinator owns measurement.

Fixture contract (all backends consume the same files, never rewritten here):
* signal.bw: chr1/chr2, each 2**22 bases; finite float32 values on sorted,
  positive-width, non-overlapping intervals, with gaps. Use exactly
  representable small values (for example [-16, 16] in steps of 1/8) so UCSC
  bedGraph text round-trips to float32 even with its text precision limit.
* features.bed and features.bb: identical records, stable chrom/start order;
  positive widths; ASCII annotation bytes after the third column, including
  literal TAB separators and empty fields. Keep a unique record name as the
  first extra field. Endpoints remain inside each chromosome.
* features.bed.gz: BGZF encoding of exactly features.bed, with .tbi index.

The recommended fixture contains at least tens of thousands of features per
chromosome, nested overlaps, equal starts, varied annotation lengths, and
gaps. No zero-width features: bedtools has special point-feature semantics.

factory includes open/parse/index building; query eagerly materializes lazy
results; normalize returns deterministic column bytes. The coordinator times
raw queries and canonical queries separately and verifies normalize outside
raw timing. Normalize each raw result before the next query: the omics out
backend deliberately reuses its caller-owned array. CLI query cost includes
process startup, I/O and stdout materialization. NCLS factory includes BED
parsing and construction of every per-chromosome index.
"""

from __future__ import annotations

import random
import struct
import subprocess
from pathlib import Path

import numpy as np

CHROM_LENGTH = 1 << 22
WIDTHS = (128, 4096, 131072)


def _queries(width):
    rng = random.Random(0x8BB1 + width)
    return [
        (f"chr{1 + i % 2}", start, start + width)
        for i in range(16)
        for start in [rng.randrange(CHROM_LENGTH - width + 1)]
    ]


def _floats(values):
    out = np.array(values, dtype="<f4", copy=True)
    # Every producer uses one quiet-NaN bit pattern, including missing bases.
    out.view("<u4")[np.isnan(out)] = 0x7FC00000
    return out.tobytes()


def _signal(raw, query):
    values = raw.values if hasattr(raw, "values") else raw
    if len(values) != query[2] - query[1]:
        raise ValueError("dense signal length differs from query width")
    return _floats(values)


def _track_rows(raw):
    if isinstance(raw, bytes):
        return [
            (int(f[1]), int(f[2]), float(f[3]))
            for line in raw.splitlines()
            if (f := line.split(b"\t"))
        ]
    if hasattr(raw, "start"):
        return list(zip(raw.start, raw.end, raw.value, strict=True))
    return list(raw or ())


def _track(raw, query):
    _, lo, hi = query
    if hasattr(raw, "start"):
        # The public result already guarantees clipped, stable-start columns.
        return (
            struct.pack("<Q", len(raw.start))
            + np.asarray(raw.start, dtype="<u4").tobytes()
            + np.asarray(raw.end, dtype="<u4").tobytes()
            + _floats(raw.value)
        )
    rows = [(max(int(s), lo), min(int(e), hi), v) for s, e, v in _track_rows(raw)]
    rows = sorted((r for r in rows if r[0] < r[1]), key=lambda row: row[0])
    return (
        struct.pack("<Q", len(rows))
        + np.asarray([r[0] for r in rows], dtype="<u4").tobytes()
        + np.asarray([r[1] for r in rows], dtype="<u4").tobytes()
        + _floats([r[2] for r in rows])
    )


def _bedgraph_signal(raw, query):
    _, lo, hi = query
    values = np.full(hi - lo, np.nan, dtype="<f4")
    for start, end, value in _track_rows(raw):
        start, end = max(start, lo), min(end, hi)
        if start < end:
            values[start - lo : end - lo] = value
    return _floats(values)


def _bed_rows(raw):
    if hasattr(raw, "annotation"):
        return [
            (int(start), int(end), raw.annotation.bytes(i))
            for i, (start, end) in enumerate(zip(raw.start, raw.end, strict=True))
        ]
    if isinstance(raw, bytes):
        raw = raw.splitlines()
    rows = []
    for record in raw or ():
        if isinstance(record, (bytes, str)):
            if isinstance(record, str):
                record = record.encode("ascii")
            fields = record.split(b"\t", 3)
            rows.append((int(fields[1]), int(fields[2]), fields[3] if len(fields) == 4 else b""))
        else:
            start, end, note = record
            rows.append(
                (
                    int(start),
                    int(end),
                    note.encode("ascii") if isinstance(note, str) else bytes(note),
                )
            )
    return rows


def _features(raw, query):
    # Full feature coordinates and every annotation byte are retained.
    if hasattr(raw, "annotation"):
        # Preserve native columns directly; no intermediate Python row table.
        return (
            struct.pack("<Q", len(raw.start))
            + np.asarray(raw.start, dtype="<u4").tobytes()
            + np.asarray(raw.end, dtype="<u4").tobytes()
            + np.asarray(raw.annotation.offsets, dtype="<u8").tobytes()
            + raw.annotation.data.tobytes()
        )
    rows = sorted(_bed_rows(raw), key=lambda row: row[0])
    offsets, data = [0], bytearray()
    for _, _, note in rows:
        data.extend(note)
        offsets.append(len(data))
    return (
        struct.pack("<Q", len(rows))
        + np.asarray([r[0] for r in rows], dtype="<u4").tobytes()
        + np.asarray([r[1] for r in rows], dtype="<u4").tobytes()
        + np.asarray(offsets, dtype="<u8").tobytes()
        + bytes(data)
    )


def _backend(name, kind, factory):
    return {"name": name, "kind": kind, "factory": factory}


def _omics_track(ow, path, mode, width):
    def factory():
        reader = ow.query.TrackReader(path)
        if mode == "intervals":
            query, normalize = reader.intervals, _track
        elif mode == "out":
            out = np.empty(width, dtype=np.float32)

            def query(region):
                return reader.signal(region, out=out)

            normalize = _signal
        else:
            query, normalize = reader.signal, _signal
        return {"query": query, "normalize": normalize, "close": reader.close, "reader": reader}

    return factory


def _omics_features(ow, path):
    def factory():
        reader = ow.query.IntervalReader(path)
        return {
            "query": reader.fetch,
            "normalize": _features,
            "close": reader.close,
            "reader": reader,
        }

    return factory


def _pybigwig(path, mode):
    def factory():
        import pyBigWig

        reader = pyBigWig.open(str(path))
        if mode == "signal":

            def query(q):
                return reader.values(*q, numpy=True)

            normalize = _signal
        elif mode == "intervals":

            def query(q):
                return list(reader.intervals(*q) or ())

            normalize = _track
        else:

            def query(q):
                return list(reader.entries(*q) or ())

            normalize = _features
        return {"query": query, "normalize": normalize, "close": reader.close}

    return factory


def _ucsc(binary, path, normalize):
    def factory():
        def query(q):
            chrom, start, end = q
            return subprocess.run(
                [
                    str(binary),
                    f"-chrom={chrom}",
                    f"-start={start}",
                    f"-end={end}",
                    str(path),
                    "/dev/stdout",
                ],
                check=True,
                capture_output=True,
                timeout=30,
            ).stdout

        return {"query": query, "normalize": normalize, "close": lambda: None}

    return factory


def _bedtools(binary, path):
    def factory():
        def query(q):
            chrom, start, end = q
            return subprocess.run(
                [str(binary), "intersect", "-a", str(path), "-b", "stdin", "-wa"],
                input=f"{chrom}\t{start}\t{end}\n".encode("ascii"),
                check=True,
                capture_output=True,
                timeout=30,
            ).stdout

        return {"query": query, "normalize": _features, "close": lambda: None}

    return factory


def _tabix(path):
    def factory():
        import pysam

        reader = pysam.TabixFile(str(path))
        return {
            "query": lambda q: list(reader.fetch(*q)),
            "normalize": _features,
            "close": reader.close,
        }

    return factory


def _ncls(path):
    def factory():
        from ncls import NCLS

        groups = {}
        with path.open("rb") as source:
            for line in source:
                fields = line.rstrip(b"\n").split(b"\t", 3)
                chrom = fields[0].decode("ascii")
                groups.setdefault(chrom, []).append(
                    (int(fields[1]), int(fields[2]), fields[3] if len(fields) == 4 else b"")
                )
        indexes = {}
        for chrom, rows in groups.items():
            indexes[chrom] = NCLS(
                np.asarray([r[0] for r in rows], dtype=np.int64),
                np.asarray([r[1] for r in rows], dtype=np.int64),
                np.arange(len(rows), dtype=np.int64),
            )

        def query(q):
            chrom, start, end = q
            ids = sorted(i for _, _, i in indexes[chrom].find_overlap(start, end))
            return [groups[chrom][i] for i in ids]

        def close():
            indexes.clear()
            groups.clear()

        return {"query": query, "normalize": _features, "close": close}

    return factory


def build_cases(ow, fixtures: Path, tools: dict) -> list[dict]:
    """Build cases without opening fixtures, installing tools, or executing queries."""
    fixtures = Path(fixtures)
    bw, bb = fixtures / "signal.bw", fixtures / "features.bb"
    bed, bgz = fixtures / "features.bed", fixtures / "features.bed.gz"
    cases = []
    for width in WIDTHS:
        queries = _queries(width)

        def add(interface, format_, backends, width=width, queries=queries):
            cases.append(
                {
                    "id": f"{interface}.{format_}.w{width}",
                    "interface": interface,
                    "format": format_,
                    "workload": f"16 fixed-seed legal queries; width={width}; chr1/chr2",
                    "queries": queries.copy(),
                    "backends": backends,
                }
            )

        add(
            "Track.signal",
            "bigwig",
            [
                _backend("omicsweave.alloc", "omics", _omics_track(ow, bw, "alloc", width)),
                _backend("omicsweave.out", "omics", _omics_track(ow, bw, "out", width)),
                _backend("pyBigWig.values.numpy", "python", _pybigwig(bw, "signal")),
                _backend(
                    "bigWigToBedGraph",
                    "cli",
                    _ucsc(tools["bigWigToBedGraph"], bw, _bedgraph_signal),
                ),
            ],
        )
        add(
            "Track.intervals",
            "bigwig",
            [
                _backend("omicsweave", "omics", _omics_track(ow, bw, "intervals", width)),
                _backend("pyBigWig.intervals", "python", _pybigwig(bw, "intervals")),
                _backend("bigWigToBedGraph", "cli", _ucsc(tools["bigWigToBedGraph"], bw, _track)),
            ],
        )
        add(
            "Interval.fetch",
            "bigbed",
            [
                _backend("omicsweave", "omics", _omics_features(ow, bb)),
                _backend("pyBigWig.entries", "python", _pybigwig(bb, "entries")),
                _backend("bigBedToBed", "cli", _ucsc(tools["bigBedToBed"], bb, _features)),
                # OmicsWeave reads BED only after converting it to bigBed; the
                # same records as text are the baselines of text-based tools.
                _backend("NCLS.inmemory", "python", _ncls(bed)),
                _backend("pysam.TabixFile", "python", _tabix(bgz)),
                _backend("bedtools.intersect", "cli", _bedtools(tools["bedtools"], bed)),
            ],
        )
    return cases
