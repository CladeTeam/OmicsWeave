"""Cached metadata methods, separate from payload query throughput."""

import argparse
import json
import statistics
import time
from pathlib import Path


def measure(fn):
    expected = fn()
    values = []
    for _ in range(5):
        count = 0
        start = time.perf_counter_ns()
        while count == 0 or time.perf_counter_ns() - start < 100_000_000:
            got = fn()
            count += 1
        assert got == expected
        values.append((time.perf_counter_ns() - start) / count)
    return {"samples_ns": values, "median_ns": statistics.median(values), "value": expected}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    import hicstraw
    import py2bit
    import pysam

    import omicsweave as ow

    f = a.fixtures
    rows = []
    with ow.query.SequenceReader(f / "ref.fa") as ours, pysam.FastaFile(str(f / "ref.fa")) as other:
        for name, fn in [
            ("OmicsWeave.length", lambda: ours.length("chr1")),
            ("pysam.get_reference_length", lambda other=other: other.get_reference_length("chr1")),
        ]:
            rows.append({"interface": "SequenceReader.length", "backend": name, **measure(fn)})
    with ow.query.SequenceReader(f / "ref.2bit") as ours:
        other = py2bit.open(str(f / "ref.2bit"))
        try:
            for name, fn in [
                ("OmicsWeave.length", lambda: ours.length("chr1")),
                ("py2bit.chroms", lambda other=other: other.chroms("chr1")),
            ]:
                rows.append(
                    {"interface": "SequenceReader.length.2bit", "backend": name, **measure(fn)}
                )
        finally:
            other.close()
    with ow.query.ContactReader(f / "round8-bp-v9.hic") as ours:
        other = hicstraw.HiCFile(str(f / "round8-bp-v9.hic"))
        for name, fn in [
            ("OmicsWeave.resolutions", lambda: tuple(ours.resolutions("BP"))),
            ("hicstraw.getResolutions", lambda other=other: tuple(other.getResolutions())),
        ]:
            rows.append({"interface": "ContactReader.resolutions", "backend": name, **measure(fn)})
        del other
    with ow.query.ContactReader(f / "round8-frag-v9.hic") as ours:
        rows.append(
            {
                "interface": "ContactReader.fragment_count",
                "backend": "OmicsWeave",
                **measure(lambda: ours.fragment_count("chr1")),
                "baseline": "N/A: no equivalent public fragment-count method in selected hic-straw/CLI",
            }
        )
    a.output.write_text(
        json.dumps(
            {
                "rows": rows,
                "scope": "warm cached metadata; no payload decoding; open costs excluded; five 100ms samples",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
