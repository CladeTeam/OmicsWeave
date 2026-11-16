import json

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import quant
from omicsweave.annotation import Annotation
from omicsweave.errors import OmicsWeaveError
from omicsweave.query import TrackReader

from ._support import CHROMS, EXONS, read_10x
from ._support import data as data


def test_coverage_and_single_cell_genes(data):
    tmp, recs, bam, ann = data
    out = tmp / "out"
    task = quant.Task(
        [bam], out, [quant.Coverage(), quant.Counts.rna_single_cell()], annotation=ann, matrix="mtx"
    )
    p = quant.plan(task, threads=2, memory="256M")
    assert [s[0] for s in p.stages] == ["coverage", "genes"]
    assert p.memory.total == p.budget.memory == 256 << 20
    assert p.budget.threads == 2
    assert [i.path for i in p.inputs] == [bam]
    report = p.run()
    assert report.records == len(recs)
    kinds = {k for k, _ in report.files}
    assert {"bigwig", "mtx", "barcodes", "features", "json"} <= kinds
    assert json.loads((out / "summary.json").read_text()) == report.summary

    depth = {c: np.zeros(n, np.float32) for c, n in CHROMS}
    for chrom, pos, *_ in recs:
        depth[chrom][pos : pos + 50] += 1
    bw = TrackReader(out / "coverage.bw")
    for c, n in CHROMS:
        got = np.nan_to_num(bw.signal((c, 0, n)).values)
        np.testing.assert_array_equal(got, depth[c])

    want = {}
    seen = set()
    for chrom, pos, flag, _, _, cb, ub, _ in recs:
        # Forward-stranded: the fragment strand is read 1's.
        rev = bool(flag & 0x10) ^ bool(flag & 0x80)
        hits = [
            g
            for g, (c, ex, strand) in EXONS.items()
            if c == chrom
            and (strand == "-") == rev
            and any(a < pos + 50 and pos < b for a, b in ex)
        ]
        if len(hits) == 1 and (cb, hits[0], ub) not in seen:
            seen.add((cb, hits[0], ub))
            want[(cb, hits[0])] = want.get((cb, hits[0]), 0) + 1
    assert read_10x(out / "genes") == want


def test_regions_count_as_their_own_features(data):
    tmp, _, bam, _ = data
    size = 2_000
    chroms = [c for c, _ in CHROMS]
    seqid, start, end = [], [], []
    for i, (_, n) in enumerate(CHROMS):
        for s in range(0, n, size):
            seqid.append(i)
            start.append(s)
            end.append(min(s + size, n))
    tiles = ow.annotation.IntervalSet.from_arrays(
        chroms, np.array(seqid), np.array(start), np.array(end)
    )
    task = quant.Task(
        [bam],
        tmp / "regions",
        [
            quant.Counts.chip(quant.Bins(size), name="bins"),
            quant.Counts.chip(quant.Regions(tiles), name="regions"),
        ],
        matrix="mtx",
    )
    quant.run(task, threads=2, memory="256M")
    bins = read_10x(tmp / "regions" / "bins")
    assert bins and read_10x(tmp / "regions" / "regions") == bins


def test_a_subset_of_the_annotation_picks_the_genes(data):
    tmp, _, bam, ann = data
    only_a = ann.subset(genes=["A"])
    task = quant.Task(
        [bam],
        tmp / "subset",
        [
            quant.Counts.rna_single_cell(name="all"),
            quant.Counts.rna_single_cell(name="a", subset=only_a, velocity=True),
        ],
        annotation=ann,
        matrix="mtx",
    )
    quant.run(task, threads=2, memory="256M")
    every, a = read_10x(tmp / "subset" / "all"), read_10x(tmp / "subset" / "a")
    assert a and {f for _, f in a} == {"A"}
    assert a == {k: v for k, v in every.items() if k[1] == "A"}
    other = Annotation.read(tmp / "genes.gtf").subset(genes=["C"])
    rows = type(other)(other.genes + 10, other.transcripts)
    with pytest.raises(OmicsWeaveError, match="rows of the annotation"):
        quant.plan(
            quant.Task([bam], tmp / "x", [quant.Counts.rna_bulk(subset=rows)], annotation=ann)
        )
