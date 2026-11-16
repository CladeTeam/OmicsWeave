from __future__ import annotations

import math

import numpy as np
import pytest

from omicsweave import motif, query, types

from ._support import DNA, random_sequences
from ._support import collection as collection


def test_statistics_steps():
    p = np.array([0.04, 0.01, 0.03, 0.2])
    m = len(p)
    order = np.argsort(p)
    bh = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    want = np.empty(m)
    want[order] = np.minimum(bh, 1)
    np.testing.assert_allclose(motif.stats.BenjaminiHochberg()(p), want)
    np.testing.assert_allclose(motif.stats.EValue(tests=10)(p), p * 10)
    np.testing.assert_allclose(motif.stats.Sidak()(p, 5), 1 - (1 - p) ** 5)
    x = np.prod(p)
    assert np.isclose(
        motif.stats.Product()(p), x * sum((-np.log(x)) ** i / math.factorial(i) for i in range(m))
    )
    assert motif.stats.MinBonferroni()(p) == pytest.approx(0.04)
    assert motif.stats.EmpiricalTail()(np.arange(10.0))(7.5) == pytest.approx(3 / 11)


def test_mast_with_python_steps(collection):
    _, motifs = collection
    pwms = motifs[:5]
    rng = np.random.default_rng(2)
    seqs = random_sequences(rng, 10, 200)
    bg = motif.Background.uniform(DNA)
    native = motif.mast(pwms, seqs, background=bg, alphabet=DNA)
    calls = []

    def product(pvalues: np.ndarray) -> float:
        calls.append(len(pvalues))
        return motif.stats.Product()(pvalues)

    python = motif.mast(
        pwms, seqs, background=bg, combine=product, adjust=lambda p: p * len(p), alphabet=DNA
    )
    assert calls
    np.testing.assert_allclose(python.pvalue, native.pvalue)
    np.testing.assert_allclose(python.evalue, native.evalue)
    np.testing.assert_array_equal(python.ranked, native.ranked)
    assert np.all(np.diff(native.hits.indptr) >= 0) and native.best.indptr[-1] == len(
        native.best.motif
    )
    text = motif.text.mast_hit_list(
        native, motif_ids=[m.id for m in pwms], names=[f"s{i}" for i in range(10)]
    )
    assert len(text.splitlines()) == len(native.hits.motif)

    def broken(_: np.ndarray) -> float:
        raise RuntimeError("from Python")

    with pytest.raises(RuntimeError, match="from Python"):
        motif.mast(pwms, seqs, background=bg, combine=broken, alphabet=DNA)


def test_mcast_clusters(collection):
    _, motifs = collection
    pwms = motifs[:5]
    rng = np.random.default_rng(9)
    site = "".join("ACGT"[i] for i in motifs[0].matrix.argmax(1))
    seq = (
        "".join(rng.choice(list("ACGT"), 300))
        + (site + "TT") * 6
        + "".join(rng.choice(list("ACGT"), 300))
    )
    bg = motif.Background.uniform(DNA)
    null = motif.stats.BackgroundNull(bg, length=50000, sequences=20, seed=3)
    stats = motif.stats.ClusterStatistics(
        null, motif.stats.Unconditional(), qvalue=motif.stats.BenjaminiHochberg()
    )
    table = motif.mcast(pwms, [seq], background=bg, statistics=stats, alphabet=DNA)
    assert len(table) >= 1 and table.statistics
    assert table.start[0] <= 300 < table.end[0]
    assert np.all((table.pvalue > 0) & (table.pvalue <= 1))
    plain = motif.mcast(pwms, [seq], background=bg, statistics=False, alphabet=DNA)
    assert not plain.statistics and np.isnan(plain.pvalue).all()
    text = motif.text.mcast(table, names=["chr"], sequences=[seq], evalue=float("inf"))
    lines = text.splitlines()
    assert lines[0].startswith("pattern_name\t") and len(lines) == len(table) + 1


def test_glam2scan_profiles():
    profile = types.Motif(
        "P",
        4,
        tuple(types.Factor("joint", (i,), np.eye(4)[c] * 9) for i, c in enumerate([0, 1, 2, 3])),
        sites=9.0,
        gaps=types.Gaps(np.zeros(4), np.zeros(3)),
    )
    seqs = ["TTTTACGTTTTT", "GGGGGGG"]
    hits = motif.glam2scan([profile], seqs, alphabet=DNA, top=3)
    assert len(hits) == 3 and hits.sequence[0] == 0 and hits.start[0] == 4
    assert np.all(np.diff(hits.score) <= 0)
    assert hits.column_indptr[-1] == len(hits.matched) == 4 * len(hits)
    text = motif.text.glam2scan(hits, alphabet=DNA, names=["a", "b"], sequences=seqs)
    assert "acgt" in text


def test_fimo_text(collection):
    path, motifs = collection
    site = "".join("ACGT"[i] for i in motifs[2].matrix.argmax(1))
    seqs = ["GG" + site + "GG"]
    hits = motif.fimo(path, seqs, background=motif.Background.uniform(DNA), pvalue=1e-3)
    assert 2 in hits.start[hits.motif == 2]
    text = motif.text.fimo(
        hits, alphabet=DNA, motif_ids=query.MotifReader(path).ids, names=["s"], sequences=seqs
    )
    lines = text.splitlines()
    assert lines[0].startswith("motif_id\tmotif_alt_id") and len(lines) == len(hits) + 1
