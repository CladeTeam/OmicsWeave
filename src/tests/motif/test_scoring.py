from __future__ import annotations

import numpy as np
import pytest

from omicsweave import motif, query, types

from ._support import DNA, random_sequences
from ._support import collection as collection


def window_scores(cm: motif.CompiledMotif, seq: str) -> dict[int, int]:
    """Each scorable window start's integer score, summed from the compiled
    integer terms (forward strand, 0-order background: no context terms)."""
    codes = DNA.lut[np.frombuffer(seq.encode(), dtype=np.uint8)]
    out = {}
    for start in range(len(seq) - cm.width + 1):
        window = codes[start : start + cm.width]
        if (window == 255).any():
            continue
        out[start] = sum(int(t.table[tuple(window[list(t.positions)])]) for t in cm.int_terms)
    return out


@pytest.mark.parametrize("rounding", [None, 2, 6])
def test_cutoffs_match_reported_survival(rounding):
    survival = np.array([1.0, 0.8131, 0.8129, 0.5, 0.1, 0.0])
    cm = motif.CompiledMotif(
        0,
        "+",
        1,
        (),
        (),
        1.0,
        0.0,
        len(survival) - 1,
        motif.NullDistribution(np.zeros_like(survival), survival, 0.0),
        rounding,
    )
    reported = cm.pvalue(np.arange(cm.max + 1))
    for p in [*survival, 0.81, -1.0, 2.0, float("nan")]:
        hits = np.flatnonzero(reported <= p)
        assert cm.cutoff(p) == (int(hits[0]) if len(hits) else None)


def test_empty_program_score_shapes_keep_sequence_lengths():
    program = motif.compile([], motif.Background.uniform(DNA), alphabet=DNA)
    scores = program.scan(["ACG", ""], motif.Scores())
    assert len(scores) == 2
    assert scores[0].shape == (0, scores.strands, 3)
    assert scores[1].shape == (0, scores.strands, 0)


def test_forbidden_windows_have_no_finite_score_probability():
    cm = motif.CompiledMotif(
        0,
        "+",
        1,
        (),
        (),
        1.0,
        0.0,
        1,
        motif.NullDistribution(np.array([0.3, 0.2]), np.array([0.5, 0.2]), 0.5),
        None,
    )
    np.testing.assert_array_equal(cm.pvalue([-1, 0, 1, 2]), [0.5, 0.5, 0.2, 0.0])


def test_scan_matches_compiled_models(collection):
    path, _ = collection
    rng = np.random.default_rng(3)
    seqs = random_sequences(rng, 6, 300)
    bg = motif.Background.frequencies([0.3, 0.2, 0.2, 0.3])
    reader = query.MotifReader(path)
    for scoring in (motif.Scoring.fimo(), motif.Scoring.exact()):
        program = motif.compile(reader, bg, scoring, strands="+")
        hits = program.scan(seqs, motif.Hits(pvalue=1e-2))
        found = set(
            zip(hits.motif.tolist(), hits.sequence.tolist(), hits.start.tolist(), strict=True)
        )
        want = set()
        for m in range(len(program)):
            cm = program.model(m)
            for s, seq in enumerate(seqs):
                for start, score in window_scores(cm, seq).items():
                    if cm.pvalue(score) <= 1e-2:
                        want.add((m, s, start))
        assert found == want
        cm = program.model(0)
        assert cm.cutoff(1e-2) is not None
        assert cm.pvalue(cm.cutoff(1e-2)) <= 1e-2 < cm.pvalue(cm.cutoff(1e-2) - 1)
        assert np.isclose(cm.null.pmf.sum() + cm.null.forbidden, 1.0)


def test_reductions_share_one_scan(collection):
    path, _ = collection
    rng = np.random.default_rng(11)
    seqs = random_sequences(rng, 8, 400)
    program = motif.compile(path, motif.Background.uniform(DNA), motif.Scoring.exact())
    hits, best, counts, scores = program.scan(
        seqs,
        motif.Hits(pvalue=1e-3),
        motif.Best(),
        motif.Counts(pvalue=1e-3),
        motif.Scores(),
        threads=4,
        tile_len=97,
    )
    alone = program.scan(seqs, motif.Hits(pvalue=1e-3), threads=1)
    for col in ("sequence", "motif", "start", "strand", "score", "pvalue"):
        np.testing.assert_array_equal(getattr(hits, col), getattr(alone, col))
    dense = np.zeros(counts.shape, dtype=np.uint32)
    np.add.at(dense, (hits.sequence, hits.motif), 1)
    np.testing.assert_array_equal(counts.to_dense(), dense)
    for s, m, start, strand, score in zip(
        hits.sequence, hits.motif, hits.start, hits.strand, hits.score, strict=True
    ):
        assert scores[s][m, 0 if strand == 1 else 1, start] == np.float32(score)
    for s, m, score in zip(best.sequence, best.motif, best.score, strict=True):
        assert np.float32(score) == np.nanmax(scores[s][m])


def test_combined_parts(collection):
    path, _ = collection
    seqs = random_sequences(np.random.default_rng(5), 4, 300)
    bg = motif.Background.uniform(DNA)
    fimo = motif.compile(path, bg, motif.Scoring.fimo())
    exact = motif.compile(path, bg, motif.Scoring.exact(), strands="+")
    both = motif.combine(fimo, exact)
    assert len(both) == 2 * len(fimo) and [p.motifs for p in both.parts] == [
        range(0, 6),
        range(6, 12),
    ]
    a, b = both.scan(seqs, motif.Hits(pvalue=1e-3, part=0), motif.Hits(pvalue=1e-3, part=1))
    np.testing.assert_array_equal(a.score, fimo.scan(seqs, motif.Hits(pvalue=1e-3)).score)
    np.testing.assert_array_equal(b.motif - 6, exact.scan(seqs, motif.Hits(pvalue=1e-3)).motif)
    with pytest.raises(IndexError):
        both.scan(seqs, motif.Hits(part=2))


def test_compile_needs_an_alphabet_for_values():
    m = types.Motif.from_pfm("A", np.eye(4))
    with pytest.raises(TypeError):
        motif.compile([m], motif.Background.uniform(DNA))
    with pytest.raises(TypeError):
        motif.compile([m], motif.Background.uniform(DNA), alphabet=DNA).scan("ACGT")
