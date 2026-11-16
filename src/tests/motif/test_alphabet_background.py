from __future__ import annotations

import numpy as np
import pytest

from omicsweave import motif, types

from ._support import DNA
from ._support import collection as collection


def test_alphabet_presets_round_trip():
    assert DNA.symbols == "ACGT" and DNA.complement == "TGCA"
    assert dict(DNA.aliases) == {"a": "A", "c": "C", "g": "G", "t": "T"}
    assert DNA.lut[ord("g")] == 2 and DNA.lut[ord("N")] == 255
    assert types.Alphabet.protein().complement is None


def test_background_large_frequencies_and_sample_lengths():
    bg = motif.Background.frequencies([1e308, 1e308, 1e308, 1e308])
    np.testing.assert_array_equal(bg.letters, np.full(4, 0.25))
    with pytest.raises(ValueError, match="nonnegative"):
        bg.sample(DNA, [-1])
    with pytest.raises(TypeError, match="integers"):
        bg.sample(DNA, [1.5])
    with pytest.raises(ValueError, match="one-dimensional"):
        bg.sample(DNA, [[1]])
    with pytest.raises(ValueError, match="offset or address space"):
        bg.sample(DNA, np.array([2**63], dtype=np.uint64))
    assert len(bg.sample(DNA, [])) == 0
    with pytest.raises(ValueError, match="nonnegative"):
        bg.truncated(-2)


def test_backgrounds():
    seqs = ["ACGTTTGA" * 20, "acgtnnACG"]
    bg = motif.Background.estimate(seqs, DNA, order=2, pseudocount=0.5)
    assert bg.order == 2 and [len(k) for k in bg.kmers] == [4, 16, 64]
    sym = bg.symmetrized(DNA)
    np.testing.assert_allclose(sym.kmers[0], (bg.kmers[0] + bg.kmers[0][::-1]) / 2)
    sample = bg.sample(DNA, [50, 0, 7], seed=1)
    assert [len(s) for s in sample] == [50, 0, 7]
    assert motif.Background.parse_meme(
        "A 0.3\nC 0.2\nG 0.2\nT 0.3\n", DNA
    ) == motif.Background.frequencies([0.3, 0.2, 0.2, 0.3])
    assert motif.Background.nrdb(DNA).order == 0
