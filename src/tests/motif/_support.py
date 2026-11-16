from __future__ import annotations

import numpy as np
import pytest

from omicsweave import formats, types

DNA = types.Alphabet.dna()


def random_pfm(rng: np.random.Generator, id: str, width: int) -> types.Motif:
    rows = rng.gamma(0.4, size=(width, 4)) + 1e-3
    return types.Motif.from_pfm(
        id, rows / rows.sum(1, keepdims=True), sites=float(rng.integers(5, 40))
    )


def dinucleotide(rng: np.random.Generator, id: str, width: int) -> types.Motif:
    factors = [
        types.Factor("joint", (i, i + 1), rng.gamma(0.5, size=16) + 1e-3) for i in range(width - 1)
    ]
    return types.Motif(id, width, tuple(factors), sites=20.0)


def random_sequences(rng: np.random.Generator, n: int, length: int) -> list[str]:
    letters = np.frombuffer(b"ACGTACGTACGTacgtN", dtype=np.uint8)
    return [
        letters[rng.integers(0, len(letters), size=length)].tobytes().decode() for _ in range(n)
    ]


@pytest.fixture
def collection(tmp_path):
    rng = np.random.default_rng(7)
    motifs = [random_pfm(rng, f"M{i}", int(rng.integers(5, 12))) for i in range(5)]
    motifs[1] = types.Motif.from_pfm(
        "M1",
        motifs[1].matrix,
        sites=motifs[1].sites,
        annotation={"name": "one", "rank": 3, "score": 0.5},
    )
    motifs.append(dinucleotide(rng, "D0", 6))
    path = tmp_path / "zoo.ffm"
    with formats.FfmWriter(path, DNA, {"name": "string", "rank": "integer", "score": "float"}) as w:
        w.add_all(motifs)
    return path, motifs
