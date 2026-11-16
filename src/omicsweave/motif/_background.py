"""Background models: the null the motif scores are contrasted with."""

from __future__ import annotations

from collections.abc import Sequence as Seq
from dataclasses import dataclass

import numpy as np

from .. import _native
from ..types import Alphabet, StringVector
from ._sequences import Sequences, as_strings

__all__ = ["Background"]


@dataclass(frozen=True, slots=True, eq=False)
class Background:
    """A Markov background of order ``len(kmers) - 1``.

    ``kmers[k - 1]`` holds the positive frequencies of every k-mer
    (``alphabet ** k`` values, first symbol most significant) for ``k = 1``
    to ``order + 1``; they need not be normalized. A symbol's probability
    after a context of ``c`` symbols is its ``(c + 1)``-mer frequency
    normalized over the last symbol, so a model reads with the longest
    context available and restarts with lower orders after a sequence
    start or a symbol outside the alphabet.

    Nothing is inferred: a background is always given explicitly, from
    frequencies, a MEME background file, a preset, or sequences.
    """

    kmers: tuple[np.ndarray, ...]

    def __post_init__(self) -> None:
        kmers = tuple(np.ascontiguousarray(k, dtype=np.float64) for k in self.kmers)
        if not kmers:
            raise ValueError("a background has at least the letter frequencies")
        object.__setattr__(self, "kmers", kmers)

    @property
    def order(self) -> int:
        return len(self.kmers) - 1

    @property
    def alphabet_len(self) -> int:
        return len(self.kmers[0])

    @property
    def letters(self) -> np.ndarray:
        """The letter (1-mer) probabilities, normalized."""
        letters = self.kmers[0]
        with np.errstate(over="ignore"):
            total = letters.sum()
        if np.isinf(total) and np.all(np.isfinite(letters)) and np.all(letters > 0):
            letters = letters / letters.max()
            total = letters.sum()
        return letters / total

    @classmethod
    def frequencies(cls, letters: Seq[float] | np.ndarray) -> Background:
        """The 0-order model of the given letter frequencies."""
        return cls((np.asarray(letters, dtype=np.float64),))

    @classmethod
    def uniform(cls, alphabet: Alphabet) -> Background:
        n = len(alphabet)
        return cls((np.full(n, 1.0 / n),))

    @classmethod
    def nrdb(cls, alphabet: Alphabet) -> Background:
        """MEME's non-redundant database frequencies (``ACGT`` or the 20
        amino acids), normalized as MEME does."""
        return cls(tuple(_native.background_nrdb(alphabet)))

    @classmethod
    def glam2(cls, alphabet: Alphabet) -> Background:
        """GLAM2's own: uniform, or Robinson & Robinson abundances for
        amino acids."""
        return cls(tuple(_native.background_glam2(alphabet)))

    @classmethod
    def parse_meme(cls, text: str, alphabet: Alphabet) -> Background:
        """The model of a MEME background file's text: lines of a symbol
        chain and its frequency, chains of each length listed completely
        before longer ones."""
        return cls(tuple(_native.background_meme(alphabet, text)))

    @classmethod
    def estimate(
        cls, sequences: Sequences, alphabet: Alphabet, *, order: int = 0, pseudocount: float = 1.0
    ) -> Background:
        """Counts every k-mer up to ``order + 1`` over runs of alphabet
        symbols, each count plus ``pseudocount``."""
        s = as_strings(sequences)
        return cls(
            tuple(_native.background_estimate(alphabet, s.offsets, s.data, order, pseudocount))
        )

    def symmetrized(self, alphabet: Alphabet) -> Background:
        """Every k-mer frequency averaged with its reverse complement's."""
        return Background(tuple(_native.background_symmetrized(list(self.kmers), alphabet)))

    def truncated(self, order: int) -> Background:
        """The same model cut to ``order`` (unchanged when already lower)."""
        if not isinstance(order, (int, np.integer)) or isinstance(order, bool):
            raise TypeError("background order must be an integer")
        if order < 0:
            raise ValueError("background order must be nonnegative")
        return Background(self.kmers[: order + 1])

    def sample(
        self, alphabet: Alphabet, lengths: Seq[int] | np.ndarray, *, seed: int = 0
    ) -> StringVector:
        """Sequences of the given lengths drawn from the model."""
        lengths = np.asarray(lengths)
        if lengths.ndim != 1:
            raise ValueError("sequence lengths must be one-dimensional")
        if lengths.size and lengths.dtype.kind not in "iu":
            raise TypeError("sequence lengths must be integers")
        if np.any(lengths < 0):
            raise ValueError("sequence lengths must be nonnegative")
        lengths = np.ascontiguousarray(lengths, dtype=np.uint64)
        offsets, data = _native.background_sample(list(self.kmers), alphabet, lengths, seed)
        return StringVector(offsets, data)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Background):
            return NotImplemented
        return len(self.kmers) == len(other.kmers) and all(
            np.array_equal(a, b) for a, b in zip(self.kmers, other.kmers, strict=True)
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"Background(order={self.order}, letters={np.round(self.letters, 4).tolist()})"
