"""A compiled motif opened up: its log-odds terms, their integer form and
the exact null distribution, for building thresholds and statistics
outside a scan."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .. import _native

__all__ = ["CompiledMotif", "NullDistribution", "Term"]


@dataclass(frozen=True, slots=True, eq=False)
class Term:
    """One additive term of a window score: a table indexed by the codes at
    ``positions`` (one axis per position; negative positions read the
    background context before the window). Real terms hold log-odds
    (``-inf`` forbids); integer terms hold their quantized cells
    (``int32``, the minimum forbids)."""

    positions: tuple[int, ...]
    table: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class NullDistribution:
    """The distribution of a window's integer score when the window and its
    context are drawn from the background: ``pmf[s] = P(S = s)`` and
    ``survival[s] = P(S >= s)`` for ``s = 0 .. max``; ``forbidden`` is the
    probability that a window has no score."""

    pmf: np.ndarray
    survival: np.ndarray
    forbidden: float


@dataclass(frozen=True, slots=True, eq=False)
class CompiledMotif:
    """Motif ``motif`` of a program on ``strand``, for windows with the full
    background context.

    A window's score is the sum of ``terms`` at its positions. The scan
    works with the integer score ``S`` (the sum of ``int_terms``), which
    stands for the real score ``S / scale + offset``; the reported score
    and p-value are those values kept to ``rounding`` significant digits.
    """

    motif: int
    strand: str
    width: int
    terms: tuple[Term, ...]
    int_terms: tuple[Term, ...]
    scale: float
    offset: float
    max: int
    null: NullDistribution
    rounding: int | None

    def _round(self, x: np.ndarray) -> np.ndarray:
        if self.rounding is None:
            return x
        flat = np.ascontiguousarray(x, dtype=np.float64).ravel()
        return _native.round_significant(flat, self.rounding).reshape(np.shape(x))

    def raw(self, s: int | np.ndarray) -> np.ndarray:
        """The real score of integer score(s) ``s``."""
        return np.asarray(s, dtype=np.float64) / self.scale + self.offset

    def score(self, s: int | np.ndarray) -> np.ndarray:
        """The reported score of integer score(s) ``s``."""
        return self._round(self.raw(s))

    def pvalue(self, s: int | np.ndarray) -> np.ndarray:
        """The reported p-value of integer score(s) ``s``: the scorable
        probability below 0, and 0 above ``max``. Forbidden windows have
        no finite score and contribute no probability."""
        s = np.asarray(s, dtype=np.int64)
        idx = np.clip(s, 0, self.max)
        p = np.where(s > self.max, 0.0, self.null.survival[idx])
        return self._round(p)

    def cutoff(self, pvalue: float) -> int | None:
        """The smallest integer score whose reported p-value is at most
        ``pvalue``; ``None`` when no score is that significant."""
        # Survival probabilities and their significant-digit rounding are
        # nonincreasing, so thresholds need no score-sized temporary array.
        lo, hi = 0, self.max + 1
        while lo < hi:
            mid = (lo + hi) // 2
            p = self.null.survival[mid]
            if self.rounding is not None:
                p = _native.round_significant(np.array([p], dtype=np.float64), self.rounding)[0]
            if p <= pvalue:
                hi = mid
            else:
                lo = mid + 1
        return lo if lo <= self.max else None

    def score_threshold(self, pvalue: float) -> float | None:
        """The reported score at :meth:`cutoff`."""
        s = self.cutoff(pvalue)
        return None if s is None else float(self.score(s))
