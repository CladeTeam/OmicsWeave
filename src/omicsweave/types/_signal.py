"""Numeric signal along a genome, such as coverage or fold change."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from ._genome import Locus, Region

__all__ = ["Signal", "SignalBins", "SignalIntervals", "Stat"]


class Stat(StrEnum):
    """A statistic of the signal in a bin, as libBigWig's ``bwStats``
    computes it. Bins without data read NaN, except :attr:`COVERAGE`."""

    MEAN = "mean"
    """Mean over the bases with data. The mean with missing bases as 0 is
    ``SUM / bin width``."""
    MIN = "min"
    MAX = "max"
    SUM = "sum"
    """Sum over the bases with data of their value."""
    STD = "std"
    """Sample standard deviation over the bases with data; 0 for one base."""
    COVERAGE = "coverage"
    """Fraction of the bin's bases with data: 0 for an empty bin."""


@dataclass(frozen=True, slots=True, eq=False)
class Signal:
    """One value per position of ``locus`` (``float32``), 5' to 3'; NaN where
    there is no data."""

    locus: Locus
    values: np.ndarray

    def __len__(self) -> int:
        return len(self.values)


@dataclass(frozen=True, slots=True, eq=False)
class SignalBins:
    """Statistics of ``locus`` in bins of ``resolution`` positions.

    Bins tile the locus 5' to 3': local position ``x`` falls in bin
    ``(x + phase) // resolution``, so with a nonzero ``phase`` or a length
    that is not a multiple of ``resolution`` the first and last bins are
    partial. ``values`` is ``float32`` of shape ``(len(stats), bins)``, one
    row per statistic in ``stats`` order; ``bins[stat]`` is a row.
    """

    locus: Locus
    resolution: int
    phase: int
    stats: tuple[Stat, ...]
    values: np.ndarray

    def __len__(self) -> int:
        return self.values.shape[1]

    def __getitem__(self, stat: Stat | str) -> np.ndarray:
        try:
            return self.values[self.stats.index(Stat(stat))]
        except ValueError:
            raise KeyError(stat) from None

    @property
    def edges(self) -> np.ndarray:
        """Local positions (``int64``) where the bins start, then the end:
        bin ``k`` covers ``[edges[k], edges[k + 1])``."""
        k = np.arange(len(self) + 1, dtype=np.int64)
        return np.clip(k * self.resolution - self.phase, 0, len(self.locus))


@dataclass(frozen=True, slots=True, eq=False)
class SignalIntervals:
    """The runs of constant signal in ``region``, in position order: cut to
    the region unless the query asked for the stored bounds.

    Columns: ``start`` and ``end`` (``uint32``), ``value`` (``float32``).
    Bases without data lie in no run.
    """

    region: Region
    start: np.ndarray
    end: np.ndarray
    value: np.ndarray

    def __len__(self) -> int:
        return len(self.start)
