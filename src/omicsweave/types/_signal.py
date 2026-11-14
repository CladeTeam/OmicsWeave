"""Numeric signal along a genome, such as coverage or fold change."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._genome import Region

__all__ = ["Signal", "SignalIntervals"]


@dataclass(frozen=True, slots=True, eq=False)
class Signal:
    """One value per base of ``region`` (``float32``); NaN where there is no data."""

    region: Region
    values: np.ndarray

    def __len__(self) -> int:
        return len(self.values)


@dataclass(frozen=True, slots=True, eq=False)
class SignalIntervals:
    """The runs of constant signal in ``region``, clipped to it, in position order.

    Columns: ``start`` and ``end`` (``uint32``), ``value`` (``float32``).
    Bases without data lie in no run.
    """

    region: Region
    start: np.ndarray
    end: np.ndarray
    value: np.ndarray

    def __len__(self) -> int:
        return len(self.start)
