"""Signal by region, from bigWig."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .. import _native
from ..formats import BigWig, BigWigSummary, PathLike
from ..types import RegionLike, Signal, SignalIntervals
from ..types._genome import as_region
from ._base import Reader, resolve

__all__ = ["TrackReader"]


class TrackReader(Reader):
    """Reads the signal of a bigWig file by region."""

    def __init__(self, source: BigWig | PathLike) -> None:
        file = resolve(source, (BigWig,), "TrackReader")
        native = _native.TrackReader(os.fspath(file))
        self._open(native, file)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        summary = native.summary()
        self._summary = None if summary is None else BigWigSummary(*summary)

    @property
    def chroms(self) -> Mapping[str, int]:
        """Chromosome lengths by name."""
        return self._chroms

    @property
    def summary(self) -> BigWigSummary | None:
        return self._summary

    def signal(self, region: RegionLike, out: np.ndarray | None = None) -> Signal:
        """One value per base of ``region``; it may run past its chromosome's end.

        ``out`` (``float32``, one element per base) receives them in place.
        """
        r = as_region(region)
        if out is None:
            out = np.empty(len(r), dtype=np.float32)
        self._native.values(r.chrom, r.start, r.end, out)
        return Signal(r, out)

    def intervals(self, region: RegionLike) -> SignalIntervals:
        """The runs of constant signal in ``region``, as stored."""
        r = as_region(region)
        return SignalIntervals(r, *self._native.intervals(r.chrom, r.start, r.end))
