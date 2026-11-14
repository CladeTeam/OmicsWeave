"""Chromatin contacts: a window of a binned contact map and what it holds."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from ..errors import InvalidError
from ._genome import U32_MAX, Region, RegionLike, as_region

__all__ = ["BinUnit", "ContactMatrix", "ContactPixels", "ContactWindow"]


class BinUnit(StrEnum):
    """What a contact-map bin coordinate counts."""

    BP = "BP"
    """Base pairs: a region is a base range."""
    FRAG = "FRAG"
    """Restriction fragments: a region is a range of fragment indices."""


def _bins(r: Region, resolution: int) -> range:
    first = r.start // resolution
    if len(r) == 0:
        return range(first, first)
    return range(first, -(-r.end // resolution))


@dataclass(frozen=True, slots=True, init=False)
class ContactWindow:
    """The bins of ``rows`` against the bins of ``cols``, ``resolution`` units wide.

    A region covers every bin it overlaps; the two may lie on different
    chromosomes.
    """

    rows: Region
    cols: Region
    resolution: int
    unit: BinUnit

    def __init__(
        self,
        rows: RegionLike,
        cols: RegionLike,
        resolution: int,
        unit: BinUnit | str = BinUnit.BP,
    ) -> None:
        if not isinstance(resolution, (int, np.integer)) or isinstance(resolution, bool):
            raise TypeError("resolution must be an integer")
        resolution = int(resolution)
        if not 0 < resolution <= U32_MAX:
            raise InvalidError(
                "contact window: resolution must be positive",
                what="contact window",
                detail="resolution must be positive",
            )
        object.__setattr__(self, "rows", as_region(rows))
        object.__setattr__(self, "cols", as_region(cols))
        object.__setattr__(self, "resolution", resolution)
        object.__setattr__(self, "unit", BinUnit(unit))

    @property
    def row_bins(self) -> range:
        return _bins(self.rows, self.resolution)

    @property
    def col_bins(self) -> range:
        return _bins(self.cols, self.resolution)

    @property
    def shape(self) -> tuple[int, int]:
        """(rows, columns) of bins."""
        return len(self.row_bins), len(self.col_bins)


@dataclass(frozen=True, slots=True, eq=False)
class ContactMatrix:
    """The contacts of ``window`` as a row-major ``float32`` matrix of
    ``window.shape``; 0 where there is no contact."""

    window: ContactWindow
    values: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class ContactPixels:
    """The nonzero pixels of ``window``, ordered by (row, col).

    Columns: ``row`` and ``col`` (``uint32``), the starts of each pixel's row
    and column bins in the window's unit, and ``value`` (``float32``).
    """

    window: ContactWindow
    row: np.ndarray
    col: np.ndarray
    value: np.ndarray

    def __len__(self) -> int:
        return len(self.row)
