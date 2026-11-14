"""Stretches of a genome: where they are, their bases, the features on them."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..errors import InvalidRegionError
from ._string import StringVector

__all__ = ["Intervals", "Region", "RegionLike", "Sequence", "as_region"]

U32_MAX = 2**32 - 1


@dataclass(frozen=True, slots=True)
class Region:
    """Bases ``[start, end)`` of chromosome ``chrom``: 0-based, half-open.

    The key every genomic lookup takes.
    """

    chrom: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if not isinstance(self.chrom, str):
            raise TypeError(f"chrom must be str, not {type(self.chrom).__name__}")
        for name in ("start", "end"):
            value = getattr(self, name)
            if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
                raise TypeError(f"{name} must be an integer")
            if not 0 <= value <= U32_MAX:
                raise ValueError(f"{name} {value} is outside [0, 2**32)")
            object.__setattr__(self, name, int(value))
        if self.start > self.end:
            raise InvalidRegionError(
                f"invalid region {self.start}..{self.end}", start=self.start, end=self.end
            )

    def __len__(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


RegionLike = Region | tuple[str, int, int]
"""A :class:`Region`, or its ``(chrom, start, end)``."""


def as_region(key: RegionLike) -> Region:
    """``key`` as a :class:`Region`."""
    if isinstance(key, Region):
        return key
    if isinstance(key, tuple) and len(key) == 3:
        return Region(*key)
    raise TypeError(f"expected a Region or (chrom, start, end), not {key!r}")


@dataclass(frozen=True, slots=True, eq=False)
class Sequence:
    """The bases of ``region`` as ASCII bytes (``uint8``)."""

    region: Region
    bases: np.ndarray

    def __len__(self) -> int:
        return len(self.bases)

    @property
    def text(self) -> str:
        return self.bases.tobytes().decode("ascii")


@dataclass(frozen=True, slots=True, eq=False)
class Intervals:
    """The features overlapping ``region``, unclipped, ordered by start.

    Columns: ``start`` and ``end`` (``uint32``), and ``annotation``: each
    feature's fields after the third, tab-separated.
    """

    region: Region
    start: np.ndarray
    end: np.ndarray
    annotation: StringVector

    def __len__(self) -> int:
        return len(self.start)
