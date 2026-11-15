"""Stretches of a genome: where they are, their bases, the features on them."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from ..errors import InvalidError, InvalidRegionError
from ._string import StringVector

__all__ = [
    "STRANDS",
    "Intervals",
    "Locus",
    "LocusLike",
    "Region",
    "RegionLike",
    "Sequence",
    "Spliced",
    "as_locus",
    "as_region",
]

U32_MAX = 2**32 - 1

STRANDS = ("+", "-", ".", "?")
"""Strand spellings: forward, reverse, not stranded, stranded but unknown."""


def _coord(name: str, value: object) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError(f"{name} must be an integer")
    if not 0 <= value <= U32_MAX:
        raise ValueError(f"{name} {value} is outside [0, 2**32)")
    return int(value)


def _strand(value: object) -> str:
    if value not in STRANDS:
        raise ValueError(f"strand must be one of + - . ?, not {value!r}")
    return str(value)


@dataclass(frozen=True, slots=True)
class Region:
    """Bases ``[start, end)`` of chromosome ``chrom``: 0-based, half-open,
    on ``strand``.

    The key every genomic lookup takes. Dense results (bases, per-base
    signal, bins, contact tiles) run 5' to 3' of the strand: on ``"-"``
    they are reversed, and bases complemented. Sparse results (features,
    reads, records, pixels) keep genomic coordinates and order whatever the
    strand. ``"."`` and ``"?"`` read as ``"+"``.
    """

    chrom: str
    start: int
    end: int
    strand: str = "."

    def __post_init__(self) -> None:
        if not isinstance(self.chrom, str):
            raise TypeError(f"chrom must be str, not {type(self.chrom).__name__}")
        object.__setattr__(self, "start", _coord("start", self.start))
        object.__setattr__(self, "end", _coord("end", self.end))
        object.__setattr__(self, "strand", _strand(self.strand))
        if self.start > self.end:
            raise InvalidRegionError(
                f"invalid region {self.start}..{self.end}", start=self.start, end=self.end
            )

    def __len__(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        tail = "" if self.strand == "." else f":{self.strand}"
        return f"{self.chrom}:{self.start}-{self.end}{tail}"

    @property
    def is_reverse(self) -> bool:
        return self.strand == "-"

    def with_strand(self, strand: str) -> Region:
        """This region on ``strand``."""
        return Region(self.chrom, self.start, self.end, strand)


RegionLike = Region | tuple[str, int, int] | tuple[str, int, int, str]
"""A :class:`Region`, or its ``(chrom, start, end)`` or ``(chrom, start, end,
strand)``."""


def as_region(key: RegionLike) -> Region:
    """``key`` as a :class:`Region`."""
    if isinstance(key, Region):
        return key
    if isinstance(key, tuple) and len(key) in (3, 4):
        return Region(*key)
    if isinstance(key, Spliced):
        raise TypeError(
            "a spliced key reads dense data only (bases, signal, bins); "
            "query each span as a Region instead"
        )
    raise TypeError(f"expected a Region or (chrom, start, end[, strand]), not {key!r}")


@dataclass(frozen=True, slots=True, init=False)
class Spliced:
    """Ascending, disjoint spans of chromosome ``chrom`` read as one
    sequence 5' to 3' of ``strand``: a transcript's exons, its coding
    part, or any set of spans.

    Dense readers return ``len(key)`` positions: the spans concatenated,
    reversed (and bases complemented) on ``"-"``. Spans may touch; empty
    spans contribute nothing. ``name`` labels the key (a transcript id).
    """

    chrom: str
    spans: tuple[tuple[int, int], ...]
    strand: str
    name: str | None

    def __init__(
        self,
        chrom: str,
        spans: Iterable[tuple[int, int]] | np.ndarray,
        strand: str = ".",
        name: str | None = None,
    ) -> None:
        if not isinstance(chrom, str):
            raise TypeError(f"chrom must be str, not {type(chrom).__name__}")
        pairs = []
        floor = 0
        for i, span in enumerate(
            np.asarray(spans).reshape(-1, 2) if isinstance(spans, np.ndarray) else spans
        ):
            start, end = (_coord("start", span[0]), _coord("end", span[1]))
            if start > end:
                raise InvalidRegionError(f"invalid span {start}..{end}", start=start, end=end)
            if i and start < floor:
                detail = f"span {start}..{end} starts before the end {floor} of the span before it"
                raise InvalidError(f"spliced key: {detail}", what="spliced key", detail=detail)
            floor = end
            pairs.append((start, end))
        object.__setattr__(self, "chrom", chrom)
        object.__setattr__(self, "spans", tuple(pairs))
        object.__setattr__(self, "strand", _strand(strand))
        object.__setattr__(self, "name", name)

    def __len__(self) -> int:
        return sum(e - s for s, e in self.spans)

    def __str__(self) -> str:
        spans = ",".join(f"{s}-{e}" for s, e in self.spans)
        tail = "" if self.strand == "." else f":{self.strand}"
        return f"{self.chrom}:[{spans}]{tail}"

    @property
    def is_reverse(self) -> bool:
        return self.strand == "-"

    @property
    def hull(self) -> Region:
        """The smallest region holding every span, on the key's strand."""
        if not self.spans:
            return Region(self.chrom, 0, 0, self.strand)
        return Region(self.chrom, self.spans[0][0], self.spans[-1][1], self.strand)

    def _layout(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        spans = np.asarray(self.spans, dtype=np.int64).reshape(-1, 2)
        starts, ends = spans[:, 0], spans[:, 1]
        offsets = np.concatenate(([0], np.cumsum(ends - starts)))
        return starts, ends, offsets

    def to_genome(self, positions: np.ndarray | Iterable[int]) -> np.ndarray:
        """Genomic coordinates (``int64``) of output positions ``0 <= i <
        len(key)``, 5' to 3' of the strand."""
        i = np.asarray(positions, dtype=np.int64)
        n = len(self)
        if np.any((i < 0) | (i >= n)):
            raise IndexError(f"positions must lie in [0, {n})")
        starts, _, offsets = self._layout()
        forward = n - 1 - i if self.is_reverse else i
        k = np.searchsorted(offsets, forward, side="right") - 1
        return starts[k] + (forward - offsets[k])

    def to_local(self, positions: np.ndarray | Iterable[int]) -> np.ndarray:
        """Output positions (``int64``) of genomic coordinates; -1 for those
        outside every span."""
        g = np.asarray(positions, dtype=np.int64)
        starts, ends, offsets = self._layout()
        k = np.searchsorted(starts, g, side="right") - 1
        inside = k >= 0
        k = np.clip(k, 0, max(len(starts) - 1, 0))
        if len(starts):
            inside &= g < ends[k]
        forward = np.where(inside, offsets[k] + (g - starts[k] if len(starts) else 0), -1)
        if self.is_reverse:
            forward = np.where(inside, len(self) - 1 - forward, -1)
        return forward


Locus = Region | Spliced
"""The key of a dense query: a region, or spliced spans."""

LocusLike = RegionLike | Spliced
"""A :data:`Locus`, or a region's tuple."""


def as_locus(key: LocusLike) -> Locus:
    """``key`` as a :data:`Locus`."""
    if isinstance(key, Spliced):
        return key
    return as_region(key)


def native_key(locus: Locus) -> tuple[str, list[tuple[int, int]], str]:
    """The ``(chrom, spans, strand)`` the native readers take."""
    if isinstance(locus, Region):
        return (locus.chrom, [(locus.start, locus.end)], locus.strand)
    return (locus.chrom, list(locus.spans), locus.strand)


@dataclass(frozen=True, slots=True, eq=False)
class Sequence:
    """The bases of ``locus`` as ASCII bytes (``uint8``), 5' to 3'."""

    locus: Locus
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
