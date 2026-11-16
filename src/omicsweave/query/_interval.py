"""Features by region, from bigBed."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

from .. import _native
from ..cache import BlockCache
from ..cache import native as native_cache
from ..errors import InvalidError
from ..formats import Bed, BigBed, BigBedSchema, File, GzipBed, PathLike, detect
from ..types import Intervals, RegionLike, StringVector
from ..types._genome import as_region
from ._base import CachedReader, conversion, resolve

__all__ = ["IntervalReader"]


def _bigbed(source: BigBed | PathLike) -> BigBed:
    """``source`` as a bigBed: BED text is a conversion source only."""
    file = source if isinstance(source, File) else detect(source)
    detail = conversion(file) if isinstance(file, (Bed, GzipBed)) else None
    if detail is not None:
        if isinstance(source, File):
            raise TypeError(detail)
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    return resolve(source, (BigBed,), "IntervalReader")


class IntervalReader(CachedReader):
    """Reads the features of a bigBed file by region.

    BED text is not read by region: convert it first (``ow convert
    peaks.bed peaks.bb --chrom-sizes genome.sizes``).
    """

    def __init__(self, source: BigBed | PathLike, *, cache: BlockCache | None = None) -> None:
        file = _bigbed(source)
        native = _native.IntervalReader(os.fspath(file), native_cache(cache))
        self._open(native, file)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        fields, defined = native.field_count
        self._schema = BigBedSchema(fields, defined, native.auto_sql)

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach decoded bigBed blocks, or detach them with None."""
        self._native.set_cache(native_cache(cache))

    @property
    def chroms(self) -> Mapping[str, int]:
        """Chromosome lengths by name."""
        return self._chroms

    @property
    def schema(self) -> BigBedSchema:
        """The record layout of the bigBed file."""
        return self._schema

    def fetch(self, region: RegionLike) -> Intervals:
        """The features overlapping ``region``."""
        r = as_region(region)
        start, end, offsets, data = self._native.query(r.chrom, r.start, r.end)
        return Intervals(r, start, end, StringVector(offsets, data))
