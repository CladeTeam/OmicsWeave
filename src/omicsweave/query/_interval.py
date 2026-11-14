"""Features by region, from bigBed or BED."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

from .. import _native
from ..formats import Bed, BigBed, BigBedSchema, GzipBed, PathLike
from ..types import Intervals, RegionLike, StringVector
from ..types._genome import as_region
from ._base import Reader, resolve

__all__ = ["IntervalReader"]


class IntervalReader(Reader):
    """Reads the features of a bigBed file, or of a BED file (plain or gzip)
    loaded whole, by region."""

    def __init__(self, source: BigBed | Bed | GzipBed | PathLike) -> None:
        file = resolve(source, (BigBed, Bed, GzipBed), "IntervalReader")
        native = _native.IntervalReader(os.fspath(file))
        self._open(native, file)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        if native.field_count is None:
            self._schema = None
        else:
            fields, defined = native.field_count
            self._schema = BigBedSchema(fields, defined, native.auto_sql)

    @property
    def chroms(self) -> Mapping[str, int | None]:
        """Chromosome lengths by name; BED files record none (``None``)."""
        return self._chroms

    @property
    def schema(self) -> BigBedSchema | None:
        """The record layout of a bigBed file; BED files declare none."""
        return self._schema

    def fetch(self, region: RegionLike) -> Intervals:
        """The features overlapping ``region``."""
        r = as_region(region)
        start, end, offsets, data = self._native.query(r.chrom, r.start, r.end)
        return Intervals(r, start, end, StringVector(offsets, data))
