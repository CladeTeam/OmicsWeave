"""Features by region, from bigBed."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

from .. import _native
from ..cache import BlockCache, BlockKey
from ..cache import native as native_cache
from ..errors import InvalidError
from ..formats import Bed, BigBed, BigBedSchema, File, GzipBed, PathLike, detect
from ..types import Intervals, RegionLike, StringVector
from ..types._genome import as_region
from ._base import CachedReader, resolve

__all__ = ["IntervalReader"]


def _bigbed(source: BigBed | PathLike) -> BigBed:
    """``source`` as a bigBed: BED text is a conversion source only."""
    file = source if isinstance(source, File) else detect(source)
    if isinstance(file, (Bed, GzipBed)):
        detail = (
            "IntervalReader reads bigBed by region, not BED text; convert it with "
            f"`ow convert {file.path} {_bigbed_name(file)} --chrom-sizes <genome.sizes>`"
        )
        if isinstance(source, File):
            raise TypeError(detail)
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    return resolve(source, (BigBed,), "IntervalReader")


def _bigbed_name(file: Bed | GzipBed) -> str:
    name = os.fspath(file.path)
    for suffix in (".bed.gz", ".bed"):
        if name.endswith(suffix):
            return name[: -len(suffix)] + ".bb"
    return name + ".bb"


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

    def describe_cache(self, region: RegionLike) -> tuple[BlockKey, ...]:
        r = as_region(region)
        return tuple(self._native.describe_cache(r.chrom, r.start, r.end))

    def prefetch_cache(self, region: RegionLike) -> None:
        r = as_region(region)
        self._native.prefetch_cache(r.chrom, r.start, r.end)

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
