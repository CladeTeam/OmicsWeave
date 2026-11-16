"""Multiple alignments by region, from bigMaf."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from types import MappingProxyType

import numpy as np

from .. import _native
from ..cache import BlockCache
from ..cache import native as native_cache
from ..errors import InvalidError
from ..formats import BigBed, File, PathLike, detect
from ..types import MultipleAlignment, RegionLike
from ..types._genome import as_region
from ._base import CachedReader, conversion, resolve

__all__ = ["MsaReader"]


def _fill_byte(fill: str | bytes | int) -> int:
    if isinstance(fill, int):
        if not 0 <= fill < 256:
            raise ValueError(f"fill {fill} is not a byte")
        return fill
    raw = fill.encode("ascii") if isinstance(fill, str) else bytes(fill)
    if len(raw) != 1:
        raise ValueError(f"fill must be one ASCII character, not {fill!r}")
    return raw[0]


def _bigmaf(source: BigBed | PathLike) -> BigBed:
    """``source`` as a bigBed: MAF text is a conversion source only."""
    file = source if isinstance(source, File) else detect(source)
    detail = conversion(file)
    if detail is not None:
        if isinstance(source, File):
            raise TypeError(detail)
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    return resolve(source, (BigBed,), "MsaReader")


class MsaReader(CachedReader):
    """Reads a bigMaf (the bigBed ``ow convert`` or UCSC's ``mafToBigMaf``
    and ``bedToBigBed`` write from MAF) by region, as a
    :class:`~omicsweave.types.MultipleAlignment`: one row per species, one
    column per base of the reference, the species whose rows open every
    block.

    A species is the ``src`` of its rows before the first dot (``hg38`` of
    ``hg38.chr1``). ``species`` fixes the rows and their order; species not
    named are skipped, and named ones a block lacks hold ``fill``. Without
    it, the rows are the species of the blocks read, in order of first
    appearance, so they can differ from query to query. Where a block
    holds more than one row of a species, the first is read. Positions no
    block covers hold ``fill`` in every row. Bases keep their case.

    MAF text is not read by region: convert it first (``ow convert
    multiz.maf.gz multiz.bigmaf``).
    """

    def __init__(
        self,
        source: BigBed | PathLike,
        species: Iterable[str] | None = None,
        *,
        fill: str | bytes | int = "-",
        cache: BlockCache | None = None,
    ) -> None:
        file = _bigmaf(source)
        names = None if species is None else tuple(species)
        if names is not None and len(set(names)) != len(names):
            raise ValueError("species are named more than once")
        fill_byte = _fill_byte(fill)
        native = _native.MsaReader(os.fspath(file), names, fill_byte, native_cache(cache))
        self._open(native, file, species=names, fill=fill_byte)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        self._species = names
        self._fill = fill_byte

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach decoded bigBed blocks, or detach them with None."""
        self._native.set_cache(native_cache(cache))

    @property
    def chroms(self) -> Mapping[str, int]:
        """Reference chromosome lengths by name."""
        return self._chroms

    @property
    def species(self) -> tuple[str, ...] | None:
        """The rows every query returns; ``None`` when they follow the blocks."""
        return self._species

    def fetch(self, region: RegionLike, out: np.ndarray | None = None) -> MultipleAlignment:
        """The alignment over ``region``.

        Optional ``out`` (writable, C-contiguous ``uint8`` of shape
        ``(len(species), len(region))``; needs ``species``) receives the
        bases in place; the result references it.
        """
        r = as_region(region)
        if out is None:
            species, bases = self._native.fetch(r.chrom, r.start, r.end, r.is_reverse)
            return MultipleAlignment(r, tuple(species), bases)
        if self._species is None:
            raise ValueError("out needs the species fixed when the reader opens")
        shape = (len(self._species), len(r))
        if (
            not isinstance(out, np.ndarray)
            or out.dtype != np.uint8
            or out.shape != shape
            or not out.flags.c_contiguous
            or not out.flags.writeable
        ):
            raise ValueError(f"out must be writable C-contiguous uint8 of shape {shape}")
        self._native.fetch_into(r.chrom, r.start, r.end, r.is_reverse, out)
        return MultipleAlignment(r, self._species, out)
