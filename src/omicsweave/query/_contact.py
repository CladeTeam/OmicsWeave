"""Contacts by window, from ``.hic``: raw counts, never normalized."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .. import _native
from ..cache import BlockCache, BlockKey
from ..cache import native as native_cache
from ..formats import Hic, PathLike
from ..types import BinUnit, ContactMatrix, ContactPixels, ContactWindow
from ._base import CachedReader, resolve

__all__ = ["ContactReader"]


def _tile(w: ContactWindow) -> tuple[object, ...]:
    return (
        w.rows.chrom,
        w.rows.start,
        w.rows.end,
        w.cols.chrom,
        w.cols.start,
        w.cols.end,
        w.resolution,
        str(w.unit),
    )


class ContactReader(CachedReader):
    """Reads the contact maps of a ``.hic`` file by window."""

    def __init__(self, source: Hic | PathLike, *, cache: BlockCache | None = None) -> None:
        file = resolve(source, (Hic,), "ContactReader")
        native = _native.ContactReader(os.fspath(file), native_cache(cache))
        self._open(native, file)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))
        self._attributes = MappingProxyType(dict(native.attributes()))

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach external decoded blocks or return to direct reads."""
        self._native.set_cache(native_cache(cache))

    def describe_cache(self, window: ContactWindow) -> tuple[BlockKey, ...]:
        """The source-bound blocks required by this window."""
        return tuple(self._native.describe_cache(_tile(window)))

    def prefetch_cache(self, window: ContactWindow) -> None:
        """Prepare reusable blocks without producing a matrix or pixels."""
        self._native.prefetch_cache(_tile(window))

    @property
    def version(self) -> int:
        return self._native.version

    @property
    def genome(self) -> str:
        return self._native.genome

    @property
    def attributes(self) -> Mapping[str, str]:
        return self._attributes

    @property
    def chroms(self) -> Mapping[str, int]:
        """Chromosome lengths (base pairs) by name, in file order."""
        return self._chroms

    def resolutions(self, unit: BinUnit | str = BinUnit.BP) -> tuple[int, ...]:
        """The bin sizes stored in ``unit``."""
        return tuple(self._native.resolutions(str(BinUnit(unit))))

    def fragment_count(self, chrom: str) -> int:
        return self._native.fragment_count(chrom)

    def sites(self, chrom: str) -> np.ndarray:
        """The restriction site positions of ``chrom`` (``uint32``)."""
        return self._native.sites(chrom)

    def matrix(self, window: ContactWindow, out: np.ndarray | None = None) -> ContactMatrix:
        """The contacts of ``window``, whole, in the orientation asked for:
        rows (columns) run 5' to 3' of the strand of ``window.rows``
        (``window.cols``), so a ``"-"`` region reverses its axis.

        By default, returns a new array independent of later queries and
        the reader's lifetime. Optional ``out`` (writable, C-contiguous
        ``float32`` of ``window.shape``) receives the values in place;
        the returned matrix references that same array.
        """
        if out is None:
            out = np.empty(window.shape, dtype=np.float32)
        self._native.matrix(_tile(window), out, window.rows.is_reverse, window.cols.is_reverse)
        return ContactMatrix(window, out)

    def pixels(self, window: ContactWindow) -> ContactPixels:
        """The nonzero pixels of ``window``, at genomic bin coordinates
        whatever the strands."""
        return ContactPixels(window, *self._native.cells(_tile(window)))
