"""Bases by locus, from FASTA (plain or bgzip, with its indexes) or 2bit."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .. import _native
from ..cache import BlockCache, BlockKey
from ..cache import native as native_cache
from ..errors import UnknownChromError
from ..formats import BgzipFasta, Fai, Fasta, Gzi, PathLike, TwoBit
from ..types import LocusLike, Sequence
from ..types._genome import as_locus, native_key
from ._base import CachedReader, resolve

__all__ = ["SequenceReader"]


def _path(file: PathLike | None) -> str | None:
    return None if file is None else os.fspath(file)


class SequenceReader(CachedReader):
    """Reads the bases of a FASTA, bgzip FASTA or 2bit file by locus: a
    :class:`~omicsweave.types.Region` or :class:`~omicsweave.types.Spliced`
    spans, 5' to 3' of the key's strand (reverse complement on ``"-"``).

    FASTA is read through its ``fai`` and bgzip FASTA also through its
    ``gzi``; when not given they are the files beside it
    (:meth:`Fai.beside <omicsweave.formats.Fai.beside>`). A missing index
    raises :class:`~omicsweave.errors.IoError`: indexes are built
    explicitly, never on open. ``soft_mask`` lowercases the soft-masked
    bases of a 2bit file; FASTA bases come as stored.
    """

    def __init__(
        self,
        source: Fasta | BgzipFasta | TwoBit | PathLike,
        fai: Fai | PathLike | None = None,
        gzi: Gzi | PathLike | None = None,
        *,
        soft_mask: bool = False,
        cache: BlockCache | None = None,
    ) -> None:
        file = resolve(source, (Fasta, BgzipFasta, TwoBit), "SequenceReader")
        fai = None if fai is None else Fai(fai)
        gzi = None if gzi is None else Gzi(gzi)
        native = _native.SequenceReader(
            os.fspath(file), _path(fai), _path(gzi), soft_mask, native_cache(cache)
        )
        self._open(native, file, fai, gzi, soft_mask=soft_mask)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach caller-owned blocks; plain FASTA rejects a nonempty cache."""
        self._native.set_cache(native_cache(cache))

    def describe_cache(self, locus: LocusLike) -> tuple[BlockKey, ...]:
        """The block keys on which this locus depends, each once."""
        return tuple(self._native.describe_cache(native_key(as_locus(locus))))

    def prefetch_cache(self, locus: LocusLike) -> None:
        """Populate the attached cache without constructing sequence output."""
        self._native.prefetch_cache(native_key(as_locus(locus)))

    @property
    def chroms(self) -> Mapping[str, int]:
        """Sequence lengths by name, in file order."""
        return self._chroms

    def length(self, chrom: str) -> int:
        """The length of sequence ``chrom``."""
        try:
            return self._chroms[chrom]
        except KeyError:
            raise UnknownChromError(f"unknown chromosome {chrom!r}", chrom=chrom) from None

    def fetch(self, locus: LocusLike, out: np.ndarray | None = None) -> Sequence:
        """The bases of ``locus``, whose spans must lie within its sequence.

        By default, returns a new array independent of later queries and
        the reader's lifetime. Optional ``out`` (writable, contiguous
        ``uint8``, one element per position) receives the bases in place;
        the returned sequence references that same array.
        """
        key = as_locus(locus)
        native = native_key(key)
        if out is None:
            length = len(key)
            if length <= 4096:
                out = self._native.fetch(native)
            else:
                # Long reads amortize the call boundary; empty storage avoids
                # a separate zero-fill before the decoder overwrites it.
                out = np.empty(length, dtype=np.uint8)
                self._native.sequence(native, out)
        else:
            self._native.sequence(native, out)
        return Sequence(key, out)
