"""Bases by region, from FASTA (plain or bgzip, with its indexes) or 2bit."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .. import _native
from ..errors import UnknownChromError
from ..formats import BgzipFasta, Fai, Fasta, Gzi, PathLike, TwoBit
from ..types import RegionLike, Sequence
from ..types._genome import as_region
from ._base import Reader, resolve

__all__ = ["SequenceReader"]


def _path(file: PathLike | None) -> str | None:
    return None if file is None else os.fspath(file)


class SequenceReader(Reader):
    """Reads the bases of a FASTA, bgzip FASTA or 2bit file by region.

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
    ) -> None:
        file = resolve(source, (Fasta, BgzipFasta, TwoBit), "SequenceReader")
        fai = None if fai is None else Fai(fai)
        gzi = None if gzi is None else Gzi(gzi)
        native = _native.SequenceReader(os.fspath(file), _path(fai), _path(gzi), soft_mask)
        self._open(native, file, fai, gzi, soft_mask=soft_mask)
        self.source = file
        self._chroms = MappingProxyType(dict(native.chroms()))

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

    def fetch(self, region: RegionLike, out: np.ndarray | None = None) -> Sequence:
        """The bases of ``region``, which must lie within its sequence.

        ``out`` (``uint8``, one element per base) receives them in place.
        """
        r = as_region(region)
        if out is None:
            out = np.empty(len(r), dtype=np.uint8)
        self._native.sequence(r.chrom, r.start, r.end, out)
        return Sequence(r, out)
