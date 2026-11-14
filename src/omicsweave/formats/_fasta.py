"""FASTA, stored three ways: plain text (``.fa``), bgzip (``.fa.gz``) and
plain gzip.

Plain and bgzip FASTA are read by region through a
:class:`~omicsweave.formats.Fai`; bgzip FASTA also needs a
:class:`~omicsweave.formats.Gzi`. Plain-gzip FASTA cannot be read by region,
only converted. :class:`FastaWriter` writes plain or bgzip FASTA, with its
indexes.
"""

from __future__ import annotations

import os

import numpy as np

from .. import _native
from ._base import File, PathLike, Writer, ensure_format
from ._fai import Fai
from ._gzi import Gzi

__all__ = ["BgzipFasta", "Fasta", "FastaWriter", "GzipFasta"]


class Fasta(File):
    """Plain-text FASTA (``.fa``, ``.fasta``, ``.fna``)."""

    __slots__ = ()
    format = "fasta"

    def build_index(self) -> Fai:
        """Builds the ``.fai`` beside the file, as ``samtools faidx`` does."""
        ensure_format(self)
        (fai,) = _native.build_fasta_index(os.fspath(self))
        return Fai(fai)


class GzipFasta(File):
    """FASTA compressed with plain gzip (not BGZF): it has no random access,
    so it is only read whole, by a conversion to bgzip FASTA or 2bit."""

    __slots__ = ()
    format = "gzip-fasta"


class BgzipFasta(File):
    """FASTA compressed with bgzip (``.fa.gz``), so it can be read by region."""

    __slots__ = ()
    format = "bgzip-fasta"

    def build_index(self) -> tuple[Fai, Gzi]:
        """Builds the ``.fai`` and ``.gzi`` beside the file in one pass, as
        ``samtools faidx`` and ``bgzip -r`` do."""
        ensure_format(self)
        fai, gzi = _native.build_fasta_index(os.fspath(self))
        return Fai(fai), Gzi(gzi)


def _as_bytes(seq: str | bytes | np.ndarray) -> bytes:
    if isinstance(seq, str):
        return seq.encode("ascii")
    if isinstance(seq, np.ndarray):
        return np.ascontiguousarray(seq, dtype=np.uint8).tobytes()
    return bytes(seq)


class FastaWriter(Writer[Fasta | BgzipFasta]):
    """Writes FASTA, building its ``.fai`` (and, compressed, its ``.gzi``)
    as it goes; :meth:`close` saves them beside the file and returns them.

    ``compress`` writes bgzip FASTA at libdeflate ``level`` (0 to 12); by
    default a path ending in ``.gz`` is compressed. ``line_width`` bases go
    on each line (0: one line per sequence).
    """

    def __init__(
        self,
        path: PathLike,
        *,
        line_width: int = 60,
        compress: bool | None = None,
        level: int = 6,
    ) -> None:
        if compress is None:
            compress = os.fspath(path).endswith(".gz")
        native = _native.FastaWriter(os.fspath(path), line_width, level if compress else None)
        super().__init__(BgzipFasta(path) if compress else Fasta(path), native)

    def add(self, name: str, seq: str | bytes | np.ndarray) -> None:
        """Appends sequence ``name``; names are unique and bases printable ASCII."""
        self._native.add(name, _as_bytes(seq))

    def _finish(self) -> tuple[File, ...]:
        fai, *gzi = self._native.close()
        return (Fai(fai), *(Gzi(p) for p in gzi))
