"""BCF (``.bcf``): VCF records encoded in binary, in BGZF."""

from __future__ import annotations

import os

from .. import _native
from .._budget import threads as _threads
from ._base import File, ensure_format
from ._csi import Csi

__all__ = ["Bcf"]


class Bcf(File):
    """A BCF file (version 2). Read by region through its CSI index
    (:class:`~omicsweave.query.VariantReader`)."""

    __slots__ = ()
    format = "bcf"

    def build_index(self, *, min_shift: int = 14, threads: int | None = None) -> Csi:
        """Builds the ``.csi`` of the sorted, BGZF-compressed file beside it
        (``<bcf>.csi``), as ``bcftools index`` does: bins of
        ``2**min_shift`` bases, as many levels as the longest contig needs.
        ``threads`` (default: every core the process may run on) inflate
        the file. Unsorted records raise
        :class:`~omicsweave.errors.InvalidError`."""
        ensure_format(self)
        return Csi(_native.build_bcf_index(os.fspath(self), min_shift, _threads(threads)))
