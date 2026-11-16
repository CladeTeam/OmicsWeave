"""BAM (``.bam``): SAM records encoded in binary, in BGZF."""

from __future__ import annotations

import os
from typing import Literal

from .. import _native
from .._budget import threads as _threads
from ._bai import Bai
from ._base import File, ensure_format
from ._csi import Csi

__all__ = ["Bam"]


class Bam(File):
    """A BAM file. Read by region through its index
    (:class:`~omicsweave.query.AlignmentReader`)."""

    __slots__ = ()
    format = "bam"

    def build_index(
        self,
        kind: Literal["auto", "bai", "csi"] = "auto",
        *,
        min_shift: int = 14,
        threads: int | None = None,
    ) -> Bai | Csi:
        """Builds the index of the coordinate-sorted file beside it, as
        ``samtools index`` does: ``<bam>.bai``, or ``<bam>.csi`` with bins
        of ``2**min_shift`` bases. ``auto`` builds a BAI when every
        reference fits one (at most 2**29 bases), else a CSI. ``threads``
        (default: every core the process may run on) inflate the file."""
        ensure_format(self)
        made, path = _native.build_bam_index(os.fspath(self), kind, min_shift, _threads(threads))
        return Bai(path) if made == "bai" else Csi(path)
