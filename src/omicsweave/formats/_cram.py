"""CRAM (``.cram``): SAM records compressed by data series, against a
reference."""

from __future__ import annotations

import os

from .. import _native
from ._base import File, ensure_format
from ._crai import Crai

__all__ = ["Cram"]


class Cram(File):
    """A CRAM file (versions 3.0 and 3.1). Read by region through its
    ``.crai`` (:class:`~omicsweave.query.AlignmentReader`, with the
    reference its slices were compressed against)."""

    __slots__ = ()
    format = "cram"

    def build_index(self) -> Crai:
        """Builds the ``.crai`` of the coordinate-sorted file beside it
        (``<cram>.crai``), as ``samtools index`` does. Only multi-reference
        slices are decoded; no reference is needed."""
        ensure_format(self)
        return Crai(_native.build_cram_index(os.fspath(self)))
