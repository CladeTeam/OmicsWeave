"""BCF (``.bcf``): VCF records encoded in binary, in BGZF."""

from __future__ import annotations

from ._base import File

__all__ = ["Bcf"]


class Bcf(File):
    """A BCF file (version 2). Read by region through its CSI index
    (:class:`~omicsweave.query.VariantReader`)."""

    __slots__ = ()
    format = "bcf"
