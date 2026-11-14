"""BAM (``.bam``): SAM records encoded in binary, in BGZF."""

from __future__ import annotations

from ._base import File

__all__ = ["Bam"]


class Bam(File):
    """A BAM file."""

    __slots__ = ()
    format = "bam"
