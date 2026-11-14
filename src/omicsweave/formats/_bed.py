"""BED, stored two ways: plain text (``.bed``) and gzip (``.bed.gz``).

Either is read whole; any gzip stream will do, BGZF included.
"""

from __future__ import annotations

from ._base import File

__all__ = ["Bed", "GzipBed"]


class Bed(File):
    """Plain-text BED (``.bed``)."""

    __slots__ = ()
    format = "bed"


class GzipBed(File):
    """BED compressed with gzip or bgzip (``.bed.gz``)."""

    __slots__ = ()
    format = "gzip-bed"
