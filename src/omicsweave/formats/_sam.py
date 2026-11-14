"""SAM, stored two ways: plain text (``.sam``) and gzip (``.sam.gz``).

Either is read whole; any gzip stream will do, BGZF included.
"""

from __future__ import annotations

from ._base import File

__all__ = ["GzipSam", "Sam"]


class Sam(File):
    """Plain-text SAM (``.sam``)."""

    __slots__ = ()
    format = "sam"


class GzipSam(File):
    """SAM compressed with gzip or bgzip (``.sam.gz``)."""

    __slots__ = ()
    format = "gzip-sam"
