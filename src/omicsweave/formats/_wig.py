"""Wiggle, stored two ways: plain text (``.wig``) and gzip (``.wig.gz``).

``fixedStep``, ``variableStep`` and bedGraph lines, as UCSC ``wigToBigWig``
reads them; wig files convert to bigWig.
"""

from __future__ import annotations

from ._base import File

__all__ = ["GzipWig", "Wig"]


class Wig(File):
    """Plain-text wiggle (``.wig``)."""

    __slots__ = ()
    format = "wig"


class GzipWig(File):
    """Wiggle compressed with gzip or bgzip (``.wig.gz``)."""

    __slots__ = ()
    format = "gzip-wig"
