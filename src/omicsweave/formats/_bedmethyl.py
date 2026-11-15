"""bedMethyl, stored two ways: plain text (``.bedmethyl``) and gzip
(``.bedmethyl.gz``).

Two dialects are read: modkit ``pileup`` output (BED9 and nine counts;
older modkit separates the counts with spaces) and ENCODE WGBS bedMethyl
(BED9, coverage and percentage). bedMethyl converts to bigBed with every
column and to bigWig with one.
"""

from __future__ import annotations

from ._base import File

__all__ = ["BedMethyl", "GzipBedMethyl"]


class BedMethyl(File):
    """Plain-text bedMethyl (``.bedmethyl``, or ``.bed`` told by content)."""

    __slots__ = ()
    format = "bedmethyl"


class GzipBedMethyl(File):
    """bedMethyl compressed with gzip or bgzip (``.bedmethyl.gz``)."""

    __slots__ = ()
    format = "gzip-bedmethyl"
