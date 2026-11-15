"""bedGraph, stored two ways: plain text (``.bedgraph``) and gzip
(``.bedgraph.gz``).

A bedGraph is BED with four fields (``chrom start end value``), so it is
also read as BED and converts to bigBed as well as to bigWig.
"""

from __future__ import annotations

from ._base import File

__all__ = ["BedGraph", "GzipBedGraph"]


class BedGraph(File):
    """Plain-text bedGraph (``.bedgraph``, ``.bdg``, ``.bg``)."""

    __slots__ = ()
    format = "bedgraph"


class GzipBedGraph(File):
    """bedGraph compressed with gzip or bgzip (``.bedgraph.gz``)."""

    __slots__ = ()
    format = "gzip-bedgraph"
