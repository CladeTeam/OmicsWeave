"""MAF multiple alignments, stored two ways: plain text (``.maf``) and gzip
(``.maf.gz``).

MAF is not read by region: it converts to bigMaf, a bigBed holding one
record per alignment block (``ow convert hg38.100way.maf.gz
hg38.100way.bigmaf``), which :class:`~omicsweave.query.MsaReader` reads.
"""

from __future__ import annotations

from ._base import File

__all__ = ["GzipMaf", "Maf"]


class Maf(File):
    """Plain-text MAF (``.maf``): a ``##maf`` line, then ``a`` blocks of
    ``s``, ``e``, ``i``, ``q`` lines separated by blank lines."""

    __slots__ = ()
    format = "maf"


class GzipMaf(File):
    """MAF compressed with gzip or bgzip (``.maf.gz``)."""

    __slots__ = ()
    format = "gzip-maf"
