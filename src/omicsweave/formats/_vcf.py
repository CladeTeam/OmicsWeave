"""VCF, stored two ways: plain text (``.vcf``) and gzip (``.vcf.gz``).

Either is read whole; any gzip stream will do, BGZF included.
"""

from __future__ import annotations

from ._base import File

__all__ = ["GzipVcf", "Vcf"]


class Vcf(File):
    """Plain-text VCF (``.vcf``)."""

    __slots__ = ()
    format = "vcf"


class GzipVcf(File):
    """VCF compressed with gzip or bgzip (``.vcf.gz``)."""

    __slots__ = ()
    format = "gzip-vcf"
