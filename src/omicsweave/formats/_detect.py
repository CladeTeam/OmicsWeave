"""Every format by name, and telling the format of a file from its content."""

from __future__ import annotations

import os
from collections.abc import Mapping
from types import MappingProxyType

from .. import _native
from ._bai import Bai
from ._bam import Bam
from ._base import File, PathLike
from ._bcf import Bcf
from ._bed import Bed, GzipBed
from ._bigbed import BigBed
from ._bigwig import BigWig
from ._cool import Cool, Mcool
from ._csi import Csi
from ._fai import Fai
from ._fasta import BgzipFasta, Fasta, GzipFasta
from ._gzi import Gzi
from ._hic import Hic
from ._sam import GzipSam, Sam
from ._sidecar import Sidecar
from ._twobit import TwoBit
from ._vcf import GzipVcf, Vcf

__all__ = ["FORMATS", "detect"]

FORMATS: Mapping[str, type[File]] = MappingProxyType(
    {
        cls.format: cls
        for cls in (
            Fasta,
            GzipFasta,
            BgzipFasta,
            TwoBit,
            Fai,
            Gzi,
            Vcf,
            GzipVcf,
            Bcf,
            Sam,
            GzipSam,
            Bam,
            Csi,
            Bai,
            Bed,
            GzipBed,
            BigWig,
            BigBed,
            Cool,
            Mcool,
            Hic,
            Sidecar,
        )
    }
)
"""Every format class by its name (:attr:`File.format`)."""


def detect(path: PathLike) -> File:
    """The file at ``path``, as the format its content shows.

    Binary formats are told by their magic numbers (BGZF ones once
    inflated). Text, inflated first when gzip compressed, is VCF when it
    starts with ``##fileformat=VCF``, SAM with a header line (``@HD``,
    ``@SQ``, ...), FASTA with ``>``, and BED otherwise. Names play no part,
    so ``.fai`` and ``.gzi`` files, which have no signature, are not told
    apart.
    """
    return FORMATS[_native.detect(os.fspath(path))](path)
