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
from ._bedgraph import BedGraph, GzipBedGraph
from ._bedmethyl import BedMethyl, GzipBedMethyl
from ._bigbed import BigBed
from ._bigwig import BigWig
from ._cool import Cool, Mcool
from ._crai import Crai
from ._cram import Cram
from ._csi import Csi
from ._fai import Fai
from ._fasta import BgzipFasta, Fasta, GzipFasta
from ._gzi import Gzi
from ._hic import Hic
from ._maf import GzipMaf, Maf
from ._motif import (
    Bamm,
    CisBp,
    Dpcm,
    Dpwm,
    Ffm,
    Glam2,
    Homer,
    Jaspar,
    JasparPfm,
    Meme,
    Pcm,
    ProBound,
    Pwm,
    Transfac,
)
from ._pairs import GzipHicProPairs, GzipJuicerShort, GzipPairs, HicProPairs, JuicerShort, Pairs
from ._sam import GzipSam, Sam
from ._sidecar import Sidecar
from ._twobit import TwoBit
from ._vcf import GzipVcf, Vcf
from ._wig import GzipWig, Wig

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
            Cram,
            Csi,
            Bai,
            Crai,
            Bed,
            GzipBed,
            BedGraph,
            GzipBedGraph,
            Wig,
            GzipWig,
            BigWig,
            BigBed,
            Cool,
            Mcool,
            Hic,
            Sidecar,
            Pairs,
            GzipPairs,
            JuicerShort,
            GzipJuicerShort,
            HicProPairs,
            GzipHicProPairs,
            BedMethyl,
            GzipBedMethyl,
            Maf,
            GzipMaf,
            Meme,
            Jaspar,
            JasparPfm,
            Transfac,
            Homer,
            CisBp,
            Pcm,
            Pwm,
            Dpcm,
            Dpwm,
            Bamm,
            ProBound,
            Glam2,
            Ffm,
        )
    }
)
"""Every format class by its name (:attr:`File.format`)."""


def detect(path: PathLike) -> File:
    """The file at ``path``, as the format its content shows.

    Binary formats are told by their magic numbers (BGZF ones once
    inflated; CRAM by ``CRAM``). A ``.crai`` is gzip text whose first line
    is six tab-separated integers. Text, inflated first when gzip compressed, is VCF when it
    starts with ``##fileformat=VCF``, SAM with a header line (``@HD``,
    ``@SQ``, ...), FASTA with ``>``, 4DN pairs with ``## pairs format``,
    MAF with ``##maf``. Motif texts: MEME with ``MEME version``, CIS-BP with
    a ``Pos A C G T`` header, TRANSFAC with line codes and a ``P0`` header,
    ProBound with a JSON object holding ``coefficients``, GLAM2 with its
    ``GLAM2: Gapped Local Alignment of Motifs`` banner; ``>``-headed
    numeric matrices are JASPAR (``A [ ... ]`` rows), HOMER (a numeric third
    header column), raw PFM (four rows per motif), PCM / PWM (four values
    per row; PWM when any is negative) or their dinucleotide forms (sixteen
    per row). Square 4 x 4 matrices cannot be told apart and need their
    format named. An FFM is told by its magic number.
    Other text is told from its first 64 KiB: wig when a ``track`` line
    says ``type=wiggle_0`` or a line declares a ``fixedStep`` /
    ``variableStep`` section, bedGraph when a ``track`` line says
    ``type=bedGraph`` or every data line has four fields, the last a
    number, HiC-Pro when every data line has 12 or 13 fields, +/- strands
    in columns 4 and 7 and non-negative u64 positions in columns 3 and 6,
    Juicer short when every data line is ``str1 chr1 pos1 frag1
    str2 chr2 pos2 frag2 [score]``, bedMethyl when every data line is BED9
    with modkit's nine counts or ENCODE's coverage and percentage, and BED
    otherwise. Names play no part,
    so ``.fai`` and ``.gzi`` files, which have no signature, are not told
    apart.
    """
    return FORMATS[_native.detect(os.fspath(path))](path)
