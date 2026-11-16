"""File formats: one class per format, each naming a file of that format.

A file object opens nothing; readers in :mod:`omicsweave.query` take it.
Formats stored more than one way have one class per storage, and index and
companion files are formats in their own right:

- FASTA: :class:`Fasta` (``.fa``), :class:`BgzipFasta` (``.fa.gz``),
  :class:`GzipFasta` (plain gzip); :class:`TwoBit` (``.2bit``)
- Indexes: :class:`Fai` (``.fai``), :class:`Gzi` (``.gzi``)
- BBI: :class:`BigWig` (``.bw``), :class:`BigBed` (``.bb``)
- BED: :class:`Bed` (``.bed``), :class:`GzipBed` (``.bed.gz``)
- Methylation: :class:`BedMethyl` (``.bedmethyl``), :class:`GzipBedMethyl`
  (``.bedmethyl.gz``)
- Multiple alignments: :class:`Maf` (``.maf``), :class:`GzipMaf`
  (``.maf.gz``); bigMaf is a :class:`BigBed`
- Signal text: :class:`BedGraph` (``.bedgraph``), :class:`GzipBedGraph`
  (``.bedgraph.gz``), :class:`Wig` (``.wig``), :class:`GzipWig`
  (``.wig.gz``)
- Contacts: :class:`Hic` (``.hic``), :class:`Cool` (``.cool``),
  :class:`Mcool` (``.mcool``); pairs text: :class:`Pairs` (``.pairs``),
  :class:`GzipPairs` (``.pairs.gz``), :class:`JuicerShort`,
  :class:`GzipJuicerShort`, :class:`HicProPairs` (``.allValidPairs`` /
  ``.validPairs``), :class:`GzipHicProPairs`
- Variants: :class:`Vcf` (``.vcf``), :class:`GzipVcf` (``.vcf.gz``),
  :class:`Bcf` (``.bcf``)
- Alignments: :class:`Sam` (``.sam``), :class:`GzipSam` (``.sam.gz``),
  :class:`Bam` (``.bam``), :class:`Cram` (``.cram``)
- Binning indexes: :class:`Csi` (``.csi``), :class:`Bai` (``.bai``);
  CRAM's :class:`Crai` (``.crai``)
- Loader metadata: :class:`Sidecar` (``.sidecar``)
- Motifs: :class:`Ffm` (``.ffm``); motif texts that convert into it:
  :class:`Meme`, :class:`Jaspar`, :class:`JasparPfm`, :class:`Transfac`,
  :class:`Homer`, :class:`CisBp`, :class:`Pcm` (HOCOMOCO PCM / PFM,
  Cluster-Buster), :class:`Pwm`, :class:`Dpcm`, :class:`Dpwm`,
  :class:`Bamm`, :class:`ProBound`, :class:`Glam2`

Writers produce files of a format with their companions:
:class:`FastaWriter`, :class:`TwoBitWriter`, :class:`BigWigWriter`,
:class:`BigBedWriter`, :class:`HicWriter`, :class:`HicBuilder`,
:class:`FfmWriter`;
:func:`read_sidecar` reads what a sidecar says of a file. :func:`detect` tells a file's
format from its content; :data:`FORMATS` maps every format name to its class.
"""

from ._bai import Bai
from ._bam import Bam
from ._base import File, PathLike, Writer
from ._bbi import BbiOptions
from ._bcf import Bcf
from ._bed import Bed, GzipBed
from ._bedgraph import BedGraph, GzipBedGraph
from ._bedmethyl import BedMethyl, GzipBedMethyl
from ._bigbed import AutoSqlField, BigBed, BigBedSchema, BigBedWriter
from ._bigwig import BigWig, BigWigSummary, BigWigWriter
from ._cool import Cool, Mcool
from ._crai import Crai
from ._cram import Cram
from ._csi import Csi
from ._detect import FORMATS, detect
from ._fai import Fai
from ._fasta import BgzipFasta, Fasta, FastaWriter, GzipFasta
from ._gzi import Gzi
from ._hic import Hic, HicBlocks, HicBuilder, HicHeader, HicWriter
from ._maf import GzipMaf, Maf
from ._motif import (
    Bamm,
    CisBp,
    Dpcm,
    Dpwm,
    Ffm,
    FfmWriter,
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
from ._sidecar import (
    BbiLevelBlocks,
    BbiSidecar,
    HicMapBlocks,
    HicSidecar,
    Sidecar,
    read_sidecar,
)
from ._twobit import TwoBit, TwoBitWriter
from ._vcf import GzipVcf, Vcf
from ._wig import GzipWig, Wig

__all__ = [
    "AutoSqlField",
    "Bai",
    "Bam",
    "Bamm",
    "BbiLevelBlocks",
    "BbiOptions",
    "BbiSidecar",
    "Bcf",
    "Bed",
    "BedGraph",
    "BedMethyl",
    "BgzipFasta",
    "BigBed",
    "BigBedSchema",
    "BigBedWriter",
    "BigWig",
    "BigWigSummary",
    "BigWigWriter",
    "CisBp",
    "Cool",
    "Crai",
    "Cram",
    "Csi",
    "Dpcm",
    "Dpwm",
    "FORMATS",
    "Fai",
    "Fasta",
    "FastaWriter",
    "Ffm",
    "FfmWriter",
    "File",
    "Glam2",
    "Gzi",
    "GzipBed",
    "GzipBedGraph",
    "GzipBedMethyl",
    "GzipFasta",
    "GzipJuicerShort",
    "GzipHicProPairs",
    "GzipMaf",
    "GzipPairs",
    "GzipSam",
    "GzipVcf",
    "GzipWig",
    "Hic",
    "HicBlocks",
    "HicBuilder",
    "HicHeader",
    "HicMapBlocks",
    "HicSidecar",
    "HicWriter",
    "Homer",
    "Jaspar",
    "JasparPfm",
    "JuicerShort",
    "HicProPairs",
    "Maf",
    "Mcool",
    "Meme",
    "Pairs",
    "PathLike",
    "Pcm",
    "ProBound",
    "Pwm",
    "Sam",
    "Sidecar",
    "Transfac",
    "TwoBit",
    "TwoBitWriter",
    "Vcf",
    "Wig",
    "Writer",
    "detect",
    "read_sidecar",
]
