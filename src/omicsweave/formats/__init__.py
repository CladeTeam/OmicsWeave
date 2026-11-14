"""File formats: one class per format, each naming a file of that format.

A file object opens nothing; readers in :mod:`omicsweave.query` take it.
Formats stored more than one way have one class per storage, and index and
companion files are formats in their own right:

- FASTA: :class:`Fasta` (``.fa``), :class:`BgzipFasta` (``.fa.gz``),
  :class:`GzipFasta` (plain gzip); :class:`TwoBit` (``.2bit``)
- Indexes: :class:`Fai` (``.fai``), :class:`Gzi` (``.gzi``)
- BBI: :class:`BigWig` (``.bw``), :class:`BigBed` (``.bb``)
- BED: :class:`Bed` (``.bed``), :class:`GzipBed` (``.bed.gz``)
- Contacts: :class:`Hic` (``.hic``), :class:`Cool` (``.cool``),
  :class:`Mcool` (``.mcool``)
- Variants: :class:`Vcf` (``.vcf``), :class:`GzipVcf` (``.vcf.gz``),
  :class:`Bcf` (``.bcf``)
- Alignments: :class:`Sam` (``.sam``), :class:`GzipSam` (``.sam.gz``),
  :class:`Bam` (``.bam``)
- Binning indexes: :class:`Csi` (``.csi``), :class:`Bai` (``.bai``)
- Loader metadata: :class:`Sidecar` (``.sidecar``)

Writers produce files of a format with their companions:
:class:`FastaWriter`, :class:`BigWigWriter`, :class:`BigBedWriter`,
:class:`HicWriter`, :class:`HicBuilder`. :func:`detect` tells a file's
format from its content; :data:`FORMATS` maps every format name to its class.
"""

from ._bai import Bai
from ._bam import Bam
from ._base import File, PathLike, Writer
from ._bbi import BbiOptions
from ._bcf import Bcf
from ._bed import Bed, GzipBed
from ._bigbed import BigBed, BigBedSchema, BigBedWriter
from ._bigwig import BigWig, BigWigSummary, BigWigWriter
from ._cool import Cool, Mcool
from ._csi import Csi
from ._detect import FORMATS, detect
from ._fai import Fai
from ._fasta import BgzipFasta, Fasta, FastaWriter, GzipFasta
from ._gzi import Gzi
from ._hic import Hic, HicBuilder, HicHeader, HicWriter
from ._sam import GzipSam, Sam
from ._sidecar import Sidecar
from ._twobit import TwoBit
from ._vcf import GzipVcf, Vcf

__all__ = [
    "FORMATS",
    "Bai",
    "Bam",
    "BbiOptions",
    "Bcf",
    "Bed",
    "BgzipFasta",
    "BigBed",
    "BigBedSchema",
    "BigBedWriter",
    "BigWig",
    "BigWigSummary",
    "BigWigWriter",
    "Cool",
    "Csi",
    "Fai",
    "Fasta",
    "FastaWriter",
    "File",
    "Gzi",
    "GzipBed",
    "GzipFasta",
    "GzipSam",
    "GzipVcf",
    "Hic",
    "HicBuilder",
    "HicHeader",
    "HicWriter",
    "Mcool",
    "PathLike",
    "Sam",
    "Sidecar",
    "TwoBit",
    "Vcf",
    "Writer",
    "detect",
]
