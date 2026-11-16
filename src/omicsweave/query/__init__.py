"""Readers: open files of a format, look entries up by key, return
:mod:`omicsweave.types` values.

- :class:`SequenceReader`: FASTA, bgzip FASTA, 2bit -> :class:`~omicsweave.types.Sequence`
- :class:`TrackReader`: bigWig -> :class:`~omicsweave.types.Signal`,
  :class:`~omicsweave.types.SignalBins`, :class:`~omicsweave.types.SignalIntervals`;
  :class:`StrandedTrackReader` reads a plus / minus pair by the key's strand
- :class:`IntervalReader`: bigBed -> :class:`~omicsweave.types.Intervals`
- :class:`MsaReader`: bigMaf (a bigBed) -> :class:`~omicsweave.types.MultipleAlignment`
- :class:`ContactReader`: ``.hic`` -> :class:`~omicsweave.types.ContactMatrix`,
  :class:`~omicsweave.types.ContactPixels`
- :class:`ConsensusReader`: BCF (+ CSI) over a reference ->
  :class:`~omicsweave.types.ConsensusSequence`, :class:`~omicsweave.types.ConsensusBatch`
- :class:`AlignmentReader`: BAM (+ BAI / CSI) -> :class:`~omicsweave.types.Alignments`
- :class:`VariantReader`: BCF (+ CSI) -> :class:`~omicsweave.types.Genotypes`
- :class:`MotifReader`: FFM -> :class:`~omicsweave.types.Motif` by number or ID

Each takes a :mod:`omicsweave.formats` file, or a path whose format is
detected; :func:`open` picks the reader for a file. Text formats (BED,
signal text, VCF, SAM, ...) are not read by region: :func:`open` refuses
them, naming the ``ow convert`` command that makes them readable. Dense reads (bases,
signal, bins, contact tiles) take a :data:`~omicsweave.types.Locus` and run
5' to 3' of its strand; sparse reads take a region and keep genomic order.
"""

from ._alignment import AlignmentReader
from ._base import Reader
from ._consensus import ConsensusOptions, ConsensusReader
from ._contact import ContactReader
from ._interval import IntervalReader
from ._motif import MotifReader
from ._msa import MsaReader
from ._open import open
from ._sequence import SequenceReader
from ._track import StrandedTrackReader, TrackReader
from ._variant import VariantReader

__all__ = [
    "AlignmentReader",
    "ConsensusOptions",
    "ConsensusReader",
    "ContactReader",
    "IntervalReader",
    "MotifReader",
    "MsaReader",
    "Reader",
    "SequenceReader",
    "StrandedTrackReader",
    "TrackReader",
    "VariantReader",
    "open",
]
