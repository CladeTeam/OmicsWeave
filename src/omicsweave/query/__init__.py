"""Readers: open files of a format, look entries up by key, return
:mod:`omicsweave.types` values.

- :class:`SequenceReader`: FASTA, bgzip FASTA, 2bit -> :class:`~omicsweave.types.Sequence`
- :class:`TrackReader`: bigWig -> :class:`~omicsweave.types.Signal`,
  :class:`~omicsweave.types.SignalIntervals`
- :class:`IntervalReader`: bigBed, BED -> :class:`~omicsweave.types.Intervals`
- :class:`ContactReader`: ``.hic`` -> :class:`~omicsweave.types.ContactMatrix`,
  :class:`~omicsweave.types.ContactPixels`
- :class:`ConsensusReader`: VCF/BCF over a reference ->
  :class:`~omicsweave.types.ConsensusSequence`, :class:`~omicsweave.types.ConsensusBatch`

Each takes a :mod:`omicsweave.formats` file, or a path whose format is
detected; :func:`open` picks the reader for a file.
"""

from ._base import Reader
from ._consensus import ConsensusOptions, ConsensusReader
from ._contact import ContactReader
from ._interval import IntervalReader
from ._open import open
from ._sequence import SequenceReader
from ._track import TrackReader

__all__ = [
    "ConsensusOptions",
    "ConsensusReader",
    "ContactReader",
    "IntervalReader",
    "Reader",
    "SequenceReader",
    "TrackReader",
    "open",
]
