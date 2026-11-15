"""The values every reader takes and returns, whatever file they come from.

- Genome: :class:`Region` (the key, with a strand), :class:`Spliced` (spans
  read as one sequence; with :class:`Region` a :data:`Locus`, the key of
  dense queries), :class:`Sequence`, :class:`Intervals`
- Signal: :class:`Signal`, :class:`SignalBins` (:class:`Stat`),
  :class:`SignalIntervals`
- Contacts: :class:`ContactWindow` (the key), :class:`BinUnit`,
  :class:`ContactMatrix`, :class:`ContactPixels`
- Variants: :class:`Genotypes`; :class:`AllAlt`, :class:`Samples`,
  :class:`Haplotype` (:data:`Alleles`), :class:`ConsensusSequence`,
  :class:`ConsensusBatch`
- Alignments: :class:`Alignments`, :class:`AlignmentField`, :class:`Cigars`
- Strings: :class:`StringVector`

Results carry their key; their data are numpy arrays, in columns.
"""

from ._alignment import AlignmentField, Alignments, Cigars
from ._contact import BinUnit, ContactMatrix, ContactPixels, ContactWindow
from ._genome import STRANDS, Intervals, Locus, LocusLike, Region, RegionLike, Sequence, Spliced
from ._signal import Signal, SignalBins, SignalIntervals, Stat
from ._string import StringVector
from ._variant import (
    GT_ABSENT,
    GT_MISSING,
    AllAlt,
    Alleles,
    ConsensusBatch,
    ConsensusSequence,
    Genotypes,
    Haplotype,
    Samples,
)

__all__ = [
    "STRANDS",
    "GT_ABSENT",
    "GT_MISSING",
    "AlignmentField",
    "Alignments",
    "AllAlt",
    "Alleles",
    "BinUnit",
    "Cigars",
    "ConsensusBatch",
    "ConsensusSequence",
    "ContactMatrix",
    "ContactPixels",
    "ContactWindow",
    "Genotypes",
    "Haplotype",
    "Intervals",
    "Locus",
    "LocusLike",
    "Region",
    "RegionLike",
    "Samples",
    "Sequence",
    "Signal",
    "SignalBins",
    "SignalIntervals",
    "Spliced",
    "Stat",
    "StringVector",
]
