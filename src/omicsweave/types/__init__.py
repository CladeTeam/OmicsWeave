"""The values every reader takes and returns, whatever file they come from.

- Genome: :class:`Region` (the key), :class:`Sequence`, :class:`Intervals`
- Signal: :class:`Signal`, :class:`SignalIntervals`
- Contacts: :class:`ContactWindow` (the key), :class:`BinUnit`,
  :class:`ContactMatrix`, :class:`ContactPixels`
- Variants: :class:`AllAlt`, :class:`Samples`, :class:`Haplotype`
  (:data:`Alleles`), :class:`ConsensusSequence`, :class:`ConsensusBatch`
- Strings: :class:`StringVector`

Results carry their key; their data are numpy arrays, in columns.
"""

from ._contact import BinUnit, ContactMatrix, ContactPixels, ContactWindow
from ._genome import Intervals, Region, RegionLike, Sequence
from ._signal import Signal, SignalIntervals
from ._string import StringVector
from ._variant import AllAlt, Alleles, ConsensusBatch, ConsensusSequence, Haplotype, Samples

__all__ = [
    "AllAlt",
    "Alleles",
    "BinUnit",
    "ConsensusBatch",
    "ConsensusSequence",
    "ContactMatrix",
    "ContactPixels",
    "ContactWindow",
    "Haplotype",
    "Intervals",
    "Region",
    "RegionLike",
    "Samples",
    "Sequence",
    "Signal",
    "SignalIntervals",
    "StringVector",
]
