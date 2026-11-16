"""The values every reader takes and returns, whatever file they come from.

- Genome: :class:`Region` (the key, with a strand), :class:`Spliced` (spans
  read as one sequence; with :class:`Region` a :data:`Locus`, the key of
  dense queries), :class:`Sequence`, :class:`Intervals`,
  :class:`MultipleAlignment`
- Signal: :class:`Signal`, :class:`SignalBins` (:class:`Stat`),
  :class:`SignalIntervals`
- Contacts: :class:`ContactWindow` (the key), :class:`BinUnit`,
  :class:`ContactMatrix`, :class:`ContactPixels`
- Variants: :class:`Variants` (:class:`Genotypes`); :class:`AllAlt`, :class:`Samples`,
  :class:`Haplotype`, :class:`Record`, :class:`Variant` (:data:`Alleles`),
  :class:`ConsensusSequence`, :class:`ConsensusBatch`, :class:`ConsensusRecords`
- Alignments: :class:`Alignments`, :class:`AlignmentField`, :class:`Cigars`,
  :class:`Modifications` (:class:`ModificationFlag`), :class:`SamFlag`
- Motifs: :class:`Alphabet`, :class:`Motif` (:class:`Factor`,
  :class:`Gaps`)
- Strings: :class:`StringVector`

Results carry their key; their data are numpy arrays, in columns.
"""

from ._alignment import (
    AlignmentField,
    Alignments,
    Cigars,
    ModificationFlag,
    Modifications,
    SamFlag,
)
from ._contact import BinUnit, ContactMatrix, ContactPixels, ContactWindow
from ._genome import (
    STRANDS,
    Intervals,
    Locus,
    LocusLike,
    MultipleAlignment,
    Region,
    RegionLike,
    Sequence,
    Spliced,
)
from ._motif import Alphabet, Factor, Gaps, Motif
from ._signal import Signal, SignalBins, SignalIntervals, Stat
from ._string import StringVector
from ._variant import (
    FLOAT_END,
    FLOAT_MISSING,
    GT_ABSENT,
    GT_MISSING,
    INT_END,
    INT_MISSING,
    AllAlt,
    Alleles,
    ConsensusBatch,
    ConsensusRecords,
    ConsensusSequence,
    Genotypes,
    Haplotype,
    Record,
    Samples,
    Variant,
    Variants,
)

__all__ = [
    "FLOAT_END",
    "FLOAT_MISSING",
    "GT_ABSENT",
    "GT_MISSING",
    "INT_END",
    "INT_MISSING",
    "STRANDS",
    "AlignmentField",
    "Alignments",
    "AllAlt",
    "Alleles",
    "Alphabet",
    "BinUnit",
    "Cigars",
    "ConsensusBatch",
    "ConsensusRecords",
    "ConsensusSequence",
    "ContactMatrix",
    "ContactPixels",
    "ContactWindow",
    "Factor",
    "Gaps",
    "Genotypes",
    "Haplotype",
    "Intervals",
    "Locus",
    "LocusLike",
    "MultipleAlignment",
    "ModificationFlag",
    "Modifications",
    "Motif",
    "Record",
    "Region",
    "RegionLike",
    "SamFlag",
    "Samples",
    "Sequence",
    "Signal",
    "SignalBins",
    "SignalIntervals",
    "Spliced",
    "Stat",
    "StringVector",
    "Variant",
    "Variants",
]
