"""Variants: the genotypes of a region, and variants applied to a reference
(which alleles to take, and the consensus sequences that result)."""

from __future__ import annotations

import collections.abc as abc
from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from .. import _native
from ._genome import Region
from ._string import StringVector

__all__ = [
    "FLOAT_END",
    "FLOAT_MISSING",
    "GT_ABSENT",
    "GT_MISSING",
    "INT_END",
    "INT_MISSING",
    "AllAlt",
    "Alleles",
    "ConsensusBatch",
    "ConsensusRecords",
    "ConsensusSequence",
    "Genotypes",
    "Haplotype",
    "Record",
    "Samples",
    "Variant",
    "Variants",
]

GT_MISSING = -1
"""A missing allele (``.``) in :attr:`Genotypes.gt`."""
GT_ABSENT = -2
"""No allele in :attr:`Genotypes.gt`: past the end of a sample's alleles
(lower ploidy than the matrix), or a record without genotypes."""


@dataclass(frozen=True, slots=True, eq=False)
class Genotypes:
    """The records overlapping a region, in file order, and the genotypes
    of the chosen samples at them.

    ``gt[i, j, k]`` is allele ``k`` of sample ``j`` at record ``i``: an
    allele index (0 = REF, 1 = the first ALT, ...), :data:`GT_MISSING` or
    :data:`GT_ABSENT`. ``phased[i, j]`` is true when every allele after the
    first is phased; a sample of one allele takes that allele's phase bit
    as htslib decodes it.
    """

    region: Region
    samples: tuple[str, ...]
    start: np.ndarray
    """int64: 0-based position of each record's first REF base."""
    end: np.ndarray
    """int64: ``start`` plus the reference span (REF, INFO/END, ...)."""
    ids: StringVector
    ref: StringVector
    alt: StringVector
    """The ALT alleles, comma-separated; empty when there are none."""
    gt: np.ndarray
    """int8, ``[records, samples, ploidy]``."""
    phased: np.ndarray
    """bool, ``[records, samples]``."""

    def __len__(self) -> int:
        return len(self.start)

    @property
    def ploidy(self) -> int:
        """Alleles per sample: the most any record of the region gives."""
        return self.gt.shape[2]


INT_MISSING = -(2**31)
"""A missing integer (``.``) in an INFO or FORMAT column of :class:`Variants`."""
INT_END = -(2**31) + 1
"""Padding past a record's (or sample's) values in an integer column."""
FLOAT_MISSING = 0x7F80_0001
"""The bits (``column.view(np.uint32)``) of a missing float; it reads as NaN."""
FLOAT_END = 0x7F80_0002
"""The bits of the padding in a float column; it reads as NaN."""


@dataclass(frozen=True, slots=True, eq=False)
class Variants:
    """The records overlapping a region, in file order (those a filter
    expression kept): :attr:`genotypes` holds their positions, IDs, alleles
    and the genotypes of the chosen samples.

    ``qual`` is ``float32``, NaN for ``.``; ``filters`` the FILTER names
    joined by ``;`` (empty for ``.``). ``info`` and ``format`` map each
    field asked for to its values:

    - INFO Flag: ``bool [records]``
    - INFO Integer / Float: ``int32`` / ``float32 [records, width]``
    - FORMAT Integer / Float: ``[records, samples, width]``
    - String: a :class:`StringVector` of a value per record (INFO) or per
      record and sample, row ``record * samples + sample`` (FORMAT)

    ``width`` is the most values any record (or sample) of the region has.
    Missing values are :data:`INT_MISSING` (floats: NaN with the bits
    :data:`FLOAT_MISSING`); padding past a record's values, and a record
    without the field, :data:`INT_END` (NaN with :data:`FLOAT_END`).
    """

    genotypes: Genotypes
    qual: np.ndarray
    filters: StringVector
    info: abc.Mapping[str, np.ndarray | StringVector]
    format: abc.Mapping[str, np.ndarray | StringVector]

    @property
    def region(self) -> Region:
        return self.genotypes.region

    def __len__(self) -> int:
        return len(self.genotypes)


@dataclass(frozen=True, slots=True)
class AllAlt:
    """The first ALT of every record, samples ignored (``bcftools consensus
    -s -``); with IUPAC codes, the codes of REF and ALT."""

    kind: ClassVar[str] = "all"


@dataclass(frozen=True, slots=True, init=False)
class Samples:
    """IUPAC codes of these samples' genotypes (``-s a,b`` without ``-H``)."""

    kind: ClassVar[str] = "samples"
    names: tuple[str, ...]

    def __init__(self, names: str | abc.Iterable[str]) -> None:
        object.__setattr__(self, "names", (names,) if isinstance(names, str) else tuple(names))


@dataclass(frozen=True, slots=True)
class Haplotype:
    """One haplotype of one sample (``-s sample -H haplotype``): ``1``, ``2``,
    ..., ``R``, ``A``, ``LR``, ``LA``, ``SR``, ``SA``, ``I`` or ``<N>pIu``."""

    kind: ClassVar[str] = "haplotype"
    sample: str
    haplotype: str = "1"

    def __post_init__(self) -> None:
        haplotype = str(self.haplotype)
        if not _native.is_haplotype(haplotype):
            raise ValueError(
                f"{haplotype!r} is not a haplotype: use 1, 2, ..., R, A, LR, LA, SR, SA, I "
                "or <N>pIu"
            )
        object.__setattr__(self, "haplotype", haplotype)


@dataclass(frozen=True, slots=True)
class Record:
    """One record of the region with its ALT ``alt`` (1 is the first),
    every other record left out: the effect of that one variant. ``index``
    counts the records the region takes in file order, as
    :meth:`~omicsweave.query.ConsensusReader.records` lists them."""

    kind: ClassVar[str] = "record"
    index: int
    alt: int = 1

    def __post_init__(self) -> None:
        if self.index < 0:
            raise ValueError(f"record index {self.index} is negative")
        if self.alt < 1:
            raise ValueError(f"alt {self.alt} is not an ALT: they are numbered from 1")


@dataclass(frozen=True, slots=True)
class Variant:
    """The record of the region at 0-based ``start`` with REF ``ref`` and
    ALT ``alt`` (one of its ALTs), every other record left out: the effect
    of that one variant. Alleles compare ignoring case; the first matching
    record in file order is taken, and a region without one is an error."""

    kind: ClassVar[str] = "variant"
    start: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start {self.start} is negative")
        if not self.ref or not self.alt:
            raise ValueError("ref and alt must not be empty")


Alleles = AllAlt | Samples | Haplotype | Record | Variant
"""Which alleles a consensus sequence takes: those of a genome
(:class:`AllAlt`, :class:`Samples`, :class:`Haplotype`), every record they
choose applied; or one variant (:class:`Record`, :class:`Variant`), every
other record left out."""


@dataclass(frozen=True, slots=True, eq=False)
class ConsensusRecords:
    """The records a consensus region takes, in file order: the order
    :class:`Record` counts them in."""

    region: Region
    start: np.ndarray
    """int64: 0-based position of each record's first REF base."""
    end: np.ndarray
    """int64: ``start`` plus the reference span (REF, INFO/END, ...)."""
    ref: StringVector
    alt: StringVector
    """The ALT alleles, comma-separated; empty when there are none."""
    masked: np.ndarray
    """bool: left out by the mask, whatever the selection."""

    def __len__(self) -> int:
        return len(self.start)


@dataclass(frozen=True, slots=True, eq=False)
class ConsensusSequence:
    """The consensus of ``region`` with ``alleles`` (``uint8`` ASCII bases),
    and its chain from reference to consensus when one was asked for."""

    region: Region
    alleles: Alleles
    bases: np.ndarray
    chain: str | None = None

    def __len__(self) -> int:
        return len(self.bases)

    @property
    def text(self) -> str:
        return self.bases.tobytes().decode("ascii")


@dataclass(frozen=True, slots=True, eq=False)
class ConsensusBatch:
    """Consensus sequences of one region, one per entry of ``alleles``,
    packed: sequence ``i`` is ``bases[offsets[i]:offsets[i + 1]]``."""

    region: Region
    alleles: tuple[Alleles, ...]
    offsets: np.ndarray
    bases: np.ndarray
    chains: tuple[str | None, ...]

    def __len__(self) -> int:
        return len(self.alleles)

    def __getitem__(self, i: int) -> ConsensusSequence:
        n = len(self)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError("batch index out of range")
        bases = self.bases[self.offsets[i] : self.offsets[i + 1]]
        return ConsensusSequence(self.region, self.alleles[i], bases, self.chains[i])

    def __iter__(self) -> abc.Iterator[ConsensusSequence]:
        return (self[i] for i in range(len(self)))
