"""Variants: the genotypes of a region, and variants applied to a reference
(which alleles to take, and the consensus sequences that result)."""

from __future__ import annotations

import collections.abc as abc
from dataclasses import dataclass

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
    "ConsensusSequence",
    "Genotypes",
    "Haplotype",
    "Samples",
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


@dataclass(frozen=True, slots=True, init=False)
class Samples:
    """IUPAC codes of these samples' genotypes (``-s a,b`` without ``-H``)."""

    names: tuple[str, ...]

    def __init__(self, names: str | abc.Iterable[str]) -> None:
        object.__setattr__(self, "names", (names,) if isinstance(names, str) else tuple(names))


@dataclass(frozen=True, slots=True)
class Haplotype:
    """One haplotype of one sample (``-s sample -H haplotype``): ``1``, ``2``,
    ..., ``R``, ``A``, ``LR``, ``LA``, ``SR``, ``SA``, ``I`` or ``<N>pIu``."""

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


Alleles = AllAlt | Samples | Haplotype
"""Which alleles a consensus sequence takes."""


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
