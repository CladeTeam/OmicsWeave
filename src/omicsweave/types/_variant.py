"""Variants applied to a reference: which alleles to take, and the
consensus sequences that result."""

from __future__ import annotations

import collections.abc as abc
from dataclasses import dataclass

import numpy as np

from ._genome import Region

__all__ = ["AllAlt", "Alleles", "ConsensusBatch", "ConsensusSequence", "Haplotype", "Samples"]


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
        object.__setattr__(self, "haplotype", str(self.haplotype))


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
