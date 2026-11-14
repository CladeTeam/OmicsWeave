"""Consensus sequences by region: ``bcftools consensus`` of a VCF or BCF
over a reference, computed in memory."""

from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from .. import _native
from ..formats import Bcf, BgzipFasta, Fasta, GzipVcf, PathLike, TwoBit, Vcf
from ..types import (
    AllAlt,
    Alleles,
    ConsensusBatch,
    ConsensusSequence,
    Haplotype,
    RegionLike,
    Samples,
)
from ..types._genome import as_region
from ._base import Reader, resolve
from ._sequence import SequenceReader

__all__ = ["ConsensusOptions", "ConsensusReader"]

Mark = Literal["upper", "lower"] | str


@dataclass(frozen=True, slots=True)
class ConsensusOptions:
    """The ``bcftools consensus`` options that shape every sequence."""

    iupac_codes: bool = False
    """``-I``: IUPAC codes instead of one allele."""
    missing: str | None = None
    """``-M``: this character for a missing genotype."""
    absent: str | None = None
    """``-a``: this character where no record covers the reference."""
    mark_del: str | None = None
    """``--mark-del``: this character for each deleted base."""
    mark_ins: Mark | None = None
    """``--mark-ins``: ``upper``, ``lower`` or a character."""
    mark_snv: Mark | None = None
    """``--mark-snv``: ``upper``, ``lower`` or a character."""
    mask: PathLike | None = None
    """``-m``: a BED file of regions to mask."""
    mask_with: Mark = "N"
    """``--mask-with``: ``upper``, ``lower`` or a character."""
    regions_overlap: Literal["start", "record", "variant"] = "record"
    """``--regions-overlap``: which records a region takes."""
    chain: bool = False
    """``-c``: also produce the chain from reference to consensus."""
    include: str | None = None
    """``-i``: only records passing this expression."""
    exclude: str | None = None
    """``-e``: leave out records passing this expression."""

    def __post_init__(self) -> None:
        if self.include is not None and self.exclude is not None:
            raise ValueError("give include or exclude, not both")


def _selector(alleles: Alleles) -> tuple[str, list[str], str | None]:
    match alleles:
        case AllAlt():
            return ("all", [], None)
        case Samples(names):
            return ("samples", list(names), None)
        case Haplotype(sample, haplotype):
            return ("haplotype", [sample], haplotype)
    raise TypeError(f"expected AllAlt, Samples or Haplotype, not {alleles!r}")


class ConsensusReader(Reader):
    """Reads personal sequences by region: the records of a VCF or BCF file
    (loaded once, whole) applied to a reference.

    ``reference`` is an open :class:`SequenceReader`, shared and left open,
    or a sequence file, which this reader opens and closes. Consecutive
    reads of one region read its reference and find its records once;
    :meth:`fetch_batch` builds a whole batch in one call.
    """

    def __init__(
        self,
        variants: Vcf | GzipVcf | Bcf | PathLike,
        reference: SequenceReader | Fasta | BgzipFasta | TwoBit | PathLike,
        options: ConsensusOptions | None = None,
    ) -> None:
        file = resolve(variants, (Vcf, GzipVcf, Bcf), "ConsensusReader")
        options = options or ConsensusOptions()
        given = reference
        self._owns_reference = not isinstance(reference, SequenceReader)
        if not isinstance(reference, SequenceReader):
            reference = SequenceReader(reference)
        o = options
        try:
            native = _native.ConsensusBuilder(
                os.fspath(file),
                reference._native,
                o.iupac_codes,
                o.missing,
                o.absent,
                o.mark_del,
                o.mark_ins,
                o.mark_snv,
                None if o.mask is None else os.fspath(o.mask),
                o.mask_with,
                o.regions_overlap,
                o.chain,
                o.include,
                o.exclude,
            )
        except BaseException:
            if self._owns_reference:
                reference.close()
            raise
        self._open(native, file, given, options)
        self.source = file
        self.reference = reference
        self.options = options
        self._samples = tuple(native.samples())

    @property
    def samples(self) -> tuple[str, ...]:
        return self._samples

    @property
    def records(self) -> int:
        """Records kept (after ``include`` / ``exclude``)."""
        return self._native.records()[0]

    @property
    def skipped_lines(self) -> int:
        """Malformed VCF lines left out while reading, as htslib skips them."""
        return self._native.records()[1]

    def fetch(self, region: RegionLike, alleles: Alleles | None = None) -> ConsensusSequence:
        """The consensus of ``region`` with ``alleles`` (default :class:`AllAlt`)."""
        r = as_region(region)
        alleles = AllAlt() if alleles is None else alleles
        bases, chain = self._native.sequence(r.chrom, r.start, r.end, _selector(alleles))
        return ConsensusSequence(r, alleles, bases, chain)

    def fetch_batch(self, region: RegionLike, alleles: Iterable[Alleles]) -> ConsensusBatch:
        """The consensus of ``region`` with each of ``alleles``."""
        r = as_region(region)
        alleles = tuple(alleles)
        offsets, bases, chains = self._native.sequences(
            r.chrom, r.start, r.end, [_selector(a) for a in alleles]
        )
        return ConsensusBatch(r, alleles, offsets, bases, tuple(chains))

    def close(self) -> None:
        """Closes the reader, and its reference if it opened it."""
        super().close()
        if self._owns_reference:
            self.reference.close()
