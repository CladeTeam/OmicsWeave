"""Consensus sequences by region: ``bcftools consensus`` of an indexed BCF
over a reference, computed in memory."""

from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from .. import _native
from ..errors import InvalidError
from ..formats import Bcf, BgzipFasta, Csi, Fasta, File, GzipVcf, PathLike, TwoBit, Vcf, detect
from ..types import (
    AllAlt,
    Alleles,
    ConsensusBatch,
    ConsensusSequence,
    Haplotype,
    RegionLike,
    Samples,
)
from ..types._genome import Region, Spliced, as_region
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


def _bcf(variants: Bcf | PathLike) -> Bcf:
    """``variants`` as a BCF: VCF text is not read by region."""
    file = variants if isinstance(variants, File) else detect(variants)
    if isinstance(file, (Vcf, GzipVcf)):
        detail = (
            "ConsensusReader reads an indexed BCF by region, not VCF text; convert it "
            f"with `ow convert {file.path} {_bcf_name(file)} --index`"
        )
        if isinstance(variants, File):
            raise TypeError(detail)
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    return resolve(variants, (Bcf,), "ConsensusReader")


def _bcf_name(file: Vcf | GzipVcf) -> str:
    name = os.fspath(file.path)
    for suffix in (".vcf.gz", ".vcf.bgz", ".vcf"):
        if name.endswith(suffix):
            return name[: -len(suffix)] + ".bcf"
    return name + ".bcf"


def _key(region: RegionLike) -> Region:
    if isinstance(region, Spliced):
        raise TypeError(
            "ConsensusReader does not read Spliced keys: a consensus edits one "
            "contiguous region; fetch each span as a Region (with its strand)"
        )
    return as_region(region)


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
    """Reads personal sequences by region: the records of a BCF, read by
    region through its CSI index, applied to a reference.

    ``index`` defaults to the ``.csi`` beside the BCF; a missing index
    raises :class:`~omicsweave.errors.IoError`. VCF text is not read by
    region: convert it first (``ow convert cohort.vcf.gz cohort.bcf
    --index``). ``reference`` is an open :class:`SequenceReader`, shared and
    left open, or a sequence file, which this reader opens and closes.
    Each region decodes the records it takes once; consecutive reads of one
    region reuse its reference and records, and :meth:`fetch_batch` builds
    a whole batch in one call.

    A region on ``"-"`` reads the reverse complement of its consensus
    (IUPAC codes complemented, case kept); the chain stays on the forward
    strand. ``"+"``, ``"."`` and ``"?"`` read forward.
    """

    def __init__(
        self,
        variants: Bcf | PathLike,
        reference: SequenceReader | Fasta | BgzipFasta | TwoBit | PathLike,
        options: ConsensusOptions | None = None,
        *,
        index: Csi | PathLike | None = None,
    ) -> None:
        file = _bcf(variants)
        options = options or ConsensusOptions()
        given = reference
        self._owns_reference = not isinstance(reference, SequenceReader)
        if not isinstance(reference, SequenceReader):
            reference = SequenceReader(reference)
        o = options
        try:
            native = _native.ConsensusBuilder(
                os.fspath(file),
                None if index is None else os.fspath(index),
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
        self._open(native, file, given, options, index=index)
        self.source = file
        self.index = detect(native.index_path)
        self.reference = reference
        self.options = options
        self._samples = tuple(native.samples())

    @property
    def samples(self) -> tuple[str, ...]:
        return self._samples

    def fetch(self, region: RegionLike, alleles: Alleles | None = None) -> ConsensusSequence:
        """The consensus of ``region`` with ``alleles`` (default
        :class:`AllAlt`), 5' to 3' of the region's strand."""
        r = _key(region)
        alleles = AllAlt() if alleles is None else alleles
        bases, chain = self._native.sequence(r.chrom, r.start, r.end, r.strand, _selector(alleles))
        return ConsensusSequence(r, alleles, bases, chain)

    def fetch_batch(self, region: RegionLike, alleles: Iterable[Alleles]) -> ConsensusBatch:
        """The consensus of ``region`` with each of ``alleles``, each 5' to
        3' of the region's strand."""
        r = _key(region)
        alleles = tuple(alleles)
        offsets, bases, chains = self._native.sequences(
            r.chrom, r.start, r.end, r.strand, [_selector(a) for a in alleles]
        )
        return ConsensusBatch(r, alleles, offsets, bases, tuple(chains))

    def clone(self) -> ConsensusReader:
        """Share the BCF header, index and reference blocks with separate work state."""
        reference = self.reference.clone()
        try:
            native = self._native.clone_context(reference._native)
        except BaseException:
            reference.close()
            raise
        cloned = object.__new__(type(self))
        cloned.__dict__ = self.__dict__.copy()
        cloned._native = native
        cloned._closed = False
        cloned.reference = reference
        cloned._owns_reference = True
        return cloned

    def close(self) -> None:
        """Closes the reader, and its reference if it opened it."""
        super().close()
        if self._owns_reference:
            self.reference.close()
