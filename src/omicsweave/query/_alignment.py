"""Alignments by region from BAM and CRAM."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from types import MappingProxyType
from typing import Literal

import numpy as np

from .. import _native
from ..cache import BlockCache, BlockKey
from ..cache import native as native_cache
from ..formats import Bai, Bam, BgzipFasta, Crai, Cram, Csi, Fasta, PathLike, TwoBit, detect
from ..types import AlignmentField, Alignments, Cigars, RegionLike
from ..types._genome import as_region
from ..types._string import StringVector
from ._base import CachedReader, resolve
from ._sequence import SequenceReader

__all__ = ["AlignmentReader"]

Fields = AlignmentField | Literal["all"] | Iterable[str]

_FIELD_NAMES = {f.name.lower(): f for f in AlignmentField}  # type: ignore[union-attr]

UNMAPPED, SECONDARY, QC_FAIL, DUPLICATE = 0x4, 0x100, 0x200, 0x400

NO_FIELDS = AlignmentField(0)
ALL_FIELDS = AlignmentField(sum(AlignmentField))


def _fields(fields: Fields) -> AlignmentField:
    if isinstance(fields, AlignmentField):
        return fields
    if fields == "all":
        return ALL_FIELDS
    if isinstance(fields, str):
        fields = (fields,)
    out = AlignmentField(0)
    for name in fields:
        try:
            out |= _FIELD_NAMES[name]
        except KeyError:
            known = ", ".join(_FIELD_NAMES)
            raise ValueError(f"unknown alignment field {name!r}: use {known} or 'all'") from None
    return out


Reference = SequenceReader | Fasta | BgzipFasta | TwoBit | PathLike


def _reference(reference: Reference | None) -> object:
    if reference is None or isinstance(reference, SequenceReader):
        return None if reference is None else reference._native
    return os.fspath(reference)


def _index(index: Bai | Csi | Crai | PathLike | None) -> str | None:
    return None if index is None else os.fspath(index)


def _strings(pair: tuple[np.ndarray, np.ndarray] | None) -> StringVector | None:
    return None if pair is None else StringVector(*pair)


class AlignmentReader(CachedReader):
    """Reads the records of a BAM file by region, through its BAI or CSI
    index, or of a CRAM file through its ``.crai``.

    ``index`` defaults to the one htslib finds beside the file
    (``<bam>.csi``, ``<name>.csi``, ``<bam>.bai``, ``<name>.bai``;
    ``<cram>.crai``, ``<name>.crai``). A missing index raises
    :class:`~omicsweave.errors.IoError`: indexes are built explicitly
    (:meth:`Bam.build_index <omicsweave.formats.Bam.build_index>`,
    :meth:`Cram.build_index <omicsweave.formats.Cram.build_index>`), never
    on open.

    ``reference`` is the sequence a CRAM was compressed against: a
    :class:`SequenceReader` or a FASTA, bgzip FASTA or 2bit file. Its
    bases rebuild the reads' sequences (and the MD / NM tags htslib
    generates), each slice checked against its MD5; fixed fields, CIGARs,
    names and qualities need no reference, and CRAMs with embedded or no
    reference need none at all. BAM takes no reference.

    ``cache`` optionally supplies caller-owned decoded BGZF blocks (BAM
    only). The file and its index must remain unchanged throughout the
    reader session.
    """

    def __init__(
        self,
        source: Bam | Cram | PathLike,
        index: Bai | Csi | Crai | PathLike | None = None,
        *,
        reference: Reference | None = None,
        cache: BlockCache | None = None,
    ) -> None:
        file = resolve(source, (Bam, Cram), "AlignmentReader")
        native = _native.AlignmentReader(
            os.fspath(file), _index(index), native_cache(cache), _reference(reference)
        )
        if reference is None:
            self._open(native, file, index)
        else:
            self._open(native, file, index, reference=reference)
        self.source = file
        chroms = native.chroms()
        self._chroms = MappingProxyType(dict(chroms))
        self._references = tuple(name for name, _ in chroms)
        self.index = detect(native.index_path)

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach externally owned BGZF blocks, or detach them with None
        (a CRAM takes none)."""
        self._native.set_cache(native_cache(cache))

    def describe_cache(self, region: RegionLike) -> tuple[BlockKey, ...]:
        """The physical BGZF dependencies of the indexed chunks."""
        r = as_region(region)
        return tuple(self._native.describe_cache(r.chrom, r.start, r.end))

    def prefetch_cache(self, region: RegionLike) -> None:
        """Load reusable BGZF blocks without constructing alignment output."""
        r = as_region(region)
        self._native.prefetch_cache(r.chrom, r.start, r.end)

    def stream_stats(self) -> dict[str, int]:
        """Physical reads, inflations and buffer movement in this context."""
        return self._native.stream_stats()

    def reset_stream_stats(self) -> None:
        """Reset physical-work counters, retaining blocks and query state."""
        self._native.reset_stream_stats()

    @property
    def chroms(self) -> Mapping[str, int]:
        """Reference lengths by name, in reference-id order."""
        return self._chroms

    @property
    def references(self) -> tuple[str, ...]:
        """Reference names by id, the order ``mate_tid`` indexes."""
        return self._references

    @property
    def header_text(self) -> bytes:
        """The SAM header text as stored."""
        return self._native.header_text()

    def fetch(self, region: RegionLike, fields: Fields = NO_FIELDS) -> Alignments:
        """The records overlapping ``region``, with the variable-length
        ``fields`` asked for (:class:`~omicsweave.types.AlignmentField`
        flags, names such as ``("name", "seq")``, or ``"all"``)."""
        r = as_region(region)
        (start, end, flag, mapq, mate_tid, mate_start, tlen, names, seqs, quals, cigars, aux) = (
            self._native.query(r.chrom, r.start, r.end, int(_fields(fields)))
        )
        return Alignments(
            r,
            start,
            end,
            flag,
            mapq,
            mate_tid,
            mate_start,
            tlen,
            _strings(names),
            _strings(seqs),
            _strings(quals),
            None if cigars is None else Cigars(*cigars),
            _strings(aux),
        )
