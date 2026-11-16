"""Alignments by region from BAM and CRAM."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from types import MappingProxyType
from typing import Literal

import numpy as np

from .. import _native
from ..cache import BlockCache
from ..cache import native as native_cache
from ..formats import Bai, Bam, BgzipFasta, Crai, Cram, Csi, Fasta, PathLike, TwoBit, detect
from ..types import AlignmentField, Alignments, Cigars, Modifications, RegionLike
from ..types._alignment import modification_code
from ..types._genome import as_region
from ..types._string import StringVector
from ._base import CachedReader, resolve
from ._sequence import SequenceReader, reference_reader

__all__ = ["AlignmentReader"]

Fields = AlignmentField | Literal["all"] | Iterable[str]

_FIELD_NAMES = {f.name.lower(): f for f in AlignmentField}  # type: ignore[union-attr]

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
    :class:`SequenceReader` (shared, left open) or a FASTA, bgzip FASTA or
    2bit file (opened with the same ``cache``). Its bases rebuild the
    reads' sequences (and the MD / NM tags htslib generates), each slice
    checked against its MD5; fixed fields, CIGARs, names and qualities
    need no reference, and CRAMs with embedded or no reference need none
    at all. BAM takes no reference.

    CRAM decoding follows htslib's defaults: ``verify_checksums`` checks
    block CRC32s and slice reference MD5s, ``decode_md`` generates the MD
    and NM tags records do not store when querying ``"aux"``, and reads
    stored without a name are named ``<name_prefix>:<n>`` (default: the file's name).

    ``cache`` optionally supplies caller-owned decoded BGZF blocks (BAM)
    or prepared containers and decoded slices (CRAM). The file, index and
    reference must remain unchanged throughout the reader session.
    """

    def __init__(
        self,
        source: Bam | Cram | PathLike,
        index: Bai | Csi | Crai | PathLike | None = None,
        *,
        reference: Reference | None = None,
        verify_checksums: bool = True,
        decode_md: bool = True,
        name_prefix: str | None = None,
        cache: BlockCache | None = None,
    ) -> None:
        file = resolve(source, (Bam, Cram), "AlignmentReader")
        opened = None
        if reference is not None and not isinstance(reference, SequenceReader):
            opened = reference = reference_reader(reference, cache)
        try:
            native = _native.AlignmentReader(
                os.fspath(file),
                _index(index),
                native_cache(cache),
                None if reference is None else reference._native,
                verify_checksums=verify_checksums,
                decode_md=decode_md,
                name_prefix=name_prefix,
            )
        finally:
            # The native reader holds its own handle on the reference.
            if opened is not None:
                opened.close()
        options = {
            k: v
            for k, v in (
                ("reference", opened.source if opened is not None else reference),
                ("verify_checksums", verify_checksums),
                ("decode_md", decode_md),
                ("name_prefix", name_prefix),
            )
            if v is not None and v is not True
        }
        self._open(native, file, index, **options)
        self.source = file
        chroms = native.chroms()
        self._chroms = MappingProxyType(dict(chroms))
        self._references = tuple(name for name, _ in chroms)
        self.index = detect(native.index_path)

    def set_cache(self, cache: BlockCache | None) -> None:
        """Attach externally owned query blocks, or detach them with None."""
        self._native.set_cache(native_cache(cache))

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

    def fetch(
        self,
        region: RegionLike,
        fields: Fields = NO_FIELDS,
        *,
        modification_codes: Iterable[str | int] | None = None,
    ) -> Alignments:
        """The records overlapping ``region``, with the variable-length
        ``fields`` asked for (:class:`~omicsweave.types.AlignmentField`
        flags, names such as ``("name", "seq")``, or ``"all"``).

        ``"modifications"`` decodes the ``MM`` / ``ML`` tags into
        :class:`~omicsweave.types.Modifications`; ``modification_codes``
        keeps only some codes (letters such as ``"m"`` / ``"h"``, ChEBI
        numbers such as ``17596``). A CRAM needs its reference for them, as
        for ``"seq"``: positions count the read's bases."""
        r = as_region(region)
        wanted = _fields(fields)
        codes: list[int] = []
        if modification_codes is not None:
            if not wanted & AlignmentField.MODIFICATIONS:
                raise ValueError("modification_codes filters the 'modifications' field: ask for it")
            if isinstance(modification_codes, str | int):
                modification_codes = (modification_codes,)
            codes = [modification_code(c) for c in modification_codes]
        fixed, names, seqs, quals, cigars, aux, mods = self._native.query(
            r.chrom, r.start, r.end, int(wanted), codes
        )
        return Alignments(
            r,
            *fixed,
            _strings(names),
            _strings(seqs),
            _strings(quals),
            None if cigars is None else Cigars(*cigars),
            _strings(aux),
            None if mods is None else Modifications(*mods),
        )
