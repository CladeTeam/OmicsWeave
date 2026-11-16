"""Aligned reads: the records of a region, column by column."""

from __future__ import annotations

import enum
from dataclasses import dataclass

import numpy as np

from ._genome import Region
from ._string import StringVector

__all__ = [
    "AlignmentField",
    "Alignments",
    "Cigars",
    "ModificationFlag",
    "Modifications",
    "SamFlag",
]

CIGAR_OPS = "MIDNSHP=X"
"""CIGAR operation letters by code (``op & 0xf``)."""


class SamFlag(enum.IntFlag):
    """The bits of a SAM record's FLAG (:attr:`Alignments.flag`, and the
    read selection of :class:`omicsweave.quant.Filter`)."""

    PAIRED = 0x1
    PROPER_PAIR = 0x2
    UNMAPPED = 0x4
    MATE_UNMAPPED = 0x8
    REVERSE = 0x10
    MATE_REVERSE = 0x20
    READ1 = 0x40
    READ2 = 0x80
    SECONDARY = 0x100
    QC_FAIL = 0x200
    DUPLICATE = 0x400
    SUPPLEMENTARY = 0x800


class AlignmentField(enum.IntFlag):
    """The variable-length fields an alignment query can return besides the
    fixed ones (positions, flag, MAPQ, mate, template length)."""

    NAME = 1
    SEQ = 2
    QUAL = 4
    CIGAR = 8
    AUX = 16
    MODIFICATIONS = 32


class ModificationFlag(enum.IntFlag):
    """The bits of :attr:`Modifications.flags`."""

    MINUS = 1
    """The modification is on the strand opposite the sequenced one (an
    ``MM`` entry such as ``G-m``)."""
    EXPLICIT = 2
    """The ``MM`` entry is explicit (``?``): bases it does not list were not
    examined, rather than unmodified."""
    UNSCORED = 4
    """The record has no ``ML`` tag; :attr:`Modifications.probability` is
    0."""


@dataclass(frozen=True, slots=True, eq=False)
class Cigars:
    """CIGAR operations of every record in BAM encoding (``len << 4 | op``,
    ``uint32``): record ``i``'s are ``ops[offsets[i]:offsets[i + 1]]``."""

    offsets: np.ndarray
    ops: np.ndarray

    def __len__(self) -> int:
        return len(self.offsets) - 1

    def __getitem__(self, i: int) -> np.ndarray:
        n = len(self)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError("record index out of range")
        return self.ops[self.offsets[i] : self.offsets[i + 1]]

    def string(self, i: int) -> str:
        """Record ``i``'s CIGAR as SAM writes it (``*`` when empty)."""
        ops = self[i]
        if len(ops) == 0:
            return "*"
        return "".join(f"{v >> 4}{CIGAR_OPS[v & 0xF]}" for v in ops.tolist())


def modification_code(code: str | int) -> int:
    """The :attr:`Modifications.code` of a modification: a one-letter
    ``MM`` code (``"m"``, ``"h"``, ``"a"``) or a ChEBI identifier (``17596``
    or ``"17596"``)."""
    if isinstance(code, bool):
        raise TypeError("a modification code is a letter or a ChEBI number")
    if isinstance(code, int):
        if code <= 0:
            raise ValueError(f"ChEBI identifiers are positive, not {code}")
        return -code
    if isinstance(code, str):
        if len(code) == 1 and code.isascii() and code.isalpha():
            return ord(code)
        if code.isascii() and code.isdigit():
            return modification_code(int(code))
        raise ValueError(f"invalid modification code {code!r}: use a letter or a ChEBI number")
    raise TypeError("a modification code is a letter or a ChEBI number")


@dataclass(frozen=True, slots=True, eq=False)
class Modifications:
    """Base modifications from the ``MM`` / ``ML`` tags (SAMtags §1.7) of
    every record, decoded as htslib's ``bam_next_basemod`` reports them:
    record ``i``'s are rows ``offsets[i]:offsets[i + 1]``, in SEQ order
    (by position, then ``MM`` entry, then code).

    Columns: ``ref_pos`` (``int64``, the reference position of the base,
    -1 when it is not aligned: inserted, soft-clipped or unmapped),
    ``qpos`` (``int32``, the index of the base in the stored SEQ, as
    :attr:`Alignments.seqs` holds it), ``code`` (``int32``, a letter code as
    its ASCII value, ``ord("m")`` = 109, or a ChEBI identifier negated;
    see :meth:`code_name`), ``canonical`` (``uint8``, the unmodified base
    of the ``MM`` entry as ASCII, on the strand sequenced), ``probability``
    (``uint8``, the raw ``ML`` value; probability ≈ ``(v + 0.5) / 256``)
    and ``flags`` (``uint8``, :class:`ModificationFlag` bits).
    """

    offsets: np.ndarray
    ref_pos: np.ndarray
    qpos: np.ndarray
    code: np.ndarray
    canonical: np.ndarray
    probability: np.ndarray
    flags: np.ndarray

    def __len__(self) -> int:
        return len(self.offsets) - 1

    def rows(self, i: int) -> slice:
        """The rows of record ``i``."""
        n = len(self)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError("record index out of range")
        return slice(int(self.offsets[i]), int(self.offsets[i + 1]))

    @property
    def record(self) -> np.ndarray:
        """The record index of every row (``int64``)."""
        return np.repeat(np.arange(len(self), dtype=np.int64), np.diff(self.offsets))

    @staticmethod
    def code_name(code: int) -> str:
        """``"m"`` for ``ord("m")``, ``"17596"`` for ``-17596``."""
        return chr(code) if code > 0 else str(-code)


@dataclass(frozen=True, slots=True, eq=False)
class Alignments:
    """The records overlapping ``region`` in file order: exactly those
    htslib's region iterator returns, unmapped records placed there
    included.

    Fixed columns: ``start`` (0-based POS, ``int64``), ``end`` (end of the
    reference span, as htslib's ``bam_endpos``, ``int64``), ``flag``
    (``uint16``), ``mapq`` (``uint8``), ``mate_tid`` (``int32``, -1 for
    none; an index into the reader's references), ``mate_start``
    (``int32``, -1 for none) and ``tlen`` (``int32``). The variable-length
    fields are present only when the query asked for them: ``names``,
    ``seqs`` (bases as ASCII), ``quals`` (raw Phred bytes, ``0xff``
    throughout when absent; with the same lengths as ``seqs``), ``cigars``,
    ``aux`` (optional fields in their binary BAM encoding) and
    ``modifications`` (:class:`Modifications`).
    """

    region: Region
    start: np.ndarray
    end: np.ndarray
    flag: np.ndarray
    mapq: np.ndarray
    mate_tid: np.ndarray
    mate_start: np.ndarray
    tlen: np.ndarray
    names: StringVector | None = None
    seqs: StringVector | None = None
    quals: StringVector | None = None
    cigars: Cigars | None = None
    aux: StringVector | None = None
    modifications: Modifications | None = None

    def __len__(self) -> int:
        return len(self.start)
