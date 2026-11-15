"""Aligned reads: the records of a region, column by column."""

from __future__ import annotations

import enum
from dataclasses import dataclass

import numpy as np

from ._genome import Region
from ._string import StringVector

__all__ = ["AlignmentField", "Alignments", "Cigars"]

CIGAR_OPS = "MIDNSHP=X"
"""CIGAR operation letters by code (``op & 0xf``)."""


class AlignmentField(enum.IntFlag):
    """The variable-length fields an alignment query can return besides the
    fixed ones (positions, flag, MAPQ, mate, template length)."""

    NAME = 1
    SEQ = 2
    QUAL = 4
    CIGAR = 8
    AUX = 16


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
    throughout when absent; with the same lengths as ``seqs``), ``cigars``
    and ``aux`` (optional fields in their binary BAM encoding).
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

    def __len__(self) -> int:
        return len(self.start)
