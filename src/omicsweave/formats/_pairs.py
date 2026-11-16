"""Contact pairs as text: 4DN pairs, Juicer short and HiC-Pro valid pairs,
each plain or compressed with gzip / BGZF.

All list one contact per row and convert to ``.hic``; a pairix index
beside a ``.pairs.gz`` is not needed.
"""

from __future__ import annotations

from ._base import File

__all__ = ["GzipHicProPairs", "GzipJuicerShort", "GzipPairs", "HicProPairs", "JuicerShort", "Pairs"]


class Pairs(File):
    """Plain-text 4DN pairs (``.pairs``): ``#chromsize:`` and ``#columns:``
    header lines, then one row per contact."""

    __slots__ = ()
    format = "pairs"


class GzipPairs(File):
    """4DN pairs compressed with gzip or bgzip (``.pairs.gz``)."""

    __slots__ = ()
    format = "gzip-pairs"


class JuicerShort(File):
    """Juicer short text: ``str1 chr1 pos1 frag1 str2 chr2 pos2 frag2``
    rows, with an optional ninth column, the score."""

    __slots__ = ()
    format = "juicer-short"


class GzipJuicerShort(File):
    """Juicer short text compressed with gzip or bgzip."""

    __slots__ = ()
    format = "gzip-juicer-short"


class HicProPairs(File):
    """HiC-Pro ``.allValidPairs`` / ``.validPairs`` contact rows.

    Conversion needs chromosome sizes and counts each row once, without
    filtering on mapping quality, restriction fragment or allele. Rows have
    12 fields, optionally a 13th allele field; positions are non-negative
    u64 integers and are binned as written, clamped to chromosome lengths.
    """

    __slots__ = ()
    format = "hicpro-pairs"


class GzipHicProPairs(File):
    """HiC-Pro valid pairs compressed with gzip or bgzip."""

    __slots__ = ()
    format = "gzip-hicpro-pairs"
