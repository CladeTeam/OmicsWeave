"""OmicsWeave sidecar (``.sidecar``): what a loader must know of a file
before reading any of its data blocks."""

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from .. import _native
from ._base import File, PathLike, companion_path, ensure_format

__all__ = ["BbiLevelBlocks", "BbiSidecar", "HicMapBlocks", "HicSidecar", "Sidecar", "read_sidecar"]


class Sidecar(File):
    """A ``.sidecar`` of a bigWig, bigBed or ``.hic`` file: every data
    block's size once inflated, and what identifies the file it was built for.

    Only loaders use sidecars; reading by key does not.
    """

    __slots__ = ()
    format = "sidecar"
    suffix = ".sidecar"

    @classmethod
    def beside(cls, source: PathLike) -> Sidecar:
        """The sidecar named after ``source``: ``<source>.sidecar``."""
        return cls(companion_path(source, cls.suffix))


def build_sidecar(source: File) -> Sidecar:
    """Builds the sidecar of ``source`` beside it, reading the whole file."""
    ensure_format(source)
    return Sidecar(_native.build_sidecar(os.fspath(source)))


def written_sidecar(path: str | None) -> tuple[File, ...]:
    """The companion files of a writer that saves its sidecar at ``path``, if any."""
    return () if path is None else (Sidecar(path),)


@dataclass(frozen=True, slots=True, eq=False)
class BbiLevelBlocks:
    """The data blocks of one level of a bigWig or bigBed file, in file
    order: ``reduction`` bases per zoom summary (0 for the base level), and
    per block the (chromosome id, base) span the index records
    (``start_chrom``, ``start``, ``end_chrom``, ``end``), its ``offset``
    and ``size`` in the file, its ``inflated`` size and its ``items``."""

    reduction: int
    blocks: Mapping[str, np.ndarray]

    def __len__(self) -> int:
        return len(self.blocks["offset"])


@dataclass(frozen=True, slots=True, eq=False)
class BbiSidecar:
    """What a bigWig or bigBed sidecar holds: the header's uncompress buffer
    size (0 for uncompressed files) and the blocks of each level, the base
    level first."""

    uncompress_buf_size: int
    levels: tuple[BbiLevelBlocks, ...]


@dataclass(frozen=True, slots=True, eq=False)
class HicMapBlocks:
    """The blocks of one ``.hic`` map: chromosome ids (header order), unit
    and bin size, the block grid (``block_bin_count`` bins per block side,
    ``block_column_count`` blocks per row), and per block its ``number``,
    ``position`` and ``size`` in the file, ``inflated`` size and
    ``records``."""

    chrom1: int
    chrom2: int
    unit: str
    resolution: int
    block_bin_count: int
    block_column_count: int
    blocks: Mapping[str, np.ndarray]

    def __len__(self) -> int:
        return len(self.blocks["number"])


@dataclass(frozen=True, slots=True, eq=False)
class HicSidecar:
    """What a ``.hic`` sidecar holds: the file version and every map's
    blocks, sorted by (chrom1, chrom2, unit, resolution)."""

    version: int
    maps: tuple[HicMapBlocks, ...]


def read_sidecar(source: File | PathLike) -> BbiSidecar | HicSidecar:
    """The sidecar of ``source`` (a bigWig, bigBed or ``.hic`` file), read
    from beside it and checked against the file it was built for."""
    kind, head, parts = _native.read_sidecar(os.fspath(source))
    if kind == "bbi":
        levels = tuple(BbiLevelBlocks(r, MappingProxyType(c)) for r, c in parts)
        return BbiSidecar(head, levels)
    maps = tuple(HicMapBlocks(*m[:6], MappingProxyType(m[6])) for m in parts)
    return HicSidecar(head, maps)
