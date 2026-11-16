"""Juicer ``.hic`` (versions 8 and 9): binned contact maps at several resolutions."""

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from .. import _native
from .._budget import threads as _threads
from ..types import BinUnit
from ._base import ChromIndex, File, PathLike, Writer
from ._sidecar import Sidecar, build_sidecar, written_sidecar

__all__ = ["Hic", "HicBlocks", "HicBuilder", "HicHeader", "HicWriter"]


class Hic(File):
    """A ``.hic`` file."""

    __slots__ = ()
    format = "hic"

    def build_sidecar(self) -> Sidecar:
        """Builds ``<file>.sidecar``, reading the whole file."""
        return build_sidecar(self)


@dataclass(frozen=True, slots=True)
class HicHeader:
    """What a new ``.hic`` file declares before any map is written."""

    genome: str
    chroms: Mapping[str, int]
    """Chromosome lengths by name, in file order."""
    bp_resolutions: tuple[int, ...]
    frag_resolutions: tuple[int, ...] = ()
    sites: Mapping[str, np.ndarray] | None = None
    """Sorted restriction site positions per chromosome; required with
    fragment resolutions."""
    version: int = 9
    attributes: Mapping[str, str] = field(default_factory=lambda: {"software": "omicsweave"})


def _native_header(header: HicHeader) -> tuple[object, ...]:
    sites = []
    if header.sites is not None:
        sites = [np.asarray(header.sites[name], dtype=np.uint32).tolist() for name in header.chroms]
    return (
        header.version,
        header.genome,
        list(header.attributes.items()),
        list(header.chroms.items()),
        list(header.bp_resolutions),
        list(header.frag_resolutions),
        sites,
    )


def _u32(a: object) -> np.ndarray:
    return np.ascontiguousarray(np.asarray(a, dtype=np.uint32).ravel())


@dataclass(frozen=True, slots=True)
class HicBlocks:
    """How a ``.hic`` writer cuts maps into blocks: ``intra`` (within a
    chromosome) and ``inter`` (between chromosomes) bins per block side, or
    ``"juicer"`` for the grid the Juicer tools choose (about 1000 bins). The
    writers' default is 128-bin blocks within chromosomes, where loaders
    read windows along the diagonal, and the Juicer grid between them. Any
    reader of the format takes a map's block size from the file."""

    intra: int | Literal["juicer"]
    inter: int | Literal["juicer"] = "juicer"


def _blocks(blocks: HicBlocks | None) -> tuple[int | None, int | None] | None:
    if blocks is None:
        return None
    return (_bins(blocks.intra), _bins(blocks.inter))


def _bins(bins: int | Literal["juicer"]) -> int | None:
    return None if bins == "juicer" else int(bins)


class HicWriter(Writer[Hic]):
    """Writes a ``.hic`` file from binned maps: blocks of ``level`` zlib
    (0-9; the Juicer tools use 6), laid out as ``blocks`` says; ``threads``
    above 1 encode and compress blocks on that many workers instead of the
    calling thread. ``sidecar`` also saves its sidecar on :meth:`close`,
    built as the file is written."""

    def __init__(
        self,
        path: PathLike,
        header: HicHeader,
        *,
        level: int = 6,
        blocks: HicBlocks | None = None,
        sidecar: bool = False,
        threads: int = 1,
    ) -> None:
        native = _native.HicWriter(
            os.fspath(path),
            *_native_header(header),
            level,
            blocks=_blocks(blocks),
            threads=_threads(threads),
        )
        super().__init__(Hic(path), native)
        self._chroms = ChromIndex(header.chroms)
        self._sidecar = sidecar

    def write_matrix(
        self,
        chrom1: str,
        chrom2: str,
        bin_size: int,
        x: np.ndarray,
        y: np.ndarray,
        value: np.ndarray,
        unit: BinUnit | str = BinUnit.BP,
    ) -> None:
        """Writes the map of ``chrom1`` x ``chrom2`` at ``bin_size``: cells at
        bin ``x`` of ``chrom1`` and bin ``y`` of ``chrom2``, in any order and
        orientation (cells of one stored position add up)."""
        self._native.write_matrix(
            self._chroms[chrom1],
            self._chroms[chrom2],
            str(BinUnit(unit)),
            bin_size,
            _u32(x),
            _u32(y),
            np.ascontiguousarray(np.asarray(value, dtype=np.float32).ravel()),
        )

    def _finish(self) -> tuple[File, ...]:
        return written_sidecar(self._native.close(self._sidecar))


class HicBuilder(Writer[Hic]):
    """Builds a ``.hic`` file from contacts, binned at every resolution of
    ``header`` when closed (as the Juicer tools' ``pre`` does). ``sidecar``
    also saves its sidecar on :meth:`close`, built as the file is written."""

    def __init__(
        self,
        path: PathLike,
        header: HicHeader,
        *,
        level: int = 6,
        blocks: HicBlocks | None = None,
        sidecar: bool = False,
        threads: int = 1,
    ) -> None:
        native = _native.HicBuilder(
            os.fspath(path),
            *_native_header(header),
            level,
            blocks=_blocks(blocks),
            threads=_threads(threads),
        )
        super().__init__(Hic(path), native)
        self._chroms = ChromIndex(header.chroms)
        self._sidecar = sidecar

    def _chrom(self, chrom: str | np.ndarray) -> np.ndarray:
        if isinstance(chrom, str):
            return np.asarray(self._chroms[chrom], dtype=np.uint32)
        return np.asarray(chrom, dtype=np.uint32)

    def add(
        self,
        chrom1: str | np.ndarray,
        pos1: np.ndarray | int,
        chrom2: str | np.ndarray,
        pos2: np.ndarray | int,
        count: np.ndarray | float = 1.0,
    ) -> None:
        """Adds ``count`` contacts between base ``pos1`` of ``chrom1`` and
        base ``pos2`` of ``chrom2``. Chromosomes are names or arrays of
        header indexes; all arguments broadcast."""
        c1, p1, c2, p2, n = np.broadcast_arrays(
            self._chrom(chrom1),
            np.asarray(pos1, dtype=np.uint32),
            self._chrom(chrom2),
            np.asarray(pos2, dtype=np.uint32),
            np.asarray(count, dtype=np.float32),
        )
        self._native.add(_u32(c1), _u32(p1), _u32(c2), _u32(p2), np.ascontiguousarray(n.ravel()))

    def _finish(self) -> tuple[File, ...]:
        return written_sidecar(self._native.close(self._sidecar))
