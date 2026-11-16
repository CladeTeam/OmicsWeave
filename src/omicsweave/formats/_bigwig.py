"""bigWig (``.bw``): signal over a genome in a BBI container."""

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .. import _native
from .._budget import threads as _threads
from ._base import ChromIndex, File, PathLike, Writer
from ._bbi import BbiOptions, native_options
from ._sidecar import Sidecar, build_sidecar, written_sidecar

__all__ = ["BigWig", "BigWigSummary", "BigWigWriter"]


class BigWig(File):
    """A bigWig file."""

    __slots__ = ()
    format = "bigwig"

    def build_sidecar(self) -> Sidecar:
        """Builds ``<file>.sidecar``, reading the whole file."""
        return build_sidecar(self)


@dataclass(frozen=True, slots=True)
class BigWigSummary:
    """Statistics over every covered base, as a bigWig header stores them."""

    bases_covered: int
    min: float
    max: float
    sum: float
    sum_squares: float


class BigWigWriter(Writer[BigWig]):
    """Writes a bigWig file over ``chroms`` (lengths by name; their order is
    the file's chromosome order). Items come sorted by chromosome and start
    and do not overlap. ``sidecar`` also saves its sidecar on :meth:`close`,
    built as the file is written. ``threads`` above 1 compress blocks on that
    many workers instead of the calling thread; the file does not depend on
    the count."""

    def __init__(
        self,
        path: PathLike,
        chroms: Mapping[str, int],
        *,
        options: BbiOptions | None = None,
        sidecar: bool = False,
        threads: int = 1,
    ) -> None:
        native = _native.BigWigWriter(
            os.fspath(path),
            list(chroms.items()),
            *native_options(options or BbiOptions()),
            threads=_threads(threads),
        )
        super().__init__(BigWig(path), native)
        self._chroms = ChromIndex(chroms)
        self._sidecar = sidecar

    def add(
        self,
        chrom: str,
        start: np.ndarray | int,
        end: np.ndarray | int,
        value: np.ndarray | float,
    ) -> None:
        """Adds ``value`` over ``[start, end)`` of ``chrom``: scalars or columns."""
        s, e, v = np.broadcast_arrays(
            np.asarray(start, dtype=np.uint32),
            np.asarray(end, dtype=np.uint32),
            np.asarray(value, dtype=np.float32),
        )
        self._native.add(
            self._chroms[chrom],
            np.ascontiguousarray(s.ravel()),
            np.ascontiguousarray(e.ravel()),
            np.ascontiguousarray(v.ravel()),
        )

    def _finish(self) -> tuple[File, ...]:
        return written_sidecar(self._native.close(self._sidecar))
