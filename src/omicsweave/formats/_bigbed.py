"""bigBed (``.bb``): BED features in a BBI container."""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np

from .. import _native
from .._budget import threads as _threads
from ..types import StringVector
from ._base import ChromIndex, File, PathLike, Writer
from ._bbi import BbiOptions, native_options
from ._sidecar import Sidecar, build_sidecar, written_sidecar

__all__ = ["BigBed", "BigBedSchema", "BigBedWriter"]


class BigBed(File):
    """A bigBed file."""

    __slots__ = ()
    format = "bigbed"

    def build_sidecar(self) -> Sidecar:
        """Builds ``<file>.sidecar``, reading the whole file."""
        return build_sidecar(self)


@dataclass(frozen=True, slots=True)
class BigBedSchema:
    """The record layout a bigBed header declares."""

    field_count: int = 3
    """Fields per record, the three coordinates included."""
    defined_field_count: int | None = None
    """Leading fields that follow the BED standard; ``None``: all of them."""
    auto_sql: str | None = None
    """The autoSql table describing every field."""


def _strings(annotation: StringVector | Iterable[str] | str | None, n: int) -> StringVector:
    if annotation is None:
        return StringVector.repeat("", n)
    if isinstance(annotation, str):
        return StringVector.repeat(annotation, n)
    if isinstance(annotation, StringVector):
        return annotation
    return StringVector.from_list(annotation)


class BigBedWriter(Writer[BigBed]):
    """Writes a bigBed file over ``chroms`` (lengths by name, in file order).
    Features come sorted by chromosome and start. ``sidecar`` also saves its
    sidecar on :meth:`close`, built as the file is written. Blocks are
    compressed on ``threads`` workers, as :class:`BigWigWriter` does."""

    def __init__(
        self,
        path: PathLike,
        chroms: Mapping[str, int],
        schema: BigBedSchema | None = None,
        *,
        options: BbiOptions | None = None,
        sidecar: bool = False,
        threads: int = 1,
    ) -> None:
        schema = schema or BigBedSchema()
        defined = schema.defined_field_count
        native = _native.BigBedWriter(
            os.fspath(path),
            list(chroms.items()),
            schema.field_count,
            schema.field_count if defined is None else defined,
            schema.auto_sql,
            *native_options(options or BbiOptions()),
            threads=_threads(threads),
        )
        super().__init__(BigBed(path), native)
        self._chroms = ChromIndex(chroms)
        self._sidecar = sidecar

    def add(
        self,
        chrom: str,
        start: np.ndarray | int,
        end: np.ndarray | int,
        annotation: StringVector | Iterable[str] | str | None = None,
    ) -> None:
        """Adds features over ``[start, end)`` of ``chrom``, with their fields
        after the third (tab-separated) as ``annotation``: one string for
        every feature, or one per feature."""
        s, e = np.broadcast_arrays(
            np.asarray(start, dtype=np.uint32), np.asarray(end, dtype=np.uint32)
        )
        s, e = np.ascontiguousarray(s.ravel()), np.ascontiguousarray(e.ravel())
        strings = _strings(annotation, len(s))
        self._native.add(
            self._chroms[chrom],
            s,
            e,
            np.ascontiguousarray(strings.offsets, dtype=np.int64),
            np.ascontiguousarray(strings.data, dtype=np.uint8),
        )

    def _finish(self) -> tuple[File, ...]:
        return written_sidecar(self._native.close(self._sidecar))
