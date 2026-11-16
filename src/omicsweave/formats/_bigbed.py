"""bigBed (``.bb``): BED features in a BBI container."""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .. import _native
from .._budget import threads as _threads
from ..types import Intervals, StringVector
from ._base import ChromIndex, File, PathLike, Writer
from ._bbi import BbiOptions, native_options
from ._sidecar import Sidecar, build_sidecar, written_sidecar

__all__ = ["AutoSqlField", "BigBed", "BigBedSchema", "BigBedWriter"]


class BigBed(File):
    """A bigBed file."""

    __slots__ = ()
    format = "bigbed"

    def build_sidecar(self) -> Sidecar:
        """Builds ``<file>.sidecar``, reading the whole file."""
        return build_sidecar(self)


@dataclass(frozen=True, slots=True)
class AutoSqlField:
    """One field of an autoSql table."""

    name: str
    type: str
    """The autoSql type as declared: ``uint``, ``float``, ``string``,
    ``char[1]``, ``int[blockCount]``, ``enum(...)``, ..."""
    comment: str = ""

    @property
    def dtype(self) -> np.dtype | None:
        """The NumPy type of a scalar number field; ``None`` for text,
        arrays, enums and sets, and for the color fields (``reserved``,
        ``itemRgb``), which hold ``r,g,b`` text."""
        if self.name in _COLOR_FIELDS:
            return None
        kind = _NUMBERS.get(self.type)
        return None if kind is None else np.dtype(kind)


_NUMBERS = {
    "byte": np.int8,
    "ubyte": np.uint8,
    "short": np.int16,
    "ushort": np.uint16,
    "int": np.int32,
    "uint": np.uint32,
    "bigint": np.int64,
    "float": np.float32,
    "double": np.float64,
}
_COLOR_FIELDS = frozenset({"reserved", "itemRgb"})
_FIELD = re.compile(
    r'^\s*(?P<type>[A-Za-z]+(?:\[[^\]]*\]|\([^)]*\))?)\s+(?P<name>\w+)\s*;\s*(?:"(?P<comment>[^"]*)")?',
    re.MULTILINE,
)
# kent's bed.as: the standard BED fields after the coordinates.
_BED = (
    ("name", "string"),
    ("score", "uint"),
    ("strand", "char[1]"),
    ("thickStart", "uint"),
    ("thickEnd", "uint"),
    ("reserved", "uint"),
    ("blockCount", "int"),
    ("blockSizes", "int[blockCount]"),
    ("chromStarts", "int[blockCount]"),
)


@lru_cache(maxsize=64)
def _parse_fields(auto_sql: str | None, field_count: int) -> tuple[AutoSqlField, ...]:
    if auto_sql is not None:
        body = auto_sql[auto_sql.find("(") + 1 : auto_sql.rfind(")")]
        fields = tuple(
            AutoSqlField(m["name"], m["type"], m["comment"] or "") for m in _FIELD.finditer(body)
        )
        if len(fields) == field_count:
            return fields
    standard = (("chrom", "string"), ("chromStart", "uint"), ("chromEnd", "uint"), *_BED)
    return tuple(
        AutoSqlField(*standard[i]) if i < len(standard) else AutoSqlField(f"field{i + 1}", "string")
        for i in range(field_count)
    )


@dataclass(frozen=True, slots=True)
class BigBedSchema:
    """The record layout a bigBed header declares."""

    field_count: int = 3
    """Fields per record, the three coordinates included."""
    defined_field_count: int | None = None
    """Leading fields that follow the BED standard; ``None``: all of them."""
    auto_sql: str | None = None
    """The autoSql table describing every field."""

    @property
    def fields(self) -> tuple[AutoSqlField, ...]:
        """Every field, the three coordinates first, as the autoSql names
        them; without autoSql (or one that does not match the field count),
        the standard BED names and ``fieldN`` for the rest."""
        return _parse_fields(self.auto_sql, self.field_count)

    def columns(
        self, features: Intervals | StringVector, names: Sequence[str] | None = None
    ) -> dict[str, np.ndarray | StringVector]:
        """The fields after the third of ``features`` (a query result, or
        its ``annotation``), one column each, by autoSql name: scalar
        numbers as NumPy arrays of their type (:attr:`AutoSqlField.dtype`),
        everything else as :class:`~omicsweave.types.StringVector`.
        ``names`` picks fields (default: all after the third)."""
        annotation = features.annotation if isinstance(features, Intervals) else features
        fields = self.fields[3:]
        index = {f.name: i for i, f in enumerate(fields)}
        wanted = [f.name for f in fields] if names is None else list(names)
        for name in wanted:
            if name not in index:
                known = ", ".join(index)
                raise KeyError(f"no field {name!r} after the coordinates; the fields are {known}")
        split = _split_fields(annotation, len(fields))
        out: dict[str, np.ndarray | StringVector] = {}
        for name in wanted:
            field = fields[index[name]]
            values = split[index[name]]
            dtype = field.dtype
            if dtype is None:
                out[name] = StringVector.from_list(values)
            else:
                try:
                    out[name] = np.array(values, dtype=np.bytes_).astype(dtype)
                except ValueError as e:
                    raise ValueError(f"field {name!r} ({field.type}): {e}") from None
        return out


def _split_fields(annotation: StringVector, count: int) -> list[list[bytes]]:
    """The ``count`` tab-separated fields of every string, field by field."""
    n = len(annotation)
    if count == 0 or n == 0:
        return [[] for _ in range(count)]
    data = np.asarray(annotation.data, dtype=np.uint8)
    offsets = np.asarray(annotation.offsets, dtype=np.int64)
    # Each row must have its own separators: a surplus in one row cannot
    # make up for a missing field in another.
    tabs = np.flatnonzero(data == ord("\t"))
    counts = np.diff(np.searchsorted(tabs, offsets)) + 1
    bad = np.flatnonzero(counts != count)
    if len(bad):
        i = int(bad[0])
        raise ValueError(f"feature {i} has {counts[i]} fields after the third, not {count}")
    # One split of the whole buffer, a tab between the strings.
    joined = np.insert(data[offsets[0] : offsets[-1]], offsets[1:-1] - offsets[0], ord("\t"))
    words = joined.tobytes().split(b"\t")
    if len(words) != n * count:
        rows = [annotation.bytes(i).split(b"\t") for i in range(n)]
        bad = next(i for i, r in enumerate(rows) if len(r) != count)
        raise ValueError(f"feature {bad} has {len(rows[bad])} fields after the third, not {count}")
    return [words[j::count] for j in range(count)]


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
