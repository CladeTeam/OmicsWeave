"""Motifs by number or ID, from an FFM collection."""

from __future__ import annotations

import os
from collections.abc import Iterator, Mapping
from types import MappingProxyType
from typing import Literal, overload

import numpy as np

from .. import _native
from ..formats import Ffm, PathLike
from ..types import Alphabet, Factor, Gaps, Motif
from ..types._motif import AnnotationValue
from ._base import Reader, resolve

__all__ = ["MotifReader"]

KeyType = Literal["string", "integer", "float"]


class MotifReader(Reader):
    """Reads an FFM motif collection (``.ffm``): motifs by number or ID as
    :class:`~omicsweave.types.Motif`, with their annotation typed by the
    collection's schema.

    The collection is read and checked whole on open and kept in memory
    (motif collections are small); lookups do no I/O. The reader is what
    :func:`omicsweave.motif.compile` takes to compile a collection without
    building Python objects for each motif.
    """

    def __init__(self, source: Ffm | PathLike) -> None:
        file = resolve(source, (Ffm,), "MotifReader")
        native = _native.MotifCollection(os.fspath(file))
        self._open(native, file)
        self.source = file
        self.alphabet = Alphabet(*native.alphabet())
        """The alphabet every motif is over."""
        self.schema: Mapping[str, KeyType] = MappingProxyType(dict(native.schema()))
        """Annotation keys and their types, in stored order."""
        self._ids: tuple[str, ...] = tuple(native.ids())
        self._index = {id: i for i, id in enumerate(self._ids)}

    @property
    def ids(self) -> tuple[str, ...]:
        """Motif IDs by number."""
        return self._ids

    def __len__(self) -> int:
        return len(self._ids)

    def __contains__(self, key: object) -> bool:
        return key in self._index

    def index(self, id: str) -> int:
        """The number of motif ``id``; :class:`KeyError` when absent."""
        return self._index[id]

    @overload
    def __getitem__(self, key: int | str) -> Motif: ...
    @overload
    def __getitem__(self, key: slice) -> list[Motif]: ...
    def __getitem__(self, key: int | str | slice) -> Motif | list[Motif]:
        if isinstance(key, slice):
            return [self._motif(i) for i in range(*key.indices(len(self)))]
        if isinstance(key, str):
            return self._motif(self._index[key])
        n = len(self)
        i = key + n if key < 0 else key
        if not 0 <= i < n:
            raise IndexError("motif index out of range")
        return self._motif(i)

    def get(self, id: str, default: Motif | None = None) -> Motif | None:
        i = self._index.get(id)
        return default if i is None else self._motif(i)

    def __iter__(self) -> Iterator[Motif]:
        for i in range(len(self)):
            yield self._motif(i)

    def _motif(self, i: int) -> Motif:
        width, support, sites, factors, gaps, values = self._native.motif(i)
        annotation: dict[str, AnnotationValue] = {}
        for (name, kind), raw in zip(self.schema.items(), values, strict=True):
            if raw == "":
                continue
            annotation[name] = (
                int(raw) if kind == "integer" else float(raw) if kind == "float" else raw
            )
        return Motif(
            self._ids[i],
            width,
            tuple(
                Factor(kind, tuple(positions.tolist()), table) for kind, positions, table in factors
            ),
            support=support,
            sites=sites,
            gaps=None if gaps is None else Gaps(*gaps),
            annotation=annotation,
        )

    def widths(self) -> np.ndarray:
        """Each motif's width (``uint32``), by number."""
        return self._native.widths()

    def __repr__(self) -> str:
        return f"MotifReader({str(self.source.path)!r}, {len(self)} motifs, {self.alphabet.name})"
