"""UCSC 2bit (``.2bit``): packed bases with N and soft-mask blocks."""

from __future__ import annotations

import os
from collections.abc import Iterable

import numpy as np

from .. import _native
from ._base import File, PathLike, Writer

__all__ = ["TwoBit", "TwoBitWriter"]


class TwoBit(File):
    """A ``.2bit`` file; it indexes itself, so it is read by region as it is."""

    __slots__ = ()
    format = "2bit"


class TwoBitWriter(Writer[TwoBit]):
    """Writes a 2bit file of the sequences ``names``, added in that order
    (the file's directory precedes its records). Runs of N and of lowercase
    (soft-masked) bases are kept as blocks; bases other than A, C, G and T
    become N, as the UCSC tools store them. ``wide_offsets`` stores 64-bit
    offsets (2bit version 1), needed when any record starts beyond
    ``2**32 - 1`` bytes. The final record may extend past 4 GiB with
    32-bit offsets."""

    def __init__(self, path: PathLike, names: Iterable[str], *, wide_offsets: bool = False) -> None:
        native = _native.TwoBitWriter(os.fspath(path), list(names), wide_offsets)
        super().__init__(TwoBit(path), native)

    @property
    def next_name(self) -> str | None:
        """The sequence to add next; ``None`` once every one is."""
        return self._native.next_name()

    def add(self, name: str, seq: str | bytes | np.ndarray) -> None:
        """Adds sequence ``name``, the next of ``names``."""
        if isinstance(seq, str):
            seq = seq.encode("ascii")
        elif isinstance(seq, np.ndarray):
            seq = np.ascontiguousarray(seq, dtype=np.uint8).tobytes()
        self._native.add(name, bytes(seq))

    def _finish(self) -> tuple[File, ...]:
        self._native.close()
        return ()
