"""A column of strings stored as Arrow stores them."""

from __future__ import annotations

import collections.abc as abc
from dataclasses import dataclass
from typing import Any, overload

import numpy as np

__all__ = ["StringVector"]


@dataclass(frozen=True, slots=True, eq=False)
class StringVector(abc.Sequence[str]):
    """Strings in Arrow ``LargeString`` layout.

    String ``i`` is ``data[offsets[i]:offsets[i + 1]]`` (UTF-8); ``offsets``
    (``int64``) has one more element than there are strings and starts at 0.
    """

    offsets: np.ndarray
    data: np.ndarray

    def __len__(self) -> int:
        return len(self.offsets) - 1

    @overload
    def __getitem__(self, i: int) -> str: ...
    @overload
    def __getitem__(self, i: slice) -> list[str]: ...
    def __getitem__(self, i: int | slice) -> str | list[str]:
        if isinstance(i, slice):
            return [self[j] for j in range(*i.indices(len(self)))]
        return self.bytes(i).decode()

    def bytes(self, i: int) -> bytes:
        """The raw bytes of string ``i``."""
        n = len(self)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError("string index out of range")
        return self.data[self.offsets[i] : self.offsets[i + 1]].tobytes()

    def to_list(self) -> list[str]:
        return [self[i] for i in range(len(self))]

    def to_arrow(self) -> Any:
        """The strings as a ``pyarrow.LargeStringArray`` sharing this memory."""
        import pyarrow as pa

        return pa.LargeStringArray.from_buffers(
            len(self), pa.py_buffer(self.offsets), pa.py_buffer(self.data)
        )

    @classmethod
    def from_list(cls, strings: abc.Iterable[str | bytes]) -> StringVector:
        encoded = [s.encode() if isinstance(s, str) else bytes(s) for s in strings]
        offsets = np.zeros(len(encoded) + 1, dtype=np.int64)
        np.cumsum([len(b) for b in encoded], out=offsets[1:])
        data = np.frombuffer(b"".join(encoded), dtype=np.uint8).copy()
        return cls(offsets, data)

    @classmethod
    def repeat(cls, string: str, n: int) -> StringVector:
        """``string``, ``n`` times; it is not encoded when ``n`` is 0."""
        encoded = string.encode() if n else b""
        if encoded:
            offsets = np.arange(n + 1, dtype=np.int64)
            offsets *= len(encoded)
        else:
            offsets = np.zeros(n + 1, dtype=np.int64)
        data = np.frombuffer(encoded * n, dtype=np.uint8).copy()
        return cls(offsets, data)

    def __repr__(self) -> str:
        shown = ", ".join(repr(s) for s in self[:3])
        more = ", ..." if len(self) > 3 else ""
        return f"StringVector([{shown}{more}], n={len(self)})"
