"""Sequences as the native layer scans them: ASCII bytes in Arrow
``LargeString`` layout."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from ..types import Sequence, StringVector

__all__ = ["SequenceLike", "Sequences", "as_strings"]

SequenceLike = str | bytes | bytearray | memoryview | np.ndarray | Sequence
"""One sequence: ASCII text, bytes, a ``uint8`` array or a
:class:`~omicsweave.types.Sequence`."""

Sequences = StringVector | SequenceLike | Iterable[SequenceLike]
"""Sequences to scan: a :class:`~omicsweave.types.StringVector` (scanned in
place), one sequence, or an iterable of them."""


def _bytes(seq: SequenceLike) -> np.ndarray:
    if isinstance(seq, Sequence):
        seq = seq.bases
    if isinstance(seq, np.ndarray):
        if seq.ndim != 1:
            raise ValueError("a sequence array is one-dimensional")
        return np.ascontiguousarray(seq, dtype=np.uint8)
    if isinstance(seq, str):
        return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    if isinstance(seq, bytes | bytearray | memoryview):
        return np.frombuffer(bytes(seq), dtype=np.uint8)
    raise TypeError(
        f"expected a sequence (str, bytes, uint8 array or Sequence), not {type(seq).__name__}"
    )


def as_strings(sequences: Sequences) -> StringVector:
    """``sequences`` as one :class:`~omicsweave.types.StringVector`."""
    if isinstance(sequences, StringVector):
        offsets = np.ascontiguousarray(sequences.offsets, dtype=np.int64)
        data = np.ascontiguousarray(sequences.data, dtype=np.uint8)
        return StringVector(offsets, data)
    if isinstance(sequences, SequenceLike):
        parts = [_bytes(sequences)]
    else:
        parts = [_bytes(s) for s in sequences]
    offsets = np.zeros(len(parts) + 1, dtype=np.int64)
    np.cumsum([len(p) for p in parts], out=offsets[1:])
    if len(parts) == 1:
        data = parts[0]
    else:
        data = np.concatenate(parts) if parts else np.zeros(0, dtype=np.uint8)
    return StringVector(offsets, data)
