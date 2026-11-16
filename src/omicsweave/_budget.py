"""The budget that conversions and quantifications run within, and how they
divide it."""

from __future__ import annotations

import operator
import os
from dataclasses import dataclass
from pathlib import Path

from . import _native

PathLike = str | os.PathLike[str]

__all__ = ["Budget", "MemoryPlan"]


@dataclass(frozen=True, slots=True)
class Budget:
    """The threads, memory and spill directory a run may use."""

    threads: int
    """Worker threads for the parallel stages."""
    memory: int
    """Bytes for buffers."""
    tmp_dir: Path | None
    """Where spill files go; ``None``: beside the output."""


@dataclass(frozen=True, slots=True)
class MemoryPlan:
    """How a run divides its memory budget, in bytes."""

    fixed: int
    """Kept for the whole run: indexes, headers, tables."""
    pipeline: int
    """Batches in flight between reading, the worker threads and writing."""
    buffer: int
    """Sort buffers and spools, spilled to disk when full."""

    @property
    def total(self) -> int:
        return self.fixed + self.pipeline + self.buffer


def threads(value: int | None) -> int | None:
    """``value`` checked as a thread count; ``None`` leaves it to the native
    layer: every core the process may run on (CPU affinity and cgroup
    limits included)."""
    if value is None:
        return None
    if isinstance(value, bool):
        raise TypeError("threads must be an integer")
    n = operator.index(value)
    if n < 1:
        raise ValueError("threads must be at least 1")
    return n


def memory(value: int | str | None) -> int | None:
    """``value`` (bytes, or a size such as ``"8G"``) in bytes; ``None``
    leaves it to the native layer's default budget."""
    if value is None:
        return None
    if isinstance(value, str):
        return _native.parse_size(value)
    if isinstance(value, bool):
        raise TypeError("memory must be a byte count or a size such as '8G'")
    n = operator.index(value)
    if n < 0:
        raise ValueError("memory must be nonnegative")
    return n


def path(value: PathLike | None) -> str | None:
    return None if value is None else os.fspath(value)


def budget(native: tuple[int, int, str | None]) -> Budget:
    t, m, tmp = native
    return Budget(t, m, None if tmp is None else Path(tmp))


def memory_plan(native: tuple[int, int, int]) -> MemoryPlan:
    return MemoryPlan(*native)
