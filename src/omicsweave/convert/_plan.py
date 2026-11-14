"""Checking a conversion against its budget before anything is read."""

from __future__ import annotations

from dataclasses import dataclass

from .. import _native
from ..formats import File, PathLike
from ._files import file_of, native_files


@dataclass(frozen=True, slots=True)
class Budget:
    """The threads and memory a conversion may use."""

    threads: int
    """Worker threads for the parallel stages."""
    memory: int
    """Bytes for buffers."""


@dataclass(frozen=True, slots=True)
class MemoryPlan:
    """How a conversion divides its memory budget, in bytes."""

    fixed: int
    """Kept for the whole run: indexes, headers, tables."""
    pipeline: int
    """Batches in flight between reading, the worker threads and writing."""
    buffer: int
    """Sort buffer, spilled to disk when full (0 for streaming conversions)."""


@dataclass(frozen=True, slots=True)
class Plan:
    """A conversion checked against its budget."""

    source: File
    target: File
    budget: Budget
    memory: MemoryPlan


def plan(
    source: File | PathLike,
    target: File | PathLike,
    *,
    threads: int | None = None,
    memory: int | str = "4G",
) -> Plan:
    """Checks the conversion of ``source`` into ``target`` against its budget
    (see :func:`run`) without reading any data."""
    files = native_files(source, target)
    src, dst, threads, budget, fixed, pipeline, buffer = _native.plan_convert(
        **files, threads=threads, memory=str(memory)
    )
    return Plan(
        file_of(src, files["src"]),
        file_of(dst, files["dst"]),
        Budget(threads, budget),
        MemoryPlan(fixed, pipeline, buffer),
    )
