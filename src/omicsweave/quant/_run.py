"""Planning and running a quantification task."""

from __future__ import annotations

import json
import os
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from .. import _native
from ..annotation import Annotation
from ..formats import PathLike
from ._products import Product

__all__ = ["MemoryPlan", "Plan", "Report", "Task", "plan", "run"]


@dataclass(frozen=True, slots=True)
class Task:
    """Coordinate-sorted BAMs, the products of one pass over them, and the
    directory they go to.

    ``annotation`` is needed by gene products (and the ATAC gene activity
    and TSS metrics); ``aliases`` match chromosome names that differ beyond
    the built-in ``chr1`` ↔ ``1``, ``chrM`` ↔ ``MT`` rules. Matrices are
    stored as ``matrix`` (``"cpk"``: CellPack, ``"mtx"``: Matrix Market,
    ``"csv"``: dense). ``reference`` is the genome the BAMs were aligned to
    (bgzip or plain FASTA with ``.fai`` / ``.gzi``, or 2bit), for labeling
    junction motifs and methylation calls; ``variants`` a BCF (with its
    CSI) of known sites, masked when counting conversions and methylation.
    """

    inputs: Sequence[PathLike]
    out_dir: PathLike
    products: Sequence[Product] = field(default_factory=tuple)
    annotation: Annotation | None = None
    aliases: Sequence[tuple[str, str]] | None = None
    matrix: Literal["cpk", "mtx", "csv"] = "cpk"
    reference: PathLike | None = None
    variants: PathLike | None = None

    def _native(self) -> tuple[tuple[Any, ...], dict[str, Any]]:
        if isinstance(self.inputs, (str, os.PathLike)):
            raise TypeError("inputs is a sequence of paths")
        args = (
            [os.fspath(p) for p in self.inputs],
            os.fspath(self.out_dir),
            [p._native() for p in self.products],
        )
        kwargs = {
            "annotation": None if self.annotation is None else self.annotation._native,
            "aliases": None if self.aliases is None else [tuple(a) for a in self.aliases],
            "matrix": self.matrix,
            "reference": None if self.reference is None else os.fspath(self.reference),
            "variants": None if self.variants is None else os.fspath(self.variants),
        }
        return args, kwargs


@dataclass(frozen=True, slots=True)
class MemoryPlan:
    """How a task divides its memory budget, in bytes."""

    fixed: int
    """Feature tables, indexes, whitelists, the annotation."""
    pipeline: int
    """Blocks and record batches in flight."""
    buffer: int
    """Sort buffers and spools of the products, spilled to disk when full."""


@dataclass(frozen=True, slots=True)
class Plan:
    """A task checked against its budget."""

    threads: int
    memory_budget: int
    memory: MemoryPlan
    needs: tuple[str, ...]
    """What the products read from records (``cigar``, ``tags``, ...)."""
    stages: tuple[tuple[str, int], ...]
    """Each accumulator with the buffer bytes it is granted."""


@dataclass(frozen=True, slots=True)
class Report:
    """What a finished run wrote."""

    files: tuple[tuple[str, Path], ...]
    """``(kind, path)`` of every file: ``bigwig``, ``cpk``, ``mtx``,
    ``barcodes``, ``features``, ``fragments``, ``tbi``, ``json``, ..."""
    records: int
    """Placed records read, every input."""
    bytes_read: int
    seconds: float
    summary: dict[str, Any]
    """The content of ``summary.json``."""


def plan(task: Task, *, threads: int | None = None, memory: int | str = "16G") -> Plan:
    """Checks ``task`` against a budget of ``threads`` (default: every core)
    and ``memory``, reading only the BAM headers (and building the feature
    spaces); raises when the smallest working set does not fit."""
    args, kwargs = task._native()
    t, budget, fixed, pipeline, buffer, needs, stages = _native.plan_quant(
        *args, **kwargs, threads=threads, memory=str(memory)
    )
    return Plan(
        t,
        budget,
        MemoryPlan(fixed, pipeline, buffer),
        tuple(needs),
        tuple((n, b) for n, b in stages),
    )


def run(
    task: Task,
    *,
    threads: int | None = None,
    memory: int | str = "16G",
    tmp_dir: PathLike | None = None,
    progress: Callable[[str, int, int], None] | None = None,
) -> Report:
    """Runs ``task``: one pass over each input, then every product is
    written; the files appear only once all of them are.

    Products spill to ``tmp_dir`` (default: the output directory) when they
    outgrow their share of ``memory``. The output does not depend on
    ``threads``. ``progress`` is called now and then with (phase, bytes
    read, bytes total).
    """
    args, kwargs = task._native()
    files, records, bytes_read, seconds, summary = _native.quantify(
        *args,
        **kwargs,
        threads=threads,
        memory=str(memory),
        tmp_dir=None if tmp_dir is None else os.fspath(tmp_dir),
        progress=progress,
    )
    return Report(
        tuple((k, Path(p)) for k, p in files),
        records,
        bytes_read,
        seconds,
        json.loads(summary),
    )
