"""A quantification task, its plan, and what a run wrote."""

from __future__ import annotations

import json
import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, NamedTuple

from .. import _budget, _native
from .._budget import Budget, MemoryPlan
from ..formats import PathLike
from ._products import Product

if TYPE_CHECKING:
    from ..annotation import Annotation

__all__ = ["Input", "OutputFile", "Plan", "Report", "Task", "plan", "run"]

Progress = Callable[[str, int, int], None]
"""Called now and then with (phase, bytes read, bytes total)."""


@dataclass(frozen=True, slots=True)
class Task:
    """Coordinate-sorted BAMs, the products of one pass over them, and the
    directory they go to.

    ``annotation`` is needed by gene products (and the ATAC gene activity
    and TSS metrics); ``aliases`` match chromosome names that differ beyond
    the built-in ``chr1`` ↔ ``1``, ``chrM`` ↔ ``MT`` rules (annotation name
    → BAM name). Matrices are stored as ``matrix`` (``"cpk"``: CellPack,
    ``"mtx"``: Matrix Market, ``"csv"``: dense). ``reference`` is the genome
    the BAMs were aligned to (bgzip or plain FASTA with ``.fai`` /
    ``.gzi``, or 2bit), for labeling, junction motifs and methylation
    calls; ``variants`` a BCF (with its CSI) of known sites, masked when
    counting conversions and methylation.
    """

    inputs: Sequence[PathLike]
    out_dir: PathLike
    products: Sequence[Product] = field(default_factory=tuple)
    _: KW_ONLY
    annotation: Annotation | None = None
    aliases: Mapping[str, str] | None = None
    matrix: Literal["cpk", "mtx", "csv"] = "cpk"
    reference: PathLike | None = None
    variants: PathLike | None = None

    def __post_init__(self) -> None:
        if isinstance(self.inputs, (str, os.PathLike)):
            raise TypeError("inputs is a sequence of paths, not one path")
        object.__setattr__(self, "inputs", tuple(self.inputs))
        object.__setattr__(self, "products", tuple(self.products))


@dataclass(frozen=True, slots=True)
class Input:
    """An input BAM as its header describes it."""

    path: Path
    references: tuple[tuple[str, int], ...]
    """``(name, length)`` by reference id."""
    read_groups: tuple[str, ...]
    """The ``ID`` of each ``@RG`` line, in header order."""
    size: int
    """Compressed bytes."""


class OutputFile(NamedTuple):
    """A file a run wrote, with its kind (``bigwig``, ``cpk``, ``mtx``,
    ``barcodes``, ``features``, ``fragments``, ``tbi``, ``json``, ...)."""

    kind: str
    path: Path


@dataclass(frozen=True, slots=True)
class Report:
    """What a finished run wrote."""

    files: tuple[OutputFile, ...]
    records: int
    """Placed records read, every input."""
    bytes_read: int
    seconds: float
    scan_seconds: float
    """Of which the pass over the inputs; the rest writes the products."""
    summary: dict[str, Any]
    """The content of ``summary.json``."""


class Plan:
    """A task checked against its budget: the BAM headers read and the
    feature spaces built. :meth:`run` makes the pass, as often as asked."""

    __slots__ = ("_native", "task")

    def __init__(self, task: Task, native: Any) -> None:
        self.task = task
        self._native = native

    @property
    def budget(self) -> Budget:
        return _budget.budget(self._native.budget)

    @property
    def memory(self) -> MemoryPlan:
        return _budget.memory_plan(self._native.memory)

    @property
    def needs(self) -> tuple[str, ...]:
        """What the products read from records (``cigar``, ``tags``, ...)."""
        return tuple(self._native.needs)

    @property
    def stages(self) -> tuple[tuple[str, int], ...]:
        """Each accumulator with the buffer bytes it is granted."""
        return tuple((n, b) for n, b in self._native.stages)

    @property
    def inputs(self) -> tuple[Input, ...]:
        return tuple(
            Input(Path(p), tuple(refs), tuple(groups), size)
            for p, refs, groups, size in self._native.inputs
        )

    def run(self, *, progress: Progress | None = None) -> Report:
        """One pass over each input, then every product is written; the
        files appear only once all of them are. The output does not depend
        on the number of threads."""
        files, records, bytes_read, seconds, scan, summary = self._native.run(progress)
        return Report(
            tuple(OutputFile(k, Path(p)) for k, p in files),
            records,
            bytes_read,
            seconds,
            scan,
            json.loads(summary),
        )

    def __reduce__(self) -> Any:
        raise TypeError("a quantification plan is process-local; plan the task again")

    def __repr__(self) -> str:
        b = self.budget
        return f"<quant.Plan {len(self.task.inputs)} inputs, {b.threads} threads, {b.memory} B>"


def plan(
    task: Task,
    *,
    threads: int | None = None,
    memory: int | str | None = None,
    tmp_dir: PathLike | None = None,
) -> Plan:
    """Checks ``task`` against a budget of ``threads`` (default: every core
    the process may run on) and ``memory`` (bytes or a size such as
    ``"16G"``; default 4 GiB), reading only the BAM headers (and building
    the feature spaces); raises when the smallest working set does not fit.
    Products spill to ``tmp_dir`` (default: the output directory) when they
    outgrow their share of memory."""
    native = _native.QuantPlan(
        task,
        threads=_budget.threads(threads),
        memory=_budget.memory(memory),
        tmp_dir=_budget.path(tmp_dir),
    )
    return Plan(task, native)


def run(
    task: Task,
    *,
    threads: int | None = None,
    memory: int | str | None = None,
    tmp_dir: PathLike | None = None,
    progress: Progress | None = None,
) -> Report:
    """Plans ``task`` (see :func:`plan`) and runs it."""
    return plan(task, threads=threads, memory=memory, tmp_dir=tmp_dir).run(progress=progress)
