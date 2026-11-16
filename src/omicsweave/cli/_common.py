"""What the commands share: budget arguments, progress, sizes for people."""

from __future__ import annotations

import argparse
import sys
from typing import Any, TextIO

from .. import _native
from ..errors import OmicsWeaveError

__all__ = ["Progress", "add_budget", "human", "progress_for"]

human = _native.human_size


def _size(text: str) -> int:
    try:
        size = _native.parse_size(text)
    except OmicsWeaveError:
        size = 0
    if size == 0:
        raise argparse.ArgumentTypeError(f"not a memory budget: {text!r} (e.g. 512M, 8G)")
    return size


def _threads(text: str) -> int:
    try:
        n = int(text)
    except ValueError:
        n = 0
    if not 1 <= n <= sys.maxsize:
        raise argparse.ArgumentTypeError(f"not a thread count: {text!r} (at least 1)")
    return n


def add_budget(p: argparse.ArgumentParser, spill: str) -> None:
    """``-t``, ``-m`` and ``--tmp-dir``; unset, they leave the budget to the
    library (every core the process may use, 4 GiB)."""
    budget = p.add_argument_group("budget")
    budget.add_argument(
        "-t",
        "--threads",
        type=_threads,
        metavar="N",
        help="threads for the parallel stages (default: every core this process may use)",
    )
    budget.add_argument(
        "-m",
        "--memory",
        type=_size,
        metavar="SIZE",
        help="memory for buffers, e.g. 512M, 8G (default: 4G)",
    )
    budget.add_argument("--tmp-dir", metavar="DIR", help=f"where buffers spill (default: {spill})")


def budget_of(args: argparse.Namespace) -> dict[str, Any]:
    return {"threads": args.threads, "memory": args.memory, "tmp_dir": args.tmp_dir}


class Progress:
    """One line on ``stream``, rewritten as the run goes."""

    def __init__(self, stream: TextIO) -> None:
        self.stream = stream
        self.shown = False

    def __call__(self, phase: str, read: int, total: int) -> None:
        share = f" {100 * read / total:5.1f}%" if total else ""
        self.stream.write(f"\r{phase:<8} {human(read)} / {human(total)}{share}\x1b[K")
        self.stream.flush()
        self.shown = True

    def clear(self) -> None:
        if self.shown:
            self.stream.write("\r\x1b[K")
            self.stream.flush()


def progress_for(args: argparse.Namespace) -> Progress | None:
    return None if args.quiet or not sys.stderr.isatty() else Progress(sys.stderr)


def print_budget(budget: Any, memory: Any, buffer: str) -> None:
    print(
        f"memory  fixed {human(memory.fixed)}, pipeline {human(memory.pipeline)}, "
        f"{buffer} {human(memory.buffer)} (budget {human(budget.memory)}, "
        f"{budget.threads} threads)"
    )
