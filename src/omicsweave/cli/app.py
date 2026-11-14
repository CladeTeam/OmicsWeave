"""The ``ow`` application: the parser of every command line, and
:func:`main`, the entry point of the ``ow`` executable.

One subcommand per module; each registers its parser and the function that
runs it. Failures of the library print as ``ow: error: ...`` and exit 1;
usage errors exit 2, an interrupt 130.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from .._native import __version__
from ..errors import OmicsWeaveError
from . import _convert

__all__ = ["main"]

_COMMANDS = (_convert,)


def parser() -> argparse.ArgumentParser:
    """The parser of every ``ow`` command line."""
    root = argparse.ArgumentParser(
        prog="ow", description="OmicsWeave: multi-omics files for model training."
    )
    root.add_argument("--version", action="version", version=f"ow {__version__}")
    commands = root.add_subparsers(title="commands", metavar="COMMAND", required=True)
    for command in _COMMANDS:
        command.register(commands)
    return root


def main(argv: Sequence[str] | None = None) -> int:
    """Runs ``ow`` with ``argv`` (default: the process's arguments); returns
    the exit status."""
    args = parser().parse_args(argv)
    try:
        return args.run(args)
    except (OmicsWeaveError, OSError) as e:
        print(f"ow: error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print(file=sys.stderr)
        return 130
