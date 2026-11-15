"""What every file format shares: a file is its path, a writer finishes one."""

from __future__ import annotations

import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any, ClassVar, Self

from .. import _native
from ..errors import InvalidError, UnknownChromError

__all__ = ["File", "PathLike", "Writer"]

PathLike = str | os.PathLike[str]


class File(os.PathLike[str]):
    """A file of one format, named by its path.

    Creating one opens nothing: it states what the file at ``path`` holds.
    Files are path-like, so they go wherever a path does, and compare equal
    when they are the same format at the same path.
    """

    __slots__ = ("_path",)

    format: ClassVar[str]
    """The format's name; for the formats :func:`~omicsweave.formats.detect`
    tells apart, the name it gives."""

    def __init__(self, path: PathLike) -> None:
        self._path = Path(path)

    @property
    def path(self) -> Path:
        return self._path

    def exists(self) -> bool:
        return self._path.is_file()

    def __fspath__(self) -> str:
        return os.fspath(self._path)

    def __eq__(self, other: object) -> bool:
        if type(other) is not type(self):
            return NotImplemented
        return self._path == other._path  # type: ignore[attr-defined]

    def __hash__(self) -> int:
        return hash((type(self), self._path))

    def __reduce__(self) -> tuple[Any, ...]:
        return (type(self), (self._path,))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({str(self._path)!r})"


def companion_path(file: PathLike, suffix: str) -> Path:
    """The path of ``file`` with ``suffix`` appended, as companion files are named."""
    return Path(f"{os.fspath(file)}{suffix}")


class ChromIndex:
    """Chromosome numbers by name, in the order a header declares them."""

    __slots__ = ("_ids",)

    def __init__(self, chroms: Mapping[str, int]) -> None:
        self._ids = {name: i for i, name in enumerate(chroms)}

    def __getitem__(self, chrom: str) -> int:
        try:
            return self._ids[chrom]
        except KeyError:
            raise UnknownChromError(f"unknown chromosome {chrom!r}", chrom=chrom) from None


def ensure_format(file: File) -> None:
    """Checks that the content of ``file`` is its format."""
    found = _native.detect(os.fspath(file))
    if found != file.format:
        detail = f"holds {found} data, not {type(file).__name__}"
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)


class Writer[F: File]:
    """A file being written: finished by :meth:`close`, or on leaving a
    ``with`` block without an exception. Leaving with an exception abandons
    the file unfinished.

    Writers are process-local. After a normal CPython fork, the child must
    create its own writers. Closing or collecting an inherited buffered sink
    does not flush the parent's data; its inherited buffers and descriptors
    remain reserved in the child until process exit.
    """

    _native: Any
    file: F
    """The file being written."""

    def __init__(self, file: F, native: Any) -> None:
        self.file = file
        self._native = native
        self._closed = False
        self._process = os.getpid()

    def _finish(self) -> tuple[File, ...]:
        raise NotImplementedError

    def close(self) -> tuple[File, ...]:
        """Finishes the file; returns the companion files written with it."""
        if self._closed:
            return ()
        self._closed = True
        if self._process != os.getpid():
            # Foreign sinks are abandoned without running finish or flushing
            # the parent's buffered data. Their native drop policy is inert.
            self._native = None
            return ()
        return self._finish()

    @property
    def closed(self) -> bool:
        return self._closed

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type: object, *exc: object) -> None:
        if exc_type is None:
            self.close()
        else:
            self._closed = True
            self._native = None

    def __repr__(self) -> str:
        state = " (closed)" if getattr(self, "_closed", False) else ""
        return f"{type(self).__name__}({getattr(self, 'file', '?')!r}){state}"
