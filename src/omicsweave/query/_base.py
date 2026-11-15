"""What every reader shares: it is opened on files, closed explicitly or on
leaving a ``with`` block, and pickled as the files it was opened on."""

from __future__ import annotations

from typing import Any, Self

from ..cache import CacheNamespace
from ..errors import InvalidError
from ..formats import File, PathLike, detect
from ..formats._base import ensure_format

__all__ = ["Reader"]


def _reopen(cls: type, args: tuple[Any, ...], kwargs: dict[str, Any]) -> Any:
    return cls(*args, **kwargs)


def resolve[F: File](source: F | PathLike, accepted: tuple[type[F], ...], reader: str) -> F:
    """``source`` as one of the ``accepted`` formats: a path is detected, a
    file object is checked against its content."""
    names = ", ".join(cls.__name__ for cls in accepted)
    if isinstance(source, File):
        if not isinstance(source, accepted):
            raise TypeError(f"{reader} reads {names}, not {type(source).__name__}")
        ensure_format(source)
        return source
    file = detect(source)
    if not isinstance(file, accepted):
        detail = f"{reader} reads {names}; this is {type(file).__name__}"
        raise InvalidError(f"{file.path}: {detail}", what=str(file.path), detail=detail)
    return file


class Reader:
    """An open reader: a context manager, closed explicitly or on exit, and
    pickled as the arguments it was opened with (it reopens on unpickling,
    e.g. in a DataLoader worker). Decoding releases the GIL, so threads can
    share one reader.

    Queries allocate result arrays automatically. Results remain valid
    after further queries or after this reader closes. Fixed-size queries
    optionally accept ``out=`` to write into caller-owned storage; those
    results share that storage and reflect any subsequent writes to it.
    Decoding workspaces are managed internally.
    """

    _native: Any
    _args: tuple[Any, ...]
    _kwargs: dict[str, Any]

    def _open(self, native: Any, *args: Any, **kwargs: Any) -> None:
        self._native = native
        self._args = args
        self._kwargs = kwargs
        self._closed = False

    def close(self) -> None:
        if not self._closed:
            self._closed = True
            self._native.close()

    @property
    def closed(self) -> bool:
        return self._closed

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __reduce__(self) -> tuple[Any, ...]:
        return (_reopen, (type(self), self._args, self._kwargs))

    def __repr__(self) -> str:
        args = ", ".join(repr(a) for a in getattr(self, "_args", ()))
        state = " (closed)" if getattr(self, "_closed", False) else ""
        return f"{type(self).__name__}({args}){state}"


class CachedReader(Reader):
    """A reader able to describe and reuse independently decoded blocks."""

    @property
    def cache_namespace(self) -> CacheNamespace:
        """Process-local source session identity, shared by execution clones."""
        return self._native.cache_namespace()

    def clone(self) -> Self:
        """Share the immutable source and cache with independent execution state."""
        native = self._native.clone_context()
        cloned = object.__new__(type(self))
        cloned.__dict__ = self.__dict__.copy()
        cloned._native = native
        cloned._closed = False
        return cloned
