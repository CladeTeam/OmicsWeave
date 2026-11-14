"""Exceptions: one class per kind of failure every reader and writer reports.

All derive from :class:`OmicsWeaveError` and carry the failure's details as
attributes. Each also derives from the built-in exception closest in meaning
(``ValueError``, ``LookupError``, ``OSError``, ...), so code that catches
those keeps working.
"""

from __future__ import annotations

import errno as _errno
from typing import Any

__all__ = [
    "BadMagicError",
    "BigEndianError",
    "FormatError",
    "InflateError",
    "InvalidError",
    "InvalidRegionError",
    "IoError",
    "NotQueriedError",
    "OmicsWeaveError",
    "OutOfRangeError",
    "OutputSizeError",
    "RegionOutOfRangeError",
    "StaleSidecarError",
    "TooDeepError",
    "TruncatedError",
    "UnknownChromError",
    "UnsupportedError",
]


class OmicsWeaveError(Exception):
    """Base of every OmicsWeave exception; details are keyword attributes."""

    def __init__(self, message: str, **details: Any) -> None:
        super().__init__(message)
        for name, value in details.items():
            setattr(self, name, value)


class TruncatedError(OmicsWeaveError):
    """A buffer or file ends before the structure it holds (``what``, ``need``, ``have``)."""

    what: str
    need: int
    have: int


class BadMagicError(OmicsWeaveError):
    """A magic number does not identify the expected structure (``what``, ``found``)."""

    what: str
    found: int


class UnsupportedError(OmicsWeaveError, NotImplementedError):
    """A valid input this implementation does not handle (``feature``)."""

    feature: str


class BigEndianError(UnsupportedError):
    """A byte-swapped (big-endian) file."""

    feature = "big-endian files"


class OutOfRangeError(OmicsWeaveError):
    """An offset outside the bytes it must address (``what``, ``offset``, ``lo``, ``hi``)."""

    what: str
    offset: int
    lo: int
    hi: int


class InvalidError(OmicsWeaveError, ValueError):
    """A structurally invalid field, record or argument (``what``, ``detail``)."""

    what: str
    detail: str


class TooDeepError(OmicsWeaveError):
    """A tree deeper than any well-formed file has, or a cycle (``depth``)."""

    depth: int


class FormatError(OmicsWeaveError, ValueError):
    """Text that cannot be parsed, at 1-based ``line`` (``detail``)."""

    line: int
    detail: str


class InflateError(OmicsWeaveError):
    """Compressed data that does not decompress or check out (``detail``)."""

    detail: str


class UnknownChromError(OmicsWeaveError, LookupError):
    """The source has no chromosome named ``chrom``."""

    chrom: str


class InvalidRegionError(OmicsWeaveError, ValueError):
    """A region with ``start > end``."""

    start: int
    end: int


class RegionOutOfRangeError(OmicsWeaveError, IndexError):
    """A region ending past its sequence of ``length`` bases where output must cover it."""

    start: int
    end: int
    length: int


class OutputSizeError(OmicsWeaveError, ValueError):
    """An output array of the wrong size (``what``, ``have``, ``need``)."""

    what: str
    have: int
    need: int


class NotQueriedError(OmicsWeaveError):
    """Output was requested before a query staged it."""


class StaleSidecarError(OmicsWeaveError):
    """A sidecar built for a different source; rebuild it from the source."""


_ERRNO = {
    "NotFound": _errno.ENOENT,
    "PermissionDenied": _errno.EACCES,
    "AlreadyExists": _errno.EEXIST,
    "IsADirectory": _errno.EISDIR,
    "NotADirectory": _errno.ENOTDIR,
    "UnexpectedEof": _errno.EIO,
    "StorageFull": _errno.ENOSPC,
    "ReadOnlyFilesystem": _errno.EROFS,
}


class IoError(OmicsWeaveError, OSError):
    """An I/O failure of the byte source (``kind`` names Rust's ``io::ErrorKind``)."""

    kind: str
    detail: str

    def __init__(self, message: str, **details: Any) -> None:
        super().__init__(message, **details)
        self.errno = _ERRNO.get(details.get("kind", ""))
        self.strerror = details.get("detail", message)
