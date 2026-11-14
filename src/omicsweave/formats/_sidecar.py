"""OmicsWeave sidecar (``.sidecar``): what a loader must know of a file
before reading any of its data blocks."""

from __future__ import annotations

import os

from .. import _native
from ._base import File, PathLike, companion_path, ensure_format

__all__ = ["Sidecar"]


class Sidecar(File):
    """A ``.sidecar`` of a bigWig, bigBed or ``.hic`` file: every data
    block's size once inflated, and what identifies the file it was built for.

    Only loaders use sidecars; reading by key does not.
    """

    __slots__ = ()
    format = "sidecar"
    suffix = ".sidecar"

    @classmethod
    def beside(cls, source: PathLike) -> Sidecar:
        """The sidecar named after ``source``: ``<source>.sidecar``."""
        return cls(companion_path(source, cls.suffix))


def build_sidecar(source: File) -> Sidecar:
    """Builds the sidecar of ``source`` beside it, reading the whole file."""
    ensure_format(source)
    return Sidecar(_native.build_sidecar(os.fspath(source)))


def written_sidecar(path: str | None) -> tuple[File, ...]:
    """The companion files of a writer that saves its sidecar at ``path``, if any."""
    return () if path is None else (Sidecar(path),)
