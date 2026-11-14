"""How sources and targets cross into the native layer and come back as files."""

from __future__ import annotations

import os

from ..formats import FORMATS, File, PathLike


def native_files(source: File | PathLike, target: File | PathLike) -> dict[str, str | None]:
    """The paths of ``source`` and ``target``, and the formats file objects
    state; a bare path leaves its format to the native layer (content for
    the source, name for the target)."""
    return {
        "src": os.fspath(source),
        "dst": os.fspath(target),
        "src_format": source.format if isinstance(source, File) else None,
        "dst_format": target.format if isinstance(target, File) else None,
    }


def file_of(name: str, path: PathLike) -> File:
    """The file at ``path`` as the format called ``name``."""
    return FORMATS[name](path)
