"""How sources and targets cross into the native layer and come back as files."""

from __future__ import annotations

import os
from collections.abc import Sequence

from ..errors import InvalidError
from ..formats import FORMATS, File, PathLike

Source = File | PathLike | Sequence[File | PathLike]
"""A file, or for motif collections a list of files."""


def native_files(source: Source, target: File | PathLike) -> dict[str, object]:
    """The paths of ``source`` and ``target``, and the formats file objects
    state; a bare path leaves its format to the native layer (content for
    the source, name for the target). Several sources form one motif
    collection."""
    sources = (
        list(source)
        if isinstance(source, Sequence) and not isinstance(source, (str, bytes))
        else [source]
    )
    if not sources:
        raise InvalidError("no source given", what="source", detail="no source given")
    formats = {s.format for s in sources if isinstance(s, File)}
    if len(formats) > 1:
        detail = f"sources state different formats: {sorted(formats)}"
        raise InvalidError(detail, what="source", detail=detail)
    return {
        "src": os.fspath(sources[0]),
        "more_src": [os.fspath(s) for s in sources[1:]],
        "dst": os.fspath(target),
        "src_format": formats.pop() if formats else None,
        "dst_format": target.format if isinstance(target, File) else None,
    }


def file_of(name: str, path: PathLike) -> File:
    """The file at ``path`` as the format called ``name``."""
    return FORMATS[name](path)
