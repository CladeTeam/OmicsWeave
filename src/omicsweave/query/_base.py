"""What every reader shares: it is opened on files, closed explicitly or on
leaving a ``with`` block, and pickled as the files it was opened on."""

from __future__ import annotations

import os
from typing import Any, Self

from ..errors import InvalidError
from ..formats import (
    Bed,
    BedGraph,
    BedMethyl,
    Cool,
    File,
    GzipBed,
    GzipBedGraph,
    GzipBedMethyl,
    GzipFasta,
    GzipHicProPairs,
    GzipJuicerShort,
    GzipMaf,
    GzipPairs,
    GzipSam,
    GzipVcf,
    GzipWig,
    HicProPairs,
    JuicerShort,
    Maf,
    Mcool,
    Pairs,
    PathLike,
    Sam,
    Vcf,
    Wig,
    detect,
)
from ..formats._base import ensure_format

__all__ = ["Reader"]

_SIZES = " --chrom-sizes <genome.sizes>"
# What each format that is not read by region converts to: the target's
# format, its suffix, and the options the conversion usually needs.
_CONVERSIONS: dict[type[File], tuple[str, str, str]] = {
    Bed: ("bigBed", ".bb", _SIZES),
    GzipBed: ("bigBed", ".bb", _SIZES),
    BedGraph: ("bigWig", ".bw", _SIZES),
    GzipBedGraph: ("bigWig", ".bw", _SIZES),
    Wig: ("bigWig", ".bw", _SIZES),
    GzipWig: ("bigWig", ".bw", _SIZES),
    BedMethyl: ("bigBed", ".bb", ""),
    GzipBedMethyl: ("bigBed", ".bb", ""),
    Vcf: ("BCF", ".bcf", " --index"),
    GzipVcf: ("BCF", ".bcf", " --index"),
    Sam: ("BAM", ".bam", " --index"),
    GzipSam: ("BAM", ".bam", " --index"),
    GzipFasta: ("bgzip FASTA", ".bgz.fa.gz", ""),
    Cool: (".hic", ".hic", ""),
    Mcool: (".hic", ".hic", ""),
    Pairs: (".hic", ".hic", ""),
    GzipPairs: (".hic", ".hic", ""),
    JuicerShort: (".hic", ".hic", _SIZES),
    GzipJuicerShort: (".hic", ".hic", _SIZES),
    HicProPairs: (".hic", ".hic", _SIZES),
    GzipHicProPairs: (".hic", ".hic", _SIZES),
    Maf: ("bigMaf", ".bigmaf", ""),
    GzipMaf: ("bigMaf", ".bigmaf", ""),
}


def _stem(path: str) -> str:
    name = path[:-3] if path.endswith(".gz") else path
    root, ext = os.path.splitext(name)
    return root if ext else name


def conversion(file: File) -> str | None:
    """Why ``file`` is not read by region and the command that makes it
    readable; ``None`` when it is read by region."""
    known = _CONVERSIONS.get(type(file))
    if known is None:
        return None
    target, suffix, options = known
    path = os.fspath(file.path)
    return (
        f"{file.format} is not read by region; convert it to {target} with "
        f"`ow convert {path} {_stem(path)}{suffix}{options}`"
    )


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
    """A reader whose decoded blocks a :class:`~omicsweave.cache.BlockCache`
    can hold, and that clones for threads."""

    def clone(self) -> Self:
        """A reader of the same open file, index and cache with its own
        decoding state, for another thread: clones decode in parallel
        instead of taking turns."""
        native = self._native.clone_context()
        cloned = object.__new__(type(self))
        cloned.__dict__ = self.__dict__.copy()
        cloned._native = native
        cloned._closed = False
        return cloned
