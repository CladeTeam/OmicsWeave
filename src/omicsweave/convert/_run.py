"""Running a conversion."""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Literal

from .. import _native
from ..formats import File, PathLike
from ._files import file_of, native_files


@dataclass(frozen=True, slots=True)
class Report:
    """What a finished conversion wrote."""

    source: File
    target: File
    companions: tuple[File, ...]
    """Indexes and sidecars written with the target."""
    records: int
    """Sequences, variants, alignments, intervals or pixels converted."""
    bytes_read: int
    """Bytes of the source read."""
    spilled: int
    """Bytes written to temporary spill files."""
    seconds: float

    @property
    def files(self) -> tuple[File, ...]:
        """Every file written, the target first."""
        return (self.target, *self.companions)


def run(
    source: File | PathLike,
    target: File | PathLike,
    *,
    threads: int | None = None,
    memory: int | str = "4G",
    tmp_dir: PathLike | None = None,
    level: int | None = None,
    index: bool | Literal["auto", "bai", "csi"] = False,
    chrom_sizes: PathLike | None = None,
    resolutions: Sequence[int] | None = None,
    sidecar: bool = True,
    add_contigs: bool = False,
    progress: Callable[[str, int, int], None] | None = None,
) -> Report:
    """Converts ``source`` into ``target``.

    A path's format is told by content for the source and by name for the
    target; a file object (:mod:`omicsweave.formats`) states it.

    ``threads`` (default: every core) run the parallel stages; ``memory``
    (bytes, or a size such as ``"8G"``) bounds the buffers. Memory that is
    not buffers (allocator slack, stacks) can take the process a little
    over. Sorts that outgrow the budget spill to ``tmp_dir`` (default: the
    target's directory).

    - ``level``: compression level (bgzip, BCF, BAM: 0-12; bigBed, hic:
      zlib 1-9).
    - ``index``: index a sorted BCF (``True``: ``.csi``) or BAM (``True`` /
      ``"auto"``: ``.bai`` when every reference fits, else ``.csi``;
      ``"bai"``; ``"csi"``).
    - ``chrom_sizes``: bigBed chromosome lengths (``chrom.sizes``, ``.fai``
      or ``.2bit``); by default each chromosome's largest interval end.
    - ``resolutions``: ``.mcool`` resolutions to convert (default: all).
    - ``sidecar``: also write the bigBed / hic ``.sidecar``.
    - ``add_contigs``: BCF: declare contigs the VCF header lacks (one more
      pass over the source); otherwise an undeclared contig is an error.
    - ``progress``: called now and then with (phase, bytes read, bytes
      total).

    An option the conversion does not have raises
    :class:`~omicsweave.errors.InvalidError`.
    """
    files = native_files(source, target)
    src, written, records, read, spilled, seconds = _native.convert(
        **files,
        threads=threads,
        memory=str(memory),
        tmp_dir=None if tmp_dir is None else os.fspath(tmp_dir),
        level=level,
        index="auto" if index is True else (index or None),
        chrom_sizes=None if chrom_sizes is None else os.fspath(chrom_sizes),
        resolutions=None if resolutions is None else list(resolutions),
        sidecar=sidecar,
        add_contigs=add_contigs,
        progress=progress,
    )
    target_file, *companions = (file_of(name, path) for name, path in written)
    return Report(
        file_of(src, files["src"]),
        target_file,
        tuple(companions),
        records,
        read,
        spilled,
        seconds,
    )
