"""Running a conversion."""

from __future__ import annotations

import os
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Literal

from .. import _native
from ..formats import File, PathLike
from ._files import file_of, native_files


class ConvertWarning(UserWarning):
    """Something a conversion passed over without failing, such as rows on
    chromosomes outside the chromosome list."""


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
    warnings: tuple[str, ...] = ()
    """What the conversion passed over, each also issued as a
    :class:`ConvertWarning`."""

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
    value: str | None = None,
    mod: str | None = None,
    strand: Literal["+", "-", "."] | None = None,
    sidecar: bool = True,
    add_contigs: bool = False,
    sort: bool = False,
    reference: PathLike | None = None,
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

    - ``level``: compression level (bgzip, BCF, BAM: 0-12; bigBed, bigWig,
      hic: zlib 1-9).
    - ``index``: index a sorted BCF (``True``: ``.csi``) or BAM (``True`` /
      ``"auto"``: ``.bai`` when every reference fits, else ``.csi``;
      ``"bai"``; ``"csi"``).
    - ``chrom_sizes``: bigBed / bigWig chromosome lengths (``chrom.sizes``,
      ``.fai`` or ``.2bit``); by default each chromosome's largest interval
      end. Pairs to hic: the chromosomes and their order, instead of the
      ``#chromsize`` header lines (Juicer short text needs it).
    - ``resolutions``: ``.mcool`` resolutions to convert (default: all);
      pairs to hic: the bin sizes, each a multiple of the finest (default:
      Juicer ``pre``'s 2.5M, 1M, 500k, 250k, 100k, 50k, 25k, 10k, 5k).
    - ``value``: bedMethyl to bigWig: the column stored, as written:
      ``percent`` (default; ``percent_modified``), ``coverage``
      (``valid_coverage``), ``score``, and modkit's ``nmod``,
      ``ncanonical``, ``nother_mod``, ``ndelete``, ``nfail``, ``ndiff``,
      ``nnocall``.
    - ``mod``: bedMethyl to bigWig: the items of this modification code
      (column 4: ``m``, ``h``, ``a``, ...); required when there are several.
    - ``strand``: bedMethyl to bigWig: the items of this strand (``+``,
      ``-``, or ``.`` for combined strands); by default all, and items of
      both strands that overlap are an error.
    - ``sidecar``: also write the bigBed / bigWig / hic ``.sidecar``.
    - ``add_contigs``: BCF: declare contigs the VCF header lacks (one more
      pass over the source); otherwise an undeclared contig is an error.
    - ``sort``: BAM to BAM: sort the records by coordinate as ``samtools
      sort`` does (the only conversion of a BAM into a BAM, so required for
      it); ``index`` then indexes the sorted BAM.
    - ``reference``: CRAM to BAM: the reference the CRAM was compressed
      against (FASTA with ``.fai``, bgzip FASTA with ``.fai`` and ``.gzi``,
      or 2bit); not needed for CRAMs with embedded or no reference.
    - ``progress``: called now and then with (phase, bytes read, bytes
      total).

    An option the conversion does not have raises
    :class:`~omicsweave.errors.InvalidError`. What the conversion passed
    over is issued as :class:`ConvertWarning` and listed in
    :attr:`Report.warnings`.
    """
    files = native_files(source, target)
    src, written, records, read, spilled, seconds, notes = _native.convert(
        **files,
        threads=threads,
        memory=str(memory),
        tmp_dir=None if tmp_dir is None else os.fspath(tmp_dir),
        level=level,
        index="auto" if index is True else (index or None),
        chrom_sizes=None if chrom_sizes is None else os.fspath(chrom_sizes),
        resolutions=None if resolutions is None else list(resolutions),
        value=value,
        mod=mod,
        strand=strand,
        sidecar=sidecar,
        add_contigs=add_contigs,
        sort=sort,
        reference=None if reference is None else os.fspath(reference),
        progress=progress,
    )
    target_file, *companions = (file_of(name, path) for name, path in written)
    for note in notes:
        warnings.warn(note, ConvertWarning, stacklevel=2)
    return Report(
        file_of(src, files["src"]),
        target_file,
        tuple(companions),
        records,
        read,
        spilled,
        seconds,
        tuple(notes),
    )
