"""Planning and running a conversion."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Literal

from .. import _budget, _native
from .._budget import Budget, MemoryPlan
from ..formats import BbiOptions, File, PathLike
from ..formats._bbi import native_options
from ._files import file_of, native_files

__all__ = ["ConvertWarning", "Plan", "Report", "plan", "run"]

Progress = Callable[[str, int, int], None]
"""Called now and then with (phase, bytes read, bytes total)."""


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


class Plan:
    """A conversion configured and checked against its budget; nothing has
    been read. :meth:`run` converts."""

    __slots__ = ("_native", "source", "target")

    def __init__(self, native: Any) -> None:
        self._native = native
        (src_format, dst_format), (src, dst) = native.formats, native.paths
        self.source: File = file_of(src_format, src)
        self.target: File = file_of(dst_format, dst)

    @property
    def budget(self) -> Budget:
        return _budget.budget(self._native.budget)

    @property
    def memory(self) -> MemoryPlan:
        return _budget.memory_plan(self._native.memory)

    def run(self, *, progress: Progress | None = None) -> Report:
        """Converts; the target and its companions appear only once the
        whole conversion succeeded. What it passed over is issued as
        :class:`ConvertWarning`."""
        written, records, read, spilled, seconds, notes = self._native.run(progress)
        target, *companions = (file_of(name, path) for name, path in written)
        for note in notes:
            warnings.warn(note, ConvertWarning, stacklevel=2)
        return Report(
            self.source, target, tuple(companions), records, read, spilled, seconds, tuple(notes)
        )

    def __reduce__(self) -> Any:
        raise TypeError("a conversion plan is process-local; plan the conversion again")

    def __repr__(self) -> str:
        return f"<convert.Plan {self.source!r} -> {self.target!r}>"


def plan(
    source: File | PathLike,
    target: File | PathLike,
    *,
    threads: int | None = None,
    memory: int | str | None = None,
    tmp_dir: PathLike | None = None,
    level: int | None = None,
    index: bool | Literal["auto", "bai", "csi"] = False,
    min_shift: int | None = None,
    chrom_sizes: PathLike | None = None,
    resolutions: Sequence[int] | None = None,
    value: str | None = None,
    mod: str | None = None,
    strand: Literal["+", "-", "."] | None = None,
    sidecar: bool = False,
    add_contigs: bool = False,
    sort: bool = False,
    reference: PathLike | None = None,
    verify_checksums: bool = True,
    bbi: BbiOptions | None = None,
) -> Plan:
    """Configures the conversion of ``source`` into ``target`` and checks it
    against its budget without reading any data.

    A path's format is told by content for the source and by name for the
    target; a file object (:mod:`omicsweave.formats`) states it.

    ``threads`` (default: every core the process may run on) run the
    parallel stages; ``memory`` (bytes, or a size such as ``"8G"``; default
    4 GiB) bounds the buffers. Memory that is not buffers (allocator slack,
    stacks) can take the process a little over. Sorts that outgrow the
    budget spill to ``tmp_dir`` (default: the target's directory).

    - ``level``: compression level (bgzip, BCF, BAM: 0-12; bigBed, bigWig,
      hic: zlib 1-9).
    - ``index``: index a sorted BCF (``True``: ``.csi``) or BAM (``True`` /
      ``"auto"``: ``.bai`` when every reference fits, else ``.csi``;
      ``"bai"``; ``"csi"``, with bins of ``2 ** min_shift`` bases, default
      14).
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
      ``verify_checksums=False`` skips the block CRC32 and slice MD5 checks.
    - ``bbi``: bigBed / bigWig block, index and zoom layout
      (:class:`~omicsweave.formats.BbiOptions`); ``level`` still sets the
      compression level when given.

    An option the conversion does not have raises
    :class:`~omicsweave.errors.InvalidError`.
    """
    native = _native.ConvertPlan(
        **native_files(source, target),
        threads=_budget.threads(threads),
        memory=_budget.memory(memory),
        tmp_dir=_budget.path(tmp_dir),
        level=level,
        index="auto" if index is True else (index or None),
        min_shift=min_shift,
        chrom_sizes=_budget.path(chrom_sizes),
        resolutions=None if resolutions is None else list(resolutions),
        value=value,
        mod=mod,
        strand=strand,
        sidecar=sidecar,
        add_contigs=add_contigs,
        sort=sort,
        reference=_budget.path(reference),
        verify_checksums=verify_checksums,
        bbi=None if bbi is None else native_options(bbi),
    )
    return Plan(native)


def run(
    source: File | PathLike,
    target: File | PathLike,
    *,
    progress: Progress | None = None,
    **options: Any,
) -> Report:
    """Converts ``source`` into ``target``: :func:`plan` with ``options``,
    then :meth:`Plan.run`."""
    return plan(source, target, **options).run(progress=progress)
