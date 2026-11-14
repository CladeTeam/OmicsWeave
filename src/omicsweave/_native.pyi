"""The native extension (``crates/python``). Its API is internal: it takes and
returns only primitives and numpy arrays, and the ``omicsweave`` packages
build the public types on it.

Every reader and writer keeps its state behind a lock taken only with the GIL
released; decoding and encoding run without the GIL. Failures raise the
exceptions of :mod:`omicsweave.errors`. Paths are ``str``; positions are
0-based, half-open and fit in ``uint32``.
"""

from collections.abc import Callable
from typing import Final, Literal, final

import numpy as np
import numpy.typing as npt

__version__: Final[str]

U8 = npt.NDArray[np.uint8]
U32 = npt.NDArray[np.uint32]
I64 = npt.NDArray[np.int64]
F32 = npt.NDArray[np.float32]

Unit = Literal["BP", "FRAG"]
"""What a contact-map bin coordinate counts."""

Tile = tuple[str, int, int, str, int, int, int, Unit]
"""A contact-map window: (row chrom, row start, row end, column chrom,
column start, column end, resolution, unit)."""

Selector = tuple[Literal["all", "samples", "haplotype"], list[str], str | None]
"""Which alleles a consensus takes: ``("all", [], None)``, ``("samples",
names, None)`` or ``("haplotype", [sample], haplotype)``."""

Progress = Callable[[str, int, int], object]
"""Called with (phase, bytes read, bytes total) between batches."""

# Formats and companion files

def detect(path: str) -> str:
    """The name of the format the file at ``path`` holds, from its content:
    ``fasta``, ``gzip-fasta``, ``bgzip-fasta``, ``2bit``, ``vcf``,
    ``gzip-vcf``, ``bcf``, ``sam``, ``gzip-sam``, ``bam``, ``csi``, ``bai``,
    ``bed``, ``gzip-bed``, ``bigwig``, ``bigbed``, ``cool``, ``mcool``,
    ``hic`` or ``sidecar``."""

def build_fasta_index(path: str) -> list[str]:
    """Builds the ``.fai`` (and, for bgzip FASTA, the ``.gzi``) beside
    ``path`` in one pass; returns their paths, ``.fai`` first."""

def build_sidecar(path: str) -> str:
    """Builds the sidecar of a bigWig, bigBed or ``.hic`` file beside it,
    reading the whole file; returns its path."""

# Conversion

def convert(
    src: str,
    dst: str,
    *,
    src_format: str | None = None,
    dst_format: str | None = None,
    threads: int | None = None,
    memory: str = "4G",
    tmp_dir: str | None = None,
    level: int | None = None,
    index: Literal["auto", "bai", "csi"] | None = None,
    chrom_sizes: str | None = None,
    resolutions: list[int] | None = None,
    sidecar: bool = True,
    add_contigs: bool = False,
    progress: Progress | None = None,
) -> tuple[str, list[tuple[str, str]], int, int, int, float]:
    """Converts ``src`` into ``dst``. A format not named is told by content
    for the source and by name for the destination; ``threads`` defaults to
    every available core, ``memory`` is a size such as ``"8G"``. An option
    the conversion does not have is an error.

    Returns (source format, (format, path) of every file written with the
    destination first, records, bytes read, bytes spilled, seconds).
    """

def plan_convert(
    src: str,
    dst: str,
    *,
    src_format: str | None = None,
    dst_format: str | None = None,
    threads: int | None = None,
    memory: str = "4G",
) -> tuple[str, str, int, int, int, int, int]:
    """Checks the conversion of ``src`` into ``dst`` against its budget
    without reading any data; returns (source format, destination format,
    threads, memory budget, fixed, pipeline, buffer bytes)."""

def conversions() -> list[tuple[str, str]]:
    """Every (source, destination) pair of format names there is a
    conversion for."""

def parse_size(text: str) -> int:
    """``text`` (``512M``, ``8G``, ``1.5GiB``, a byte count) in bytes."""

def human_size(bytes: int) -> str:
    """``bytes`` for people: ``512 B``, ``1.5 GiB``."""

# Sequences: FASTA, bgzip FASTA, 2bit

@final
class SequenceReader:
    """A FASTA, bgzip FASTA or 2bit file, told by content. FASTA is read
    through its ``fai`` and bgzip FASTA also its ``gzi`` (by default the
    files beside it); they must exist. ``soft_mask`` lowercases the
    soft-masked bases of a 2bit file."""

    def __new__(
        cls,
        path: str,
        fai: str | None = None,
        gzi: str | None = None,
        soft_mask: bool = False,
    ) -> SequenceReader: ...
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length) of every sequence, in file order."""
    def sequence(self, chrom: str, start: int, end: int, out: U8) -> None:
        """Writes the ASCII bases of ``chrom:[start, end)``, which must lie
        within the sequence, into ``out`` (one element per base)."""
    def close(self) -> None:
        """Releases the file; later calls raise ``ValueError``."""

@final
class FastaWriter:
    """FASTA being written, its ``.fai`` (and ``.gzi``) built as it goes.
    ``level`` (0-12) compresses as BGZF; ``None`` writes plain text.
    ``line_width`` bases go on each line (0: one line per sequence)."""

    def __new__(cls, path: str, line_width: int = 60, level: int | None = None) -> FastaWriter: ...
    def add(self, name: str, seq: bytes) -> None:
        """Appends sequence ``name``; names are unique and bases printable ASCII."""
    def close(self) -> list[str]:
        """Ends the file and saves its indexes beside it; returns their
        paths, ``.fai`` first (``[]`` once closed)."""

# Signal: bigWig

@final
class TrackReader:
    """A bigWig file."""

    def __new__(cls, path: str) -> TrackReader: ...
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length) of every chromosome."""
    def summary(self) -> tuple[int, float, float, float, float] | None:
        """(bases covered, min, max, sum, sum of squares), as the header
        stores them."""
    def values(self, chrom: str, start: int, end: int, out: F32) -> None:
        """Writes one value per base of ``chrom:[start, end)`` into ``out``;
        NaN where there is no data. The range may run past the chromosome."""
    def intervals(self, chrom: str, start: int, end: int) -> tuple[U32, U32, F32]:
        """The intervals with data in ``chrom:[start, end)``, clipped to it,
        as (start, end, value) columns."""
    def close(self) -> None:
        """Releases the file; later calls raise ``ValueError``."""

@final
class BigWigWriter:
    """A bigWig file being written over ``chroms`` ((name, length) in file
    order) with the BBI options of :class:`omicsweave.formats.BbiOptions`."""

    def __new__(
        cls,
        path: str,
        chroms: list[tuple[str, int]],
        items_per_slot: int | None,
        block_size: int,
        compress: bool,
        level: int,
        zoom_levels: int,
        zoom_reductions: list[int] | None,
    ) -> BigWigWriter: ...
    def add(self, chrom: int, start: U32, end: U32, value: F32) -> None:
        """Adds items of chromosome ``chrom`` (its index in ``chroms``),
        sorted and not overlapping; the columns have one length."""
    def close(self, sidecar: bool) -> str | None:
        """Ends the file; with ``sidecar``, saves its sidecar beside it and
        returns its path."""

# Features: bigBed, BED

@final
class IntervalReader:
    """A bigBed file, or a BED file (plain or gzip) loaded whole; told by
    content."""

    auto_sql: str | None
    """The autoSql schema of a bigBed file."""
    field_count: tuple[int, int] | None
    """(fields, standard BED fields) of a bigBed file."""

    def __new__(cls, path: str) -> IntervalReader: ...
    def chroms(self) -> list[tuple[str, int | None]]:
        """(name, length) of every chromosome; BED files record no lengths."""
    def query(self, chrom: str, start: int, end: int) -> tuple[U32, U32, I64, U8]:
        """The features overlapping ``chrom:[start, end)``, unclipped, as
        columns: starts, ends, and the annotations (fields after the third)
        as Arrow ``LargeString`` offsets and bytes."""
    def close(self) -> None:
        """Releases the file; later calls raise ``ValueError``."""

@final
class BigBedWriter:
    """A bigBed file being written over ``chroms`` ((name, length) in file
    order), with records of ``field_count`` fields of which the first
    ``defined_field_count`` are standard BED, and the BBI options of
    :class:`omicsweave.formats.BbiOptions`."""

    def __new__(
        cls,
        path: str,
        chroms: list[tuple[str, int]],
        field_count: int,
        defined_field_count: int,
        auto_sql: str | None,
        items_per_slot: int | None,
        block_size: int,
        compress: bool,
        level: int,
        zoom_levels: int,
        zoom_reductions: list[int] | None,
    ) -> BigBedWriter: ...
    def add(self, chrom: int, start: U32, end: U32, offsets: I64, data: U8) -> None:
        """Adds features of chromosome ``chrom`` (its index in ``chroms``),
        sorted by start; the annotations come as Arrow ``LargeString``
        offsets (one more than features) and bytes. The whole batch is
        checked before any of it is written."""
    def close(self, sidecar: bool) -> str | None:
        """Ends the file; with ``sidecar``, saves its sidecar beside it and
        returns its path (``None`` once closed)."""

# Contacts: .hic

@final
class ContactReader:
    """A ``.hic`` file, versions 8 and 9: raw contact counts, never
    normalized."""

    version: int
    genome: str

    def __new__(cls, path: str) -> ContactReader: ...
    def attributes(self) -> list[tuple[str, str]]:
        """The header's (key, value) attributes."""
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length in base pairs) of every chromosome, in file order."""
    def resolutions(self, unit_name: Unit) -> list[int]:
        """The bin sizes stored in ``unit_name``."""
    def fragment_count(self, chrom: str) -> int:
        """The restriction fragments of ``chrom``."""
    def sites(self, chrom: str) -> U32:
        """The restriction site positions of ``chrom``."""
    def matrix(self, tile: Tile, out: npt.NDArray[np.float32]) -> None:
        """Writes the contacts of ``tile`` row-major into ``out``: a
        C-contiguous ``float32`` array of (row bins, column bins), 0 where
        there is no contact. ``out`` is left intact on error."""
    def cells(self, tile: Tile) -> tuple[U32, U32, F32]:
        """The nonzero pixels of ``tile``, ordered by (row, col), as (row
        bin start, column bin start, value) columns."""
    def close(self) -> None:
        """Releases the file; later calls raise ``ValueError``."""

@final
class HicWriter:
    """A ``.hic`` file being written from binned maps. ``sites`` holds the
    sorted restriction sites of each chromosome of ``chroms`` (empty without
    fragment resolutions); ``level`` is the zlib level (1-9)."""

    def __new__(
        cls,
        path: str,
        version: int,
        genome: str,
        attributes: list[tuple[str, str]],
        chroms: list[tuple[str, int]],
        bp_resolutions: list[int],
        frag_resolutions: list[int],
        sites: list[list[int]],
        level: int,
    ) -> HicWriter: ...
    def write_matrix(
        self,
        chr1: int,
        chr2: int,
        unit_name: Unit,
        bin_size: int,
        x: U32,
        y: U32,
        value: F32,
    ) -> None:
        """Writes the map of chromosomes ``chr1`` x ``chr2`` (indexes in
        ``chroms``) at ``bin_size``: cells at bin ``x`` of ``chr1`` and bin
        ``y`` of ``chr2``, in any order and orientation."""
    def close(self, sidecar: bool) -> str | None:
        """Ends the file; with ``sidecar``, saves its sidecar beside it and
        returns its path."""

@final
class HicBuilder:
    """A ``.hic`` file built from contacts, binned at every resolution when
    closed (as the Juicer tools' ``pre`` does). Arguments as for
    :class:`HicWriter`."""

    def __new__(
        cls,
        path: str,
        version: int,
        genome: str,
        attributes: list[tuple[str, str]],
        chroms: list[tuple[str, int]],
        bp_resolutions: list[int],
        frag_resolutions: list[int],
        sites: list[list[int]],
        level: int,
    ) -> HicBuilder: ...
    def add(self, chr1: U32, pos1: U32, chr2: U32, pos2: U32, count: F32) -> None:
        """Adds ``count[i]`` contacts between base ``pos1[i]`` of ``chr1[i]``
        and base ``pos2[i]`` of ``chr2[i]`` (indexes in ``chroms``)."""
    def close(self, sidecar: bool) -> str | None:
        """Bins every contact at every resolution and writes the file; with
        ``sidecar``, saves its sidecar beside it and returns its path."""

# Consensus: VCF/BCF over a reference

@final
class ConsensusBuilder:
    """Personal sequences from a VCF/BCF file (loaded once, whole) over
    ``reference``, as ``bcftools consensus`` 1.23 builds them. The options
    are ``bcftools consensus``'s: ``-I``, ``-M``, ``-a``, ``--mark-del``,
    ``--mark-ins``, ``--mark-snv``, ``-m``, ``--mask-with``,
    ``--regions-overlap``, ``-c``, ``-i`` and ``-e`` (``include`` and
    ``exclude`` are filter expressions; at most one is given)."""

    def __new__(
        cls,
        vcf: str,
        reference: SequenceReader,
        iupac_codes: bool,
        missing: str | None,
        absent: str | None,
        mark_del: str | None,
        mark_ins: str | None,
        mark_snv: str | None,
        mask: str | None,
        mask_with: str,
        regions_overlap: Literal["start", "record", "variant"],
        chain: bool,
        include: str | None,
        exclude: str | None,
    ) -> ConsensusBuilder: ...
    def samples(self) -> list[str]:
        """The samples of the VCF, in header order."""
    def records(self) -> tuple[int, int]:
        """(records kept, malformed lines skipped)."""
    def sequence(
        self, chrom: str, start: int, end: int, selector: Selector
    ) -> tuple[U8, str | None]:
        """The consensus of ``chrom:[start, end)`` with ``selector``: its
        ASCII bases, and its chain with ``chain``."""
    def sequences(
        self, chrom: str, start: int, end: int, selectors: list[Selector]
    ) -> tuple[I64, U8, list[str | None]]:
        """The consensus of ``chrom:[start, end)`` with each selector,
        packed: sequence ``i`` is ``bases[offsets[i]:offsets[i + 1]]``;
        returns (offsets, bases, chains)."""
    def close(self) -> None:
        """Releases the records; the reference stays open."""
