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
I8 = npt.NDArray[np.int8]
BOOL = npt.NDArray[np.bool_]
U16 = npt.NDArray[np.uint16]
U32 = npt.NDArray[np.uint32]
I32 = npt.NDArray[np.int32]
I64 = npt.NDArray[np.int64]
F32 = npt.NDArray[np.float32]

Unit = Literal["BP", "FRAG"]
"""What a contact-map bin coordinate counts."""

LocusKey = tuple[str, list[tuple[int, int]], str]
"""A dense-query key: chromosome, ascending disjoint ``(start, end)`` spans
and strand (``+ - . ?``). Output runs 5' to 3' of the strand over the
concatenated spans."""

Tile = tuple[str, int, int, str, int, int, int, Unit]
"""A contact-map window: (row chrom, row start, row end, column chrom,
column start, column end, resolution, unit)."""

Selector = tuple[Literal["all", "samples", "haplotype"], list[str], str | None]
"""Which alleles a consensus takes: ``("all", [], None)``, ``("samples",
names, None)`` or ``("haplotype", [sample], haplotype)``."""

Progress = Callable[[str, int, int], object]
"""Called with (phase, bytes read, bytes total) between batches."""

Strings = tuple[I64, U8] | None
"""A string column (offsets, data), or ``None`` when not asked for."""

@final
class CacheNamespace:
    """Opaque process-local identity of one immutable source session."""
    def __hash__(self) -> int: ...

@final
class BlockKey:
    """Opaque source-bound block identity produced by a reader."""
    @property
    def source(self) -> CacheNamespace: ...
    def __hash__(self) -> int: ...

@final
class BlockCache:
    def __new__(cls, max_bytes: int, max_entries: int) -> BlockCache: ...
    def stats(self) -> dict[str, int]: ...
    def keys(self) -> list[BlockKey]: ...
    def evict(self, key: BlockKey) -> bool: ...
    def clear(self, source: CacheNamespace | None = None) -> None: ...
    def resize(self, max_bytes: int, max_entries: int) -> None: ...

# Formats and companion files

def detect(path: str) -> str:
    """The name of the format the file at ``path`` holds, from its content:
    ``fasta``, ``gzip-fasta``, ``bgzip-fasta``, ``2bit``, ``vcf``,
    ``gzip-vcf``, ``bcf``, ``sam``, ``gzip-sam``, ``bam``, ``csi``, ``bai``,
    ``bed``, ``gzip-bed``, ``bedgraph``, ``gzip-bedgraph``, ``wig``,
    ``gzip-wig``, ``bigwig``, ``bigbed``, ``cool``, ``mcool``, ``hic``,
    ``sidecar``, ``pairs``, ``gzip-pairs``, ``juicer-short``,
    ``gzip-juicer-short``, ``bedmethyl`` or ``gzip-bedmethyl``."""

def build_fasta_index(path: str) -> list[str]:
    """Builds the ``.fai`` (and, for bgzip FASTA, the ``.gzi``) beside
    ``path`` in one pass; returns their paths, ``.fai`` first."""

def build_sidecar(path: str) -> str:
    """Builds the sidecar of a bigWig, bigBed or ``.hic`` file beside it,
    reading the whole file; returns its path."""

def build_bam_index(
    path: str, kind: Literal["auto", "bai", "csi"], min_shift: int, threads: int
) -> tuple[Literal["bai", "csi"], str]:
    """Builds the index of the sorted BAM at ``path`` beside it (``.bai`` or
    ``.csi``) as ``samtools index`` builds it, inflating on ``threads``;
    returns its kind and path."""

def build_cram_index(path: str) -> str:
    """Builds the ``.crai`` of the sorted CRAM at ``path`` beside it as
    ``samtools index`` builds it; returns its path."""

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
    value: str | None = None,
    mod: str | None = None,
    strand: Literal["+", "-", "."] | None = None,
    sidecar: bool = True,
    add_contigs: bool = False,
    sort: bool = False,
    reference: str | None = None,
    progress: Progress | None = None,
) -> tuple[str, list[tuple[str, str]], int, int, int, float, list[str]]:
    """Converts ``src`` into ``dst``. A format not named is told by content
    for the source and by name for the destination; ``threads`` defaults to
    every available core, ``memory`` is a size such as ``"8G"``. An option
    the conversion does not have is an error.

    Returns (source format, (format, path) of every file written with the
    destination first, records, bytes read, bytes spilled, seconds,
    warnings).
    """

def plan_convert(
    src: str,
    dst: str,
    *,
    src_format: str | None = None,
    dst_format: str | None = None,
    threads: int | None = None,
    memory: str = "4G",
    chrom_sizes: str | None = None,
    resolutions: list[int] | None = None,
    value: str | None = None,
    mod: str | None = None,
    strand: Literal["+", "-", "."] | None = None,
    sort: bool = False,
) -> tuple[str, str, int, int, int, int, int]:
    """Checks the conversion of ``src`` into ``dst`` (sorting, given
    ``sort``; with the options that shape its plan) against its budget
    without reading any data; returns (source
    format, destination format, threads, memory budget, fixed, pipeline,
    buffer bytes)."""

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
        cache: BlockCache | None = None,
    ) -> SequenceReader: ...
    def fetch(self, key: LocusKey) -> U8: ...
    def set_cache(self, cache: BlockCache | None) -> None: ...
    def cache_namespace(self) -> CacheNamespace: ...
    def describe_cache(self, key: LocusKey) -> list[BlockKey]: ...
    def prefetch_cache(self, key: LocusKey) -> None: ...
    def clone_context(self) -> SequenceReader: ...
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length) of every sequence, in file order."""
    def sequence(self, key: LocusKey, out: U8) -> None:
        """Writes the ASCII bases of ``key``, whose spans must lie within the
        sequence, into ``out`` (one element per position); reverse
        complemented on ``-``."""
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

    def __new__(cls, path: str, cache: BlockCache | None = None) -> TrackReader: ...
    def set_cache(self, cache: BlockCache | None) -> None: ...
    def cache_namespace(self) -> CacheNamespace: ...
    def describe_cache(self, key: LocusKey) -> list[BlockKey]: ...
    def prefetch_cache(self, key: LocusKey) -> None: ...
    def clone_context(self) -> TrackReader: ...
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length) of every chromosome."""
    def summary(self) -> tuple[int, float, float, float, float] | None:
        """(bases covered, min, max, sum, sum of squares), as the header
        stores them."""
    def values(self, key: LocusKey, out: F32) -> None:
        """Writes one value per position of ``key`` into ``out``, 5' to 3';
        NaN where there is no data. Spans may run past the chromosome."""
    def bins(
        self,
        key: LocusKey,
        resolution: int,
        phase: int,
        stats: list[str],
        exact: bool,
        out: npt.NDArray[np.float32],
    ) -> None:
        """Writes ``stats`` (``mean min max sum std coverage``) over the bins
        of ``key`` into ``out``, C-contiguous ``(len(stats), bins)``: bin
        ``(x + phase) // resolution`` holds local position ``x``. ``exact``
        reads the full-resolution data, otherwise a zoom level."""
    def zoom_reductions(self) -> list[int]:
        """Bases per record of each zoom level, in file order."""
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
    """A bigBed file, read by region."""

    auto_sql: str | None
    """The autoSql schema, when the file has one."""
    field_count: tuple[int, int]
    """(fields, standard BED fields)."""

    def __new__(cls, path: str, cache: BlockCache | None = None) -> IntervalReader: ...
    def set_cache(self, cache: BlockCache | None) -> None: ...
    def cache_namespace(self) -> CacheNamespace: ...
    def describe_cache(self, chrom: str, start: int, end: int) -> list[BlockKey]: ...
    def prefetch_cache(self, chrom: str, start: int, end: int) -> None: ...
    def clone_context(self) -> IntervalReader: ...
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length) of every chromosome."""
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

    def __new__(cls, path: str, cache: BlockCache | None = None) -> ContactReader: ...
    def set_cache(self, cache: BlockCache | None) -> None: ...
    def cache_namespace(self) -> CacheNamespace: ...
    def describe_cache(self, tile: Tile) -> list[BlockKey]: ...
    def prefetch_cache(self, tile: Tile) -> None: ...
    def clone_context(self) -> ContactReader: ...
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
    def matrix(
        self,
        tile: Tile,
        out: npt.NDArray[np.float32],
        flip_rows: bool = False,
        flip_cols: bool = False,
    ) -> None:
        """Writes the contacts of ``tile`` row-major into ``out``: a
        C-contiguous ``float32`` array of (row bins, column bins), 0 where
        there is no contact; rows (columns) run from the high coordinate
        down with ``flip_rows`` (``flip_cols``). ``out`` is left intact on
        error."""
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

# Consensus: BCF + CSI over a reference

@final
class ConsensusBuilder:
    """Personal sequences from a BCF read by region through its CSI
    (``index``, or the one beside the file, which must exist) over
    ``reference``, as ``bcftools consensus`` 1.23 builds them. The options
    are ``bcftools consensus``'s: ``-I``, ``-M``, ``-a``, ``--mark-del``,
    ``--mark-ins``, ``--mark-snv``, ``-m``, ``--mask-with``,
    ``--regions-overlap``, ``-c``, ``-i`` and ``-e`` (``include`` and
    ``exclude`` are filter expressions; at most one is given)."""

    def __new__(
        cls,
        bcf: str,
        index: str | None,
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
    def clone_context(self, reference: SequenceReader) -> ConsensusBuilder: ...
    def samples(self) -> list[str]:
        """The samples of the BCF, in header order."""
    @property
    def index_path(self) -> str: ...
    def sequence(
        self, chrom: str, start: int, end: int, strand: str, selector: Selector
    ) -> tuple[U8, str | None]:
        """The consensus of ``chrom:[start, end)`` with ``selector``, 5' to
        3' of ``strand`` (``+ - . ?``; reverse complement on ``-``): its
        ASCII bases, and its chain (forward strand) with ``chain``."""
    def sequences(
        self, chrom: str, start: int, end: int, strand: str, selectors: list[Selector]
    ) -> tuple[I64, U8, list[str | None]]:
        """The consensus of ``chrom:[start, end)`` with each selector, 5' to
        3' of ``strand``, packed: sequence ``i`` is
        ``bases[offsets[i]:offsets[i + 1]]``; returns (offsets, bases, chains)."""
    def close(self) -> None:
        """Releases the BCF; the reference stays open."""

# Alignments: BAM

@final
class AlignmentReader:
    """A BAM file with its BAI or CSI index, or a CRAM file with its
    ``.crai`` (``index``, or the one htslib finds beside the file, which
    must exist) and the reference its slices need (a ``SequenceReader`` or
    the path of a FASTA / bgzip FASTA / 2bit file; CRAM only)."""

    def __new__(
        cls,
        path: str,
        index: str | None = None,
        cache: BlockCache | None = None,
        reference: SequenceReader | str | None = None,
    ) -> AlignmentReader: ...
    def set_cache(self, cache: BlockCache | None) -> None: ...
    def cache_namespace(self) -> CacheNamespace: ...
    def describe_cache(self, chrom: str, start: int, end: int) -> list[BlockKey]: ...
    def prefetch_cache(self, chrom: str, start: int, end: int) -> None: ...
    def clone_context(self) -> AlignmentReader: ...
    def stream_stats(self) -> dict[str, int]: ...
    def reset_stream_stats(self) -> None: ...
    def chroms(self) -> list[tuple[str, int]]:
        """(name, length) of every reference, in reference-id order."""
    def header_text(self) -> bytes:
        """The SAM header text as stored."""
    @property
    def index_path(self) -> str:
        """The index file read."""
    def query(
        self, chrom: str, start: int, end: int, fields: int
    ) -> tuple[
        I64, I64, U16, U8, I32, I32, I32, Strings, Strings, Strings, tuple[I64, U32] | None, Strings
    ]:
        """The records overlapping ``chrom:[start, end)``: (start, end, flag,
        mapq, mate_tid, mate_start, tlen), then names, sequences, qualities,
        CIGARs and aux fields, each ``None`` unless the ``fields`` bits ask
        for it."""
    def close(self) -> None:
        """Releases the file; further calls raise ``ValueError``."""

# Variants: BCF

@final
class VariantReader:
    """A BCF file with its CSI index (``index``, or ``<bcf>.csi``), read by
    region."""

    def __new__(
        cls, path: str, index: str | None = None, cache: BlockCache | None = None
    ) -> VariantReader: ...
    def set_cache(self, cache: BlockCache | None) -> None: ...
    def cache_namespace(self) -> CacheNamespace: ...
    def describe_cache(self, chrom: str, start: int, end: int) -> list[BlockKey]: ...
    def prefetch_cache(self, chrom: str, start: int, end: int) -> None: ...
    def clone_context(self) -> VariantReader: ...
    def stream_stats(self) -> dict[str, int]: ...
    def reset_stream_stats(self) -> None: ...
    def chroms(self) -> list[tuple[str, int | None]]:
        """(name, length or None) of every contig the header declares."""
    def samples(self) -> list[str]:
        """Sample names in header order."""
    @property
    def index_path(self) -> str:
        """The index file read."""
    def genotypes(
        self, chrom: str, start: int, end: int, samples: list[int] | None = None
    ) -> tuple[I64, I64, tuple[I64, U8], tuple[I64, U8], tuple[I64, U8], int, I8, BOOL]:
        """The records overlapping ``chrom:[start, end)``: (start, end, ids,
        REF, ALT, ploidy, gt, phased); ``gt`` is ``[records, samples,
        ploidy]`` and ``phased`` ``[records, samples]``, flattened."""
    def close(self) -> None:
        """Releases the file; further calls raise ``ValueError``."""

# Annotation: GTF / GFF3

Columns = dict[str, npt.NDArray[np.generic] | tuple[I64, U8] | None]
"""Named numpy columns; string columns are (offsets, data)."""

@final
class IntervalSet:
    """Stranded intervals with an owner row and an optional name, indexed."""

    @staticmethod
    def from_columns(
        chroms: list[str],
        seqid: U32,
        start: U32,
        end: U32,
        strand: I8,
        owner: U32 | None = None,
        names: list[str] | None = None,
        owner_level: Literal["gene", "transcript", "exon"] | None = None,
    ) -> IntervalSet: ...
    def __len__(self) -> int: ...
    def chroms(self) -> list[str]: ...
    @property
    def owner_level(self) -> Literal["gene", "transcript", "exon"] | None: ...
    def columns(self) -> Columns:
        """seqid, start, end, strand (int8), owner, name (or ``None``)."""
    def overlaps(self, chrom: str, start: int, end: int) -> U32: ...
    def overlaps_batch(
        self, chroms: list[str], starts: U32, ends: U32
    ) -> tuple[npt.NDArray[np.uint64], U32]: ...

@final
class Annotation:
    """A GTF or GFF3 file read whole (plain, gzip or BGZF). ``attributes``
    is ``"all"``, ``"none"`` or ``"only"`` (with ``keys``)."""

    def __new__(
        cls,
        path: str,
        attributes: Literal["all", "none", "only"] | None = None,
        keys: list[str] | None = None,
        aliases: list[tuple[str, str]] | None = None,
        strict: bool = False,
        threads: int = 0,
        dialect: Literal["gtf", "gff3"] | None = None,
    ) -> Annotation: ...
    @property
    def dialect(self) -> Literal["gtf", "gff3"]: ...
    def chroms(self) -> list[str]: ...
    def heap_bytes(self) -> int: ...
    def genes(self) -> Columns: ...
    def transcripts(self) -> Columns: ...
    def features(self, exons_only: bool = False) -> Columns: ...
    def exons(self) -> U32:
        """Feature rows of exons."""
    def report(self) -> dict[str, object]: ...
    def gene_rows(self, key: str, by: Literal["id", "name"]) -> U32: ...
    def transcript_rows(self, key: str, by: Literal["id", "name"]) -> U32: ...
    def transcript_parts(
        self, rows: U32, part: Literal["exons", "cds", "utr5", "utr3"]
    ) -> tuple[U32, I8, I64, U32, U32]:
        """``part`` of each transcript row as (seqid, strand sign, CSR
        offsets, span starts, span ends)."""
    def attributes(
        self, level: Literal["gene", "transcript", "feature"], row: int
    ) -> list[tuple[str, str]]: ...
    def overlaps(self, level: str, chrom: str, start: int, end: int) -> U32:
        """Rows at ``level`` overlapping the region; exon rows are positions
        in ``exons()``."""
    def overlaps_batch(
        self, level: str, chroms: list[str], starts: U32, ends: U32
    ) -> tuple[npt.NDArray[np.uint64], U32]: ...
    def nearest(self, level: str, chrom: str, start: int, end: int) -> tuple[U32, int]: ...
    def select(
        self,
        gene_biotypes: list[str] | None = None,
        transcript_biotypes: list[str] | None = None,
        tags: list[str] | None = None,
        tag_values: list[str] | None = None,
        genes: list[str] | None = None,
        chroms: list[str] | None = None,
    ) -> tuple[U32, U32]:
        """(gene rows, transcript rows)."""
    def introns(self, genes: U32 | None = None, transcripts: U32 | None = None) -> IntervalSet: ...
    def exonic_regions(
        self, genes: U32 | None = None, transcripts: U32 | None = None
    ) -> IntervalSet: ...
    def junctions(
        self, genes: U32 | None = None, transcripts: U32 | None = None
    ) -> tuple[IntervalSet, U32, U32, U32, U32, U32, U32]:
        """(set, gene offsets, genes, transcript offsets, transcripts, donor,
        acceptor)."""
    def exon_bins(
        self, aggregate: bool = False, genes: U32 | None = None, transcripts: U32 | None = None
    ) -> tuple[IntervalSet, U32, U32, U32, U32, U32]:
        """(set, part, transcript offsets, transcripts, gene offsets, genes)."""
    def windows(
        self,
        kind: Literal["promoter", "tes", "body"],
        level: Literal["gene", "transcript"],
        upstream: int,
        downstream: int,
        chrom_lengths: list[tuple[str, int]] | None = None,
        genes: U32 | None = None,
        transcripts: U32 | None = None,
    ) -> IntervalSet: ...
    def match_chroms(
        self, names: list[str], aliases: list[tuple[str, str]] | None = None
    ) -> I64: ...
    def write_gtf(self, path: str, layout: Literal["auto", "original", "model"]) -> None: ...
    def write_bed12(self, path: str, transcripts: U32 | None = None) -> None: ...

# Quantification (omicsweave.quant). Products are dicts of primitives with a
# "kind" of "coverage", "counts", "atac", "junctions", "labeling",
# "methylation", "hic" or "alleles".

def plan_quant(
    inputs: list[str],
    out_dir: str,
    products: list[dict[str, object]],
    *,
    annotation: Annotation | None = None,
    aliases: list[tuple[str, str]] | None = None,
    matrix: Literal["cpk", "mtx", "csv"] = "cpk",
    reference: str | None = None,
    variants: str | None = None,
    threads: int | None = None,
    memory: str = "16G",
) -> tuple[int, int, int, int, int, list[str], list[tuple[str, int]]]:
    """Checks a quantification against its budget, reading only the BAM
    headers; returns (threads, memory budget, fixed, pipeline, buffer bytes,
    needs, [(product stage, buffer bytes)])."""

def quantify(
    inputs: list[str],
    out_dir: str,
    products: list[dict[str, object]],
    *,
    annotation: Annotation | None = None,
    aliases: list[tuple[str, str]] | None = None,
    matrix: Literal["cpk", "mtx", "csv"] = "cpk",
    reference: str | None = None,
    variants: str | None = None,
    threads: int | None = None,
    memory: str = "16G",
    tmp_dir: str | None = None,
    progress: Progress | None = None,
) -> tuple[list[tuple[str, str]], int, int, float, str]:
    """Runs a quantification; returns ((kind, path) of every file written,
    records, bytes read, seconds, summary JSON)."""
