"""Gene annotation (GTF 2.2 / GFF3) as columnar tables.

:meth:`Annotation.read` reads GENCODE, Ensembl, RefSeq and UCSC GTF and
GFF3 files (plain, gzip or BGZF) into genes, transcripts and features with
0-based, half-open ``uint32`` coordinates. Columns are numpy arrays and
strings are :class:`~omicsweave.types.StringVector`; rows refer to each
other by index (a transcript's ``gene``, a gene's ``transcript_start`` /
``transcript_end`` range). Overlap queries return row indices; derived views
(introns, exonic regions, junctions, exon bins, promoter / TES / gene-body
windows) are :class:`IntervalSet` objects whose ``owner`` column names the
gene or transcript of each row.

Strands are ``int8``: ``1`` (``+``), ``-1`` (``-``), ``0`` (``.`` or
``?``).
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Literal

import numpy as np

from . import _native
from .formats import PathLike
from .types import Region, RegionLike, Spliced, StringVector
from .types._genome import as_region

__all__ = [
    "FEATURE_KINDS",
    "TAGS",
    "Annotation",
    "AnnotationReport",
    "ExonBins",
    "FeatureTable",
    "Gene",
    "GeneTable",
    "IntervalSet",
    "Junctions",
    "Level",
    "Part",
    "Problem",
    "Selection",
    "Transcript",
    "TranscriptTable",
]

Level = Literal["gene", "transcript", "exon"]

Part = Literal["exons", "cds", "utr5", "utr3"]
"""A part of a transcript read as a spliced key: ``exons`` (the mature
transcript), ``cds`` (the coding blocks, stop codon included, so a complete
ORF is a multiple of three long), ``utr5`` and ``utr3`` (the exonic bases
5' and 3' of them)."""

_STRAND = {1: "+", -1: "-", 0: "."}

FEATURE_KINDS = (
    "exon",
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "UTR",
    "start_codon",
    "stop_codon",
    "Selenocysteine",
    "other",
)
"""Names of the :attr:`FeatureTable.kind` codes."""

TAGS = (
    "basic",
    "Ensembl_canonical",
    "MANE_Select",
    "MANE_Plus_Clinical",
    "CCDS",
    "GENCODE_Primary",
    "appris_principal",
    "appris_alternative",
    "RefSeq_Select",
    "mRNA_start_NF",
    "mRNA_end_NF",
    "cds_start_NF",
    "cds_end_NF",
    "readthrough_transcript",
)
"""Known transcript tags: bit ``i`` of :attr:`TranscriptTable.tag_bits` is
``TAGS[i]`` (spelling and case variants, such as RefSeq's ``MANE Select``,
map to the same bit)."""


def _strings(column: tuple[np.ndarray, np.ndarray]) -> StringVector:
    return StringVector(*column)


@dataclass(frozen=True, slots=True, eq=False)
class GeneTable:
    """Genes sorted by (chromosome, start, end); gene ``g``'s transcripts are
    rows ``transcript_start[g]:transcript_end[g]`` of the transcript table.
    ``synthesized`` genes had no record and span their transcripts."""

    chroms: tuple[str, ...]
    seqid: np.ndarray
    start: np.ndarray
    end: np.ndarray
    strand: np.ndarray
    id: StringVector
    name: StringVector
    biotype: StringVector
    source: StringVector
    type: StringVector
    transcript_start: np.ndarray
    transcript_end: np.ndarray
    synthesized: np.ndarray

    def __len__(self) -> int:
        return len(self.start)


@dataclass(frozen=True, slots=True, eq=False)
class TranscriptTable:
    """Transcripts, contiguous per gene and sorted by position within it.
    Transcript ``t``'s features are rows ``feature_start[t]:feature_end[t]``
    of the feature table; its tags are
    ``tag_values[tag_offsets[t]:tag_offsets[t + 1]]``, and the known ones
    are also bits of ``tag_bits`` (see :data:`TAGS`). ``length`` is the
    spliced length, ``cds_length`` the summed CDS."""

    chroms: tuple[str, ...]
    seqid: np.ndarray
    start: np.ndarray
    end: np.ndarray
    strand: np.ndarray
    gene: np.ndarray
    id: StringVector
    name: StringVector
    biotype: StringVector
    source: StringVector
    type: StringVector
    tag_bits: np.ndarray
    tag_values: StringVector
    tag_offsets: np.ndarray
    feature_start: np.ndarray
    feature_end: np.ndarray
    length: np.ndarray
    cds_length: np.ndarray
    synthesized: np.ndarray

    def __len__(self) -> int:
        return len(self.start)

    def has_tag(self, tag: str) -> np.ndarray:
        """Which transcripts carry known tag ``tag`` (a name of :data:`TAGS`)."""
        if tag not in TAGS:
            raise ValueError(f"unknown tag {tag!r}; known: {', '.join(TAGS)}")
        return (self.tag_bits & np.uint32(1 << TAGS.index(tag))) != 0

    def tags(self, t: int) -> list[str]:
        """Every tag value of transcript ``t``."""
        return self.tag_values[int(self.tag_offsets[t]) : int(self.tag_offsets[t + 1])]


@dataclass(frozen=True, slots=True, eq=False)
class FeatureTable:
    """Sub-transcript features (exons, CDS, UTRs, codons, ...), contiguous
    per transcript and sorted by (start, end, kind). ``kind`` codes name
    :data:`FEATURE_KINDS`; ``phase`` is 0-2, or 255 for ``.``; ``feature``
    is the row in the full feature table (for :attr:`Annotation.exons`)."""

    chroms: tuple[str, ...]
    feature: np.ndarray
    transcript: np.ndarray
    seqid: np.ndarray
    kind: np.ndarray
    start: np.ndarray
    end: np.ndarray
    strand: np.ndarray
    phase: np.ndarray
    type: StringVector

    def __len__(self) -> int:
        return len(self.start)


@dataclass(frozen=True, slots=True)
class Problem:
    """A problem found while reading: its kind (``orphan``-style names such
    as ``unknown_parent``, ``strand_conflict``, ``outside_parent``), the
    1-based line, and the id involved (``""`` when none)."""

    kind: str
    line: int
    id: str


@dataclass(frozen=True, slots=True)
class AnnotationReport:
    """What reading found: counts, and every problem in line order."""

    records: int
    comments: int
    genes: int
    transcripts: int
    features: int
    others: int
    synthesized_genes: int
    synthesized_transcripts: int
    problems: tuple[Problem, ...]

    @property
    def is_clean(self) -> bool:
        return not self.problems

    def count(self, kind: str) -> int:
        return sum(p.kind == kind for p in self.problems)


@dataclass(frozen=True, slots=True)
class Selection:
    """Gene and transcript rows, ascending; every selected transcript's gene
    is selected."""

    genes: np.ndarray
    transcripts: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class Gene:
    """One gene: its row, fields, transcript rows and kept attributes."""

    index: int
    id: str
    name: str
    biotype: str
    region: Region
    strand: int
    transcripts: np.ndarray
    attributes: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True, eq=False)
class Transcript:
    """One transcript: its row, gene row, fields, feature rows and kept
    attributes."""

    index: int
    gene: int
    id: str
    name: str
    biotype: str
    region: Region
    strand: int
    tags: tuple[str, ...]
    features: np.ndarray
    attributes: tuple[tuple[str, str], ...]


class IntervalSet:
    """Stranded intervals with an owner (a gene or transcript row, per
    :attr:`owner_level`) and an optional name per row, indexed for overlap
    queries. Derived views of an annotation are interval sets; so can be
    peaks or bins (:meth:`from_arrays`)."""

    def __init__(self, native: Any) -> None:
        self._native = native

    @classmethod
    def from_arrays(
        cls,
        chroms: Sequence[str],
        seqid: np.ndarray,
        start: np.ndarray,
        end: np.ndarray,
        strand: np.ndarray | None = None,
        owner: np.ndarray | None = None,
        names: Sequence[str] | None = None,
    ) -> IntervalSet:
        """A set from columns; ``seqid`` numbers ``chroms``."""
        start = np.ascontiguousarray(start, dtype=np.uint32)
        if strand is None:
            strand = np.zeros(len(start), dtype=np.int8)
        native = _native.IntervalSet.from_columns(
            list(chroms),
            np.ascontiguousarray(seqid, dtype=np.uint32),
            start,
            np.ascontiguousarray(end, dtype=np.uint32),
            np.ascontiguousarray(strand, dtype=np.int8),
            None if owner is None else np.ascontiguousarray(owner, dtype=np.uint32),
            None if names is None else list(names),
        )
        return cls(native)

    def __len__(self) -> int:
        return len(self._native)

    @cached_property
    def chroms(self) -> tuple[str, ...]:
        return tuple(self._native.chroms())

    @cached_property
    def _columns(self) -> dict[str, Any]:
        return self._native.columns()

    @property
    def seqid(self) -> np.ndarray:
        return self._columns["seqid"]

    @property
    def start(self) -> np.ndarray:
        return self._columns["start"]

    @property
    def end(self) -> np.ndarray:
        return self._columns["end"]

    @property
    def strand(self) -> np.ndarray:
        return self._columns["strand"]

    @property
    def owner(self) -> np.ndarray:
        """Owner rows (``2**32 - 1`` for none)."""
        return self._columns["owner"]

    @property
    def owner_level(self) -> Level | None:
        return self._native.owner_level

    @cached_property
    def name(self) -> StringVector | None:
        names = self._columns["name"]
        return None if names is None else _strings(names)

    def overlaps(self, region: RegionLike) -> np.ndarray:
        """Rows overlapping ``region``, in start order."""
        r = as_region(region)
        return self._native.overlaps(r.chrom, r.start, r.end)

    def overlaps_batch(self, regions: Iterable[RegionLike]) -> tuple[np.ndarray, np.ndarray]:
        """Rows overlapping each region: ``(offsets, rows)``, region ``i``'s
        rows being ``rows[offsets[i]:offsets[i + 1]]``."""
        return _batch(self._native.overlaps_batch, regions)

    def __repr__(self) -> str:
        return f"IntervalSet(n={len(self)}, owner_level={self.owner_level!r})"


def _batch(call: Any, regions: Iterable[RegionLike], *lead: Any) -> tuple[np.ndarray, np.ndarray]:
    rs = [as_region(r) for r in regions]
    starts = np.fromiter((r.start for r in rs), dtype=np.uint32, count=len(rs))
    ends = np.fromiter((r.end for r in rs), dtype=np.uint32, count=len(rs))
    return call(*lead, [r.chrom for r in rs], starts, ends)


@dataclass(frozen=True, slots=True, eq=False)
class Junctions:
    """Annotated splice junctions (distinct introns): ``intervals`` sorted by
    (chromosome, start, end, strand), owner = first gene; ``donor`` and
    ``acceptor`` are the first and last intron bases on the transcript
    strand; junction ``i``'s genes are
    ``genes[gene_offsets[i]:gene_offsets[i + 1]]``, likewise transcripts."""

    intervals: IntervalSet
    donor: np.ndarray
    acceptor: np.ndarray
    gene_offsets: np.ndarray
    genes: np.ndarray
    transcript_offsets: np.ndarray
    transcripts: np.ndarray

    def __len__(self) -> int:
        return len(self.intervals)


@dataclass(frozen=True, slots=True, eq=False)
class ExonBins:
    """DEXSeq-like disjoint exon bins: the exonic bases of a gene (or of an
    aggregate of genes whose exons overlap on one strand) cut where the set
    of covering transcripts changes. ``part`` numbers each gene's bins from
    1; bin ``i``'s transcripts are
    ``transcripts[transcript_offsets[i]:transcript_offsets[i + 1]]``, its
    aggregate's genes likewise."""

    intervals: IntervalSet
    part: np.ndarray
    transcript_offsets: np.ndarray
    transcripts: np.ndarray
    gene_offsets: np.ndarray
    genes: np.ndarray

    def __len__(self) -> int:
        return len(self.intervals)


def _reopen(path: str, kwargs: dict[str, Any]) -> Annotation:
    return Annotation.read(path, **kwargs)


class Annotation:
    """A gene annotation read whole; immutable, shared safely by threads,
    pickled as its path and read options.

    Native queries are process-local. After fork, reopen the annotation
    (or unpickle it) in the child; existing exported numpy columns can
    still be read as ordinary arrays.
    """

    def __init__(self, native: Any, path: str, kwargs: dict[str, Any]) -> None:
        self._native = native
        self._path = path
        self._kwargs = kwargs

    @classmethod
    def read(
        cls,
        path: PathLike,
        attributes: Literal["all"] | Sequence[str] | None = "all",
        aliases: Mapping[str, str] | None = None,
        strict: bool = False,
        threads: int = 0,
        dialect: Literal["gtf", "gff3"] | None = None,
    ) -> Annotation:
        """Reads a GTF or GFF3 file (plain, gzip or BGZF, told by content).

        ``attributes`` keeps every attribute (``"all"``: lossless, needed to
        write the original GTF back), only the keys listed, or none
        (``None``); ids, names, biotypes and tags are read regardless.
        ``aliases`` renames chromosomes (``{"1": "chr1"}``). Problems
        (orphans, duplicate ids, strand conflicts, ...) are collected in
        :attr:`report`, or raised as
        :class:`~omicsweave.errors.FormatError` when ``strict``.
        """
        if attributes == "all":
            mode, keys = "all", None
        elif attributes is None:
            mode, keys = "none", None
        elif isinstance(attributes, str):
            raise ValueError('attributes must be "all", None or a list of keys')
        else:
            mode, keys = "only", [str(k) for k in attributes]
        pairs = None if aliases is None else [(str(a), str(b)) for a, b in aliases.items()]
        fspath = os.fspath(path)
        native = _native.Annotation(fspath, mode, keys, pairs, strict, threads, dialect)
        kwargs = dict(
            attributes=attributes
            if isinstance(attributes, str) or attributes is None
            else list(attributes),
            aliases=None if aliases is None else dict(aliases),
            strict=strict,
            threads=threads,
            dialect=dialect,
        )
        return cls(native, fspath, kwargs)

    def __reduce__(self) -> tuple[Any, ...]:
        return (_reopen, (self._path, self._kwargs))

    def __repr__(self) -> str:
        return (
            f"Annotation({self._path!r}, dialect={self.dialect!r}, genes={len(self.genes)}, "
            f"transcripts={len(self.transcripts)})"
        )

    @property
    def path(self) -> str:
        return self._path

    @property
    def dialect(self) -> Literal["gtf", "gff3"]:
        return self._native.dialect

    @cached_property
    def chroms(self) -> tuple[str, ...]:
        """Chromosome names in order of first appearance; ``seqid`` columns
        index this."""
        return tuple(self._native.chroms())

    @property
    def heap_bytes(self) -> int:
        """Memory held by the native tables, strings and built caches."""
        return self._native.heap_bytes()

    @cached_property
    def genes(self) -> GeneTable:
        c = self._native.genes()
        return GeneTable(
            self.chroms,
            c["seqid"],
            c["start"],
            c["end"],
            c["strand"],
            _strings(c["id"]),
            _strings(c["name"]),
            _strings(c["biotype"]),
            _strings(c["source"]),
            _strings(c["type"]),
            c["transcript_start"],
            c["transcript_end"],
            c["synthesized"],
        )

    @cached_property
    def transcripts(self) -> TranscriptTable:
        c = self._native.transcripts()
        return TranscriptTable(
            self.chroms,
            c["seqid"],
            c["start"],
            c["end"],
            c["strand"],
            c["gene"],
            _strings(c["id"]),
            _strings(c["name"]),
            _strings(c["biotype"]),
            _strings(c["source"]),
            _strings(c["type"]),
            c["tag_bits"],
            _strings(c["tag_values"]),
            c["tag_offsets"],
            c["feature_start"],
            c["feature_end"],
            c["length"],
            c["cds_length"],
            c["synthesized"],
        )

    def _features(self, exons_only: bool) -> FeatureTable:
        c = self._native.features(exons_only)
        return FeatureTable(
            self.chroms,
            c["feature"],
            c["transcript"],
            c["seqid"],
            c["kind"],
            c["start"],
            c["end"],
            c["strand"],
            c["phase"],
            _strings(c["type"]),
        )

    @cached_property
    def features(self) -> FeatureTable:
        return self._features(False)

    @cached_property
    def exons(self) -> FeatureTable:
        """The exon rows of :attr:`features`; exon-level queries return rows
        of this table."""
        return self._features(True)

    @cached_property
    def report(self) -> AnnotationReport:
        r = self._native.report()
        problems = tuple(
            Problem(k, int(line), i)
            for k, line, i in zip(
                r["problem_kinds"], r["problem_lines"], r["problem_ids"], strict=True
            )
        )
        return AnnotationReport(
            r["records"],
            r["comments"],
            r["genes"],
            r["transcripts"],
            r["features"],
            r["others"],
            r["synthesized_genes"],
            r["synthesized_transcripts"],
            problems,
        )

    def _attrs(self, level: str, row: int) -> tuple[tuple[str, str], ...]:
        return tuple(self._native.attributes(level, row))

    def gene_rows(self, key: str, by: Literal["id", "name"] = "id") -> np.ndarray:
        """Rows of the genes with this id or name (an id may repeat on
        several chromosomes, a name across genes)."""
        return self._native.gene_rows(key, by)

    def transcript_rows(self, key: str, by: Literal["id", "name"] = "id") -> np.ndarray:
        return self._native.transcript_rows(key, by)

    def gene(self, key: str) -> Gene:
        """The first gene with id ``key``, else with name ``key``."""
        rows = self.gene_rows(key)
        if len(rows) == 0:
            rows = self.gene_rows(key, "name")
        if len(rows) == 0:
            raise KeyError(key)
        g = int(rows[0])
        t = self.genes
        return Gene(
            g,
            t.id[g],
            t.name[g],
            t.biotype[g],
            Region(self.chroms[int(t.seqid[g])], int(t.start[g]), int(t.end[g])),
            int(t.strand[g]),
            np.arange(t.transcript_start[g], t.transcript_end[g], dtype=np.uint32),
            self._attrs("gene", g),
        )

    def transcript(self, key: str) -> Transcript:
        """The first transcript with id ``key``, else with name ``key``."""
        rows = self.transcript_rows(key)
        if len(rows) == 0:
            rows = self.transcript_rows(key, "name")
        if len(rows) == 0:
            raise KeyError(key)
        i = int(rows[0])
        t = self.transcripts
        return Transcript(
            i,
            int(t.gene[i]),
            t.id[i],
            t.name[i],
            t.biotype[i],
            Region(self.chroms[int(t.seqid[i])], int(t.start[i]), int(t.end[i])),
            int(t.strand[i]),
            tuple(t.tags(i)),
            np.arange(t.feature_start[i], t.feature_end[i], dtype=np.uint32),
            self._attrs("transcript", i),
        )

    def locus(self, transcript: str | int, part: Part = "exons") -> Spliced:
        """``part`` of a transcript (by id, else name, or by row) as a
        :class:`~omicsweave.types.Spliced` key on its strand, named by its
        id: reading it gives the spliced bases or signal 5' to 3'."""
        if isinstance(transcript, (int, np.integer)):
            row = int(transcript)
        else:
            row = self.transcript(transcript).index
        return self.loci([row], part)[0]

    def loci(
        self, transcripts: Iterable[int] | np.ndarray | None = None, part: Part = "exons"
    ) -> list[Spliced]:
        """:meth:`locus` of every transcript row given (default: all), in
        one native call."""
        if transcripts is None:
            rows = np.arange(len(self.transcripts), dtype=np.uint32)
        else:
            rows = np.ascontiguousarray(np.fromiter(transcripts, dtype=np.uint32))
        seqid, strand, offsets, starts, ends = self._native.transcript_parts(rows, part)
        chroms = self.chroms
        ids = self.transcripts.id
        return [
            Spliced(
                chroms[int(seqid[k])],
                zip(
                    starts[offsets[k] : offsets[k + 1]].tolist(),
                    ends[offsets[k] : offsets[k + 1]].tolist(),
                    strict=True,
                ),
                _STRAND[int(strand[k])],
                ids[int(r)],
            )
            for k, r in enumerate(rows)
        ]

    def feature_attributes(self, row: int) -> tuple[tuple[str, str], ...]:
        """The kept attributes of feature row ``row``."""
        return self._attrs("feature", row)

    def overlaps(self, region: RegionLike, level: Level = "gene") -> np.ndarray:
        """Rows of the gene, transcript or exon table overlapping ``region``,
        in start order; an unknown chromosome has none."""
        r = as_region(region)
        return self._native.overlaps(level, r.chrom, r.start, r.end)

    def overlaps_batch(
        self, regions: Iterable[RegionLike], level: Level = "gene"
    ) -> tuple[np.ndarray, np.ndarray]:
        """Rows overlapping each region: ``(offsets, rows)``."""
        return _batch(self._native.overlaps_batch, regions, level)

    def nearest(self, region: RegionLike, level: Level = "gene") -> tuple[np.ndarray, int]:
        """The rows nearest to ``region`` (all ties) and their distance in
        bases: 0 when overlapping, ``-1`` when the chromosome has none."""
        r = as_region(region)
        return self._native.nearest(level, r.chrom, r.start, r.end)

    def subset(
        self,
        biotypes: Iterable[str] | None = None,
        transcript_biotypes: Iterable[str] | None = None,
        tags: Iterable[str] | None = None,
        tag_values: Iterable[str] | None = None,
        genes: Iterable[str] | None = None,
        chroms: Iterable[str] | None = None,
    ) -> Selection:
        """Genes and transcripts passing every filter given: gene biotypes,
        transcript biotypes, all of the known ``tags`` (:data:`TAGS`), any of
        the ``tag_values``, gene ids or names, chromosomes. A gene passes
        transcript filters when one of its transcripts does."""

        def lst(v: Iterable[str] | None) -> list[str] | None:
            return None if v is None else [str(x) for x in v]

        g, t = self._native.select(
            lst(biotypes),
            lst(transcript_biotypes),
            lst(tags),
            lst(tag_values),
            lst(genes),
            lst(chroms),
        )
        return Selection(g, t)

    @staticmethod
    def _scope(subset: Selection | None) -> tuple[Any, Any]:
        if subset is None:
            return None, None
        return (
            np.ascontiguousarray(subset.genes, dtype=np.uint32),
            np.ascontiguousarray(subset.transcripts, dtype=np.uint32),
        )

    def introns(self, subset: Selection | None = None) -> IntervalSet:
        """Introns of each transcript (owner: transcript row)."""
        return IntervalSet(self._native.introns(*self._scope(subset)))

    def exonic_regions(self, subset: Selection | None = None) -> IntervalSet:
        """The union of each gene's exons (owner: gene row)."""
        return IntervalSet(self._native.exonic_regions(*self._scope(subset)))

    def junctions(self, subset: Selection | None = None) -> Junctions:
        s, go, gv, to, tv, donor, acceptor = self._native.junctions(*self._scope(subset))
        return Junctions(IntervalSet(s), donor, acceptor, go, gv, to, tv)

    def exon_bins(self, aggregate: bool = False, subset: Selection | None = None) -> ExonBins:
        """Disjoint exon bins per gene, or per aggregate of genes whose exons
        overlap on one strand (DEXSeq's default)."""
        s, part, to, tv, go, gv = self._native.exon_bins(aggregate, *self._scope(subset))
        return ExonBins(IntervalSet(s), part, to, tv, go, gv)

    def _windows(
        self,
        kind: str,
        level: Level,
        upstream: int,
        downstream: int,
        chrom_lengths: Mapping[str, int] | None,
        subset: Selection | None,
    ) -> IntervalSet:
        lengths = (
            None if chrom_lengths is None else [(str(k), int(v)) for k, v in chrom_lengths.items()]
        )
        return IntervalSet(
            self._native.windows(kind, level, upstream, downstream, lengths, *self._scope(subset))
        )

    def promoters(
        self,
        upstream: int = 2000,
        downstream: int = 200,
        level: Level = "transcript",
        chrom_lengths: Mapping[str, int] | None = None,
        subset: Selection | None = None,
    ) -> IntervalSet:
        """``upstream`` bases before each TSS and ``downstream`` from it on,
        strand-aware (``GenomicRanges::promoters``), clipped at 0 and at the
        ``chrom_lengths`` given (matched as :meth:`match_chroms` does)."""
        return self._windows("promoter", level, upstream, downstream, chrom_lengths, subset)

    def tss(self, level: Level = "transcript", subset: Selection | None = None) -> IntervalSet:
        """The first transcribed base of each gene or transcript."""
        return self._windows("promoter", level, 0, 1, None, subset)

    def tes(self, level: Level = "transcript", subset: Selection | None = None) -> IntervalSet:
        """The last transcribed base of each gene or transcript."""
        return self._windows("tes", level, 1, 0, None, subset)

    def tes_windows(
        self,
        upstream: int = 0,
        downstream: int = 0,
        level: Level = "transcript",
        chrom_lengths: Mapping[str, int] | None = None,
        subset: Selection | None = None,
    ) -> IntervalSet:
        """``upstream`` bases inside each body before its TES and
        ``downstream`` past it."""
        return self._windows("tes", level, upstream, downstream, chrom_lengths, subset)

    def gene_bodies(
        self,
        upstream: int = 2000,
        downstream: int = 0,
        level: Level = "gene",
        chrom_lengths: Mapping[str, int] | None = None,
        subset: Selection | None = None,
    ) -> IntervalSet:
        """Each body extended ``upstream`` before its TSS and ``downstream``
        past its TES: the gene-activity regions of Signac (2 kb upstream,
        the default) or ArchR (``upstream=5000``)."""
        return self._windows("body", level, upstream, downstream, chrom_lengths, subset)

    def match_chroms(
        self, names: Sequence[str], aliases: Mapping[str, str] | None = None
    ) -> np.ndarray:
        """For each name (a BAM header's references, say), the index in
        :attr:`chroms` it refers to, ``-1`` for none: the same name, an
        alias, or the built-in ``chr1`` ↔ ``1`` and ``chrM`` ↔ ``MT``
        spellings."""
        pairs = None if aliases is None else [(str(a), str(b)) for a, b in aliases.items()]
        return self._native.match_chroms([str(n) for n in names], pairs)

    def write_gtf(
        self, path: PathLike, layout: Literal["auto", "original", "model"] = "auto"
    ) -> None:
        """Writes GTF: ``original`` reproduces a GTF read with all
        attributes (every record and comment, in file order); ``model``
        writes gene, transcript and feature lines from the tables; ``auto``
        picks ``original`` when possible."""
        self._native.write_gtf(os.fspath(path), layout)

    def write_bed12(self, path: PathLike, subset: Selection | None = None) -> None:
        """Writes transcripts (all, or the subset's) as BED12."""
        rows = None if subset is None else np.ascontiguousarray(subset.transcripts, dtype=np.uint32)
        self._native.write_bed12(os.fspath(path), rows)
