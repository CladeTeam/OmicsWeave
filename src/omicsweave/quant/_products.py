"""What one pass over the BAMs computes: product options as dataclasses."""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, Literal

from ..formats import PathLike

__all__ = [
    "AlleleGenes",
    "Alleles",
    "Atac",
    "Contacts",
    "Counts",
    "Coverage",
    "Filter",
    "Genes",
    "JunctionMatrix",
    "Junctions",
    "Labeled",
    "Labeling",
    "MethylMatrix",
    "Methylation",
    "Product",
]

UNMAPPED, PAIRED, PROPER_PAIR, MATE_UNMAPPED = 0x4, 0x1, 0x2, 0x8
SECONDARY, QC_FAIL, DUPLICATE, SUPPLEMENTARY = 0x100, 0x200, 0x400, 0x800

Library = Literal[0, 1, 2]
"""featureCounts' ``-s``: 0 unstranded, 1 read 1 on the transcript's strand,
2 read 1 on the opposite strand."""


@dataclass(frozen=True, slots=True)
class Filter:
    """Which records a product counts."""

    require: int = 0
    """FLAG bits a record must all have."""
    exclude: int = UNMAPPED | SECONDARY | QC_FAIL | DUPLICATE
    """FLAG bits a record must have none of (default ``0x704``, as
    ``samtools depth``)."""
    min_mapq: int = 0
    unique: bool = False
    """Only records whose ``NH`` tag is 1 (none counts as 1)."""

    def _native(self) -> tuple[int, int, int, bool]:
        return (self.require, self.exclude, self.min_mapq, self.unique)


def _path(p: PathLike | None) -> str | None:
    return None if p is None else os.fspath(p)


@dataclass(frozen=True, slots=True)
class Labeled:
    """A labeled layer of a gene matrix: units (reads, pairs, UMIs) with at
    least ``min_conversions`` conversions.

    ``conversion`` is in transcript sense (``"T>C"`` for 4sU: on the genome
    T>C for forward transcripts, A>G for reverse ones); bases below
    ``min_quality`` and positions of the task's ``variants`` are not counted.
    The reference is the task's ``reference``, else each read's ``MD`` tag.
    A pair's conversions are those of both mates; a molecule is labeled when
    any of its reads is. Under an unstranded library, the strand is the
    assigned gene's.
    """

    conversion: str = "T>C"
    min_quality: int = 20
    min_conversions: int = 1

    def _native(self) -> tuple[str, int, int]:
        return (self.conversion, self.min_quality, self.min_conversions)


@dataclass(frozen=True, slots=True)
class Coverage:
    """Coverage as bigWig: ``<name>.bw``, or ``<name>.forward.bw`` and
    ``<name>.reverse.bw`` when strands are split.

    ``unit``: ``"bases"`` (aligned ``M``/``=``/``X`` bases, as ``samtools
    depth``), ``"bases+deletions"`` (``samtools depth -J``), ``"span"``,
    ``"five-prime"``, ``"three-prime"``, or ``"fragment"`` (a pair once,
    from its leftmost mate, ``[pos, pos + |TLEN|)``, up to
    ``max_fragment`` bases when not 0). ``scale``: ``"none"``, ``"cpm"``
    (per million counted reads or fragments) or a factor.
    """

    name: str = "coverage"
    unit: Literal["bases", "bases+deletions", "span", "five-prime", "three-prime", "fragment"] = (
        "bases"
    )
    max_fragment: int = 0
    strands: Literal["merged", "split"] = "merged"
    library: Library = 0
    """With split strands: the strand of the transcript a read stands for
    (read 2 counts as read 1's strand; ``2`` flips)."""
    filter: Filter = field(default_factory=Filter)
    scale: Literal["none", "cpm"] | float = "none"
    sidecar: bool = False

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "coverage",
            "name": self.name,
            "unit": self.unit,
            "max_fragment": self.max_fragment,
            "strands": self.strands,
            "library": self.library,
            "filter": self.filter._native(),
            "scale": self.scale,
            "sidecar": self.sidecar,
        }


@dataclass(frozen=True, slots=True)
class Counts:
    """A count matrix, rows × features, written as a 10x directory
    ``<name>/`` (``barcodes.tsv.gz``, ``features.tsv.gz`` and the matrix).

    - ``features``: ``"genes"`` (merged exons), ``"gene-full"`` (gene spans),
      ``"promoters"`` / ``"gene-activity"`` (TSS windows / bodies extended
      by ``upstream`` and ``downstream``), ``"peaks"`` (every feature of the
      bigBed ``peaks``), ``"bins"`` (tiles of ``bin_size``).
    - ``unit``: ``"read"``, ``"molecule"`` (reads deduplicated by UMI),
      ``"pair"`` (both mates together, featureCounts ``-p
      --countReadPairs``), ``"fragment"`` / ``"tn5-fragment"`` (a pair's span
      from its leftmost mate), ``"cut-site"`` (Tn5 insertions).
    - ``rows``: a tag (``"CB"``), ``"@RG"`` (read groups) or ``"@file"``
      (one row per input BAM); ``whitelist`` restricts and orders them.
    - ``dedup``: ``"none"``, ``"umi"`` or ``"position"`` (fragments of a
      row counted once).
    - A unit overlapping more than one feature is ambiguous unless
      ``multi_overlap``.
    - ``table``: also write a featureCounts-style ``counts.tsv``.
    - ``velocity``: add the ``spliced``, ``unspliced`` and ``ambiguous`` layers
      (genes; velocyto / STARsolo ``Velocyto`` rules), written as
      ``spliced.cpk`` etc. beside ``matrix.cpk``.
    - ``labeled``: add a ``labeled`` layer (:class:`Labeled`).

    UMI and molecule identity follows the encoding and fingerprint rules
    documented in :mod:`omicsweave.quant`.
    """

    name: str
    features: Literal["genes", "gene-full", "promoters", "gene-activity", "peaks", "bins"]
    unit: Literal["read", "molecule", "pair", "fragment", "tn5-fragment", "cut-site"] = "read"
    rows: str = "@file"
    whitelist: PathLike | None = None
    umi_tag: str = "UB"
    dedup: Literal["none", "umi", "position"] | None = None
    library: Library = 0
    multi_overlap: bool = False
    filter: Filter = field(default_factory=Filter)
    max_fragment: int = 0
    table: bool = False
    upstream: int = 0
    downstream: int = 0
    peaks: PathLike | None = None
    bin_size: int | None = None
    biotypes: Sequence[str] | None = None
    velocity: bool = False
    labeled: Labeled | None = None

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "counts",
            "name": self.name,
            "features": self.features,
            "unit": self.unit,
            "rows": self.rows,
            "whitelist": _path(self.whitelist),
            "umi_tag": self.umi_tag,
            "dedup": self.dedup,
            "library": self.library,
            "multi_overlap": self.multi_overlap,
            "filter": self.filter._native(),
            "max_fragment": self.max_fragment,
            "table": self.table,
            "upstream": self.upstream,
            "downstream": self.downstream,
            "peaks": _path(self.peaks),
            "bin_size": self.bin_size,
            "biotypes": None if self.biotypes is None else list(self.biotypes),
            "velocity": self.velocity,
            "labeled": None if self.labeled is None else self.labeled._native(),
        }


def Genes(
    mode: Literal["single-cell", "bulk"] = "single-cell",
    *,
    name: str = "genes",
    span: Literal["exons", "full"] = "exons",
    pairs: bool = False,
    library: Library | None = None,
    barcode_tag: str = "CB",
    whitelist: PathLike | None = None,
    biotypes: Sequence[str] | None = None,
    multi_overlap: bool = False,
    table: bool | None = None,
    velocity: bool = False,
    labeled: Labeled | None = None,
) -> Counts:
    """Gene counts.

    ``"single-cell"`` (STARsolo ``Gene`` / Cell Ranger): genes × UMIs ×
    ``barcode_tag``, the raw matrix over every barcode, unique alignments,
    forward-stranded unless ``library`` says otherwise. ``"bulk"``
    (featureCounts defaults): genes × reads (``pairs``: read pairs) × input
    file, unstranded unless ``library``, every mapped record with ``NH`` 1
    counted, with a ``counts.tsv`` table. ``span="full"`` counts introns
    too (STARsolo ``GeneFull``). ``velocity`` and ``labeled`` add layers (see
    :class:`Counts`).
    """
    features = "genes" if span == "exons" else "gene-full"
    if mode == "single-cell":
        return Counts(
            name=name,
            features=features,
            unit="molecule",
            rows=barcode_tag,
            whitelist=whitelist,
            library=1 if library is None else library,
            multi_overlap=multi_overlap,
            filter=Filter(unique=True),
            table=bool(table),
            biotypes=biotypes,
            velocity=velocity,
            labeled=labeled,
        )
    return Counts(
        name=name,
        features=features,
        unit="pair" if pairs else "read",
        rows="@file",
        library=0 if library is None else library,
        multi_overlap=multi_overlap,
        filter=Filter(exclude=UNMAPPED, unique=True),
        table=True if table is None else table,
        biotypes=biotypes,
        velocity=velocity,
        labeled=labeled,
    )


@dataclass(frozen=True, slots=True)
class Atac:
    """Single-cell ATAC under ``<name>/``: ``fragments.tsv.gz`` (+ ``.tbi``),
    ``barcode_metrics.csv``, and cut-site matrices over ``peaks`` (bigBed)
    and ``tiles``, and fragments over gene bodies with ``gene_activity``
    (upstream, downstream) flanks when an annotation is given.

    A fragment is a proper pair from its leftmost mate, Tn5-shifted (``pos +
    4`` to ``pos + |TLEN| - 5``); duplicates of a barcode collapse, the fifth
    column counting them.
    """

    name: str = "atac"
    barcode_tag: str = "CB"
    whitelist: PathLike | None = None
    filter: Filter = field(
        default_factory=lambda: Filter(
            require=PAIRED | PROPER_PAIR,
            exclude=UNMAPPED | MATE_UNMAPPED | SECONDARY | QC_FAIL | SUPPLEMENTARY,
            min_mapq=30,
        )
    )
    max_fragment: int = 0
    fragments: bool = True
    peaks: PathLike | None = None
    tiles: int | None = 5000
    gene_activity: tuple[int, int] | None = (2000, 0)

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "atac",
            "name": self.name,
            "barcode_tag": self.barcode_tag,
            "whitelist": _path(self.whitelist),
            "filter": self.filter._native(),
            "max_fragment": self.max_fragment,
            "fragments": self.fragments,
            "peaks": _path(self.peaks),
            "tiles": self.tiles,
            "gene_activity": self.gene_activity,
        }


@dataclass(frozen=True, slots=True)
class JunctionMatrix:
    """Junctions × rows, beside the junction table: one column per junction
    of the table. ``rows`` is a tag (``"CB"``), ``"@RG"`` or ``"@file"``;
    ``dedup`` ``"umi"`` counts UMIs (STARsolo ``SJ``), ``"none"`` reads or
    pairs; ``filter`` applies on top of the table's."""

    rows: str = "CB"
    whitelist: PathLike | None = None
    umi_tag: str = "UB"
    dedup: Literal["none", "umi"] = "umi"
    filter: Filter = field(default_factory=lambda: Filter(unique=True))

    @staticmethod
    def bulk() -> JunctionMatrix:
        """Reads (pairs) per input file, unique alignments."""
        return JunctionMatrix(
            rows="@file", dedup="none", filter=Filter(exclude=UNMAPPED, unique=True)
        )


@dataclass(frozen=True, slots=True)
class Junctions:
    """Splice junctions as ``<name>/junctions.tab``, the nine columns of
    STAR's ``SJ.out.tab``: chromosome, first and last intron base (1-based),
    strand (0, 1 ``+``, 2 ``-``), motif (0 non-canonical, 1 GT/AG, 2 CT/AC,
    3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT; from the task's ``reference``),
    annotated, unique reads, multi-mapping reads, largest overhang.

    A pair crossing a junction counts once (``pairs``); every alignment of a
    multi-mapper counts. With ``matrix``, ``<name>/`` is also a 10x directory
    of junctions × rows.
    """

    name: str = "junctions"
    filter: Filter = field(
        default_factory=lambda: Filter(exclude=UNMAPPED | QC_FAIL | SUPPLEMENTARY)
    )
    pairs: bool = True
    matrix: JunctionMatrix | None = None

    def _native(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "kind": "junctions",
            "name": self.name,
            "filter": self.filter._native(),
            "pairs": self.pairs,
        }
        if self.matrix is not None:
            m = self.matrix
            out.update(
                matrix_rows=m.rows,
                matrix_whitelist=_path(m.whitelist),
                matrix_umi_tag=m.umi_tag,
                matrix_dedup=m.dedup,
                matrix_filter=m.filter._native(),
            )
        return out


@dataclass(frozen=True, slots=True)
class Labeling:
    """Metabolic labeling statistics under ``<name>/``:
    ``conversions.tsv.gz`` (reads by coverage of the conversion's source base
    and number of conversions) and ``mismatches.tsv`` (per mate, sequencing
    cycle and transcript-sense reference base, the read bases), with overall
    substitution rates in ``summary.json``.

    A read's transcript strand follows ``library`` (read 2 as read 1's
    strand; ``2`` flips; ``0`` is read as ``1``). Reference bases come from
    the task's ``reference``, else the ``MD`` tag; bases below
    ``min_quality`` and positions of the task's ``variants`` are skipped.
    """

    name: str = "labeling"
    filter: Filter = field(
        default_factory=lambda: Filter(
            exclude=UNMAPPED | SECONDARY | QC_FAIL | DUPLICATE | SUPPLEMENTARY, unique=True
        )
    )
    library: Library = 1
    conversion: str = "T>C"
    min_quality: int = 20

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "labeling",
            "name": self.name,
            "filter": self.filter._native(),
            "library": self.library,
            "conversion": self.conversion,
            "min_quality": self.min_quality,
        }


Features = Literal["genes", "gene-full", "promoters", "gene-activity", "peaks", "bins"]
Context = Literal["CpG", "CHG", "CHH"]

JUICER_RESOLUTIONS = (
    2_500_000,
    1_000_000,
    500_000,
    250_000,
    100_000,
    50_000,
    25_000,
    10_000,
    5_000,
)
"""The base-pair resolutions the Juicer tools' ``pre`` bins at."""


@dataclass(frozen=True, slots=True)
class MethylMatrix:
    """Methylation calls of ``contexts`` per region × row, written as a 10x
    directory ``<methylation>/<name>/`` with the layers ``methylated`` and
    ``total`` (``methylated.cpk``, ``total.cpk``).

    ``features`` and their options are those of :class:`Counts` (genes,
    promoters, gene activity, peaks, bins); ``rows`` is a tag (``"CB"``),
    ``"@RG"`` or ``"@file"`` (one row per BAM: single-cell bisulfite with a
    BAM per cell).
    """

    name: str
    features: Features
    rows: str = "@file"
    whitelist: PathLike | None = None
    contexts: Sequence[Context] = ("CpG",)
    upstream: int = 0
    downstream: int = 0
    peaks: PathLike | None = None
    bin_size: int | None = None
    biotypes: Sequence[str] | None = None

    def _native(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "features": self.features,
            "rows": self.rows,
            "whitelist": _path(self.whitelist),
            "contexts": list(self.contexts),
            "upstream": self.upstream,
            "downstream": self.downstream,
            "peaks": _path(self.peaks),
            "bin_size": self.bin_size,
            "biotypes": None if self.biotypes is None else list(self.biotypes),
        }


@dataclass(frozen=True, slots=True)
class Methylation:
    """DNA methylation (WGBS, EM-seq, single-cell bisulfite, long reads)
    under ``<name>/``: ``cytosines.tsv.gz`` (+ ``.tbi``; chrom, start, end,
    context, strand, methylated, unmethylated), ``beta.bw``,
    ``coverage.bw``, ``mbias.tsv`` (also in ``summary.json``) and the
    region matrices of ``matrices``.

    ``source``: ``"bismark"`` (``XM`` / ``XG`` tags), ``"modifications"``
    (``MM`` / ``ML``, 5mC ``C+m``, called at ``min_probability``),
    ``"reference"`` (the read against the task's ``reference`` or its ``MD``
    tag: C>T on the original top strand, G>A on the bottom; the strand from
    ``YD``, ``XG`` or read 1's orientation), or ``"auto"`` (per record, in
    that order of tags). Bases below ``min_quality`` and positions of the
    task's ``variants`` are not called. With ``overlap``, read 2's calls
    inside read 1's span are dropped. ``merge_cpg`` writes one line per CpG
    (both strands). The per-cytosine outputs need a single BAM; turn
    ``table`` and ``bigwig`` off to count matrices over several.
    """

    name: str = "methylation"
    source: Literal["auto", "bismark", "reference", "modifications"] = "auto"
    contexts: Sequence[Context] = ("CpG",)
    merge_cpg: bool = False
    filter: Filter = field(
        default_factory=lambda: Filter(
            exclude=UNMAPPED | SECONDARY | QC_FAIL | DUPLICATE | SUPPLEMENTARY, min_mapq=10
        )
    )
    min_quality: int = 5
    min_probability: float = 0.5
    overlap: bool = True
    table: bool = True
    bigwig: bool = True
    matrices: Sequence[MethylMatrix] = ()

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "methylation",
            "name": self.name,
            "source": self.source,
            "contexts": list(self.contexts),
            "merge_cpg": self.merge_cpg,
            "filter": self.filter._native(),
            "min_quality": self.min_quality,
            "min_probability": self.min_probability,
            "overlap": self.overlap,
            "table": self.table,
            "bigwig": self.bigwig,
            "matrices": [m._native() for m in self.matrices],
        }


@dataclass(frozen=True, slots=True)
class Contacts:
    """Hi-C / Micro-C contacts under ``<name>/``: ``contacts.pairs.gz`` (4DN
    pairs, header with ``#chromsize`` lines) and ``contacts.hic`` (Juicer,
    raw counts) at ``resolutions`` (none: no ``.hic``), with per-barcode
    ``cell_metrics.csv`` (and ``summary.json`` entries) by ``barcode_tag``.

    Needs a coordinate-sorted BAM with mate tags (``MC``, ``MQ``: ``samtools
    fixmate -m``): each pair is taken from its leftmost mate, its ends the
    two 5' positions. Both mates need MAPQ ``filter.min_mapq``; within a
    chromosome, contacts closer than ``min_distance`` are dropped, and with
    ``self_ligation`` inward and outward ones closer than it; ``dedup``
    removes contacts with the same ends, strands and barcode.
    """

    name: str = "hic"
    filter: Filter = field(
        default_factory=lambda: Filter(
            require=PAIRED,
            exclude=UNMAPPED | MATE_UNMAPPED | SECONDARY | QC_FAIL | SUPPLEMENTARY | DUPLICATE,
            min_mapq=30,
        )
    )
    min_distance: int = 0
    self_ligation: int = 0
    dedup: bool = True
    pairs: bool = True
    resolutions: Sequence[int] = JUICER_RESOLUTIONS
    genome: str = "unknown"
    barcode_tag: str | None = None
    whitelist: PathLike | None = None
    sidecar: bool = False

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "hic",
            "name": self.name,
            "filter": self.filter._native(),
            "min_distance": self.min_distance,
            "self_ligation": self.self_ligation,
            "dedup": self.dedup,
            "pairs": self.pairs,
            "resolutions": list(self.resolutions),
            "genome": self.genome,
            "barcode_tag": self.barcode_tag,
            "whitelist": _path(self.whitelist),
            "sidecar": self.sidecar,
        }


@dataclass(frozen=True, slots=True)
class AlleleGenes:
    """Gene-level allele-specific expression: units counted per gene as
    ``hap1``, ``hap2`` or ``both`` by the alleles they show at the phased
    heterozygous sites of ``sample`` (default: the panel's first)."""

    sample: str | None = None
    library: Library = 1

    def _native(self) -> tuple[str | None, int]:
        return (self.sample, self.library)


@dataclass(frozen=True, slots=True)
class Alleles:
    """Allele counts at the biallelic SNVs of ``sites`` (BCF with its CSI;
    cellsnp-lite / vartrix): a 10x directory ``<name>/`` of sites × rows with
    the layers ``ref``, ``alt`` and ``other`` (``ref.cpk``, ...),
    ``variants.tsv.gz`` (per site: REF / ALT / other units, ALT fraction),
    and with ``genes`` (needs an annotation) ``<name>/genes/``.

    ``unit``: ``"molecule"`` (UMI-based) or
    ``"fragment"`` (read-name-based: a pair once), with the identity
    encoding documented in :mod:`omicsweave.quant`. Bases below
    ``min_quality`` are not counted. :meth:`bulk` gives the per-BAM
    defaults.
    """

    sites: PathLike
    name: str = "alleles"
    rows: str = "CB"
    whitelist: PathLike | None = None
    unit: Literal["molecule", "fragment"] = "molecule"
    umi_tag: str = "UB"
    filter: Filter = field(
        default_factory=lambda: Filter(exclude=UNMAPPED | SECONDARY | QC_FAIL, min_mapq=20)
    )
    min_quality: int = 20
    table: bool = True
    genes: AlleleGenes | None = None

    @staticmethod
    def bulk(
        sites: PathLike, *, name: str = "alleles", genes: AlleleGenes | None = None
    ) -> Alleles:
        """Read pairs per input BAM, duplicates excluded."""
        return Alleles(
            sites,
            name=name,
            rows="@file",
            unit="fragment",
            filter=Filter(exclude=UNMAPPED | SECONDARY | QC_FAIL | DUPLICATE, min_mapq=20),
            genes=genes,
        )

    def _native(self) -> dict[str, Any]:
        return {
            "kind": "alleles",
            "name": self.name,
            "sites": os.fspath(self.sites),
            "rows": self.rows,
            "whitelist": _path(self.whitelist),
            "unit": self.unit,
            "umi_tag": self.umi_tag,
            "filter": self.filter._native(),
            "min_quality": self.min_quality,
            "table": self.table,
            "genes": None if self.genes is None else self.genes._native(),
        }


Product = Coverage | Counts | Atac | Junctions | Labeling | Methylation | Contacts | Alleles
