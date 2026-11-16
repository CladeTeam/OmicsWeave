"""What one pass over the BAMs computes: each product a frozen dataclass.

Defaults and presets are the native layer's, read once at import; the
native layer reads a product back attribute by attribute and checks every
value, so a dataclass holds options, not rules.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from typing import TYPE_CHECKING, Any, ClassVar, Literal, Self

from .. import _native
from ..formats import PathLike
from ._features import Features, GeneActivity, from_native

if TYPE_CHECKING:
    from ..annotation import Selection

__all__ = [
    "AlleleGenes",
    "Alleles",
    "Atac",
    "Contacts",
    "Context",
    "Counts",
    "Coverage",
    "Filter",
    "JunctionMatrix",
    "Junctions",
    "Labeled",
    "Labeling",
    "MethylMatrix",
    "Methylation",
    "Product",
    "Strandedness",
]

Strandedness = Literal["unstranded", "forward", "reverse"]
"""How reads relate to the transcript's strand (featureCounts ``-s`` 0, 1,
2): ``"forward"``: read 1 on the transcript's strand (read 2 counts as read
1's); ``"reverse"``: read 1 on the opposite strand."""

Context = Literal["CpG", "CHG", "CHH"]

_D: Mapping[str, Mapping[str, Any]] = _native.quant_defaults()


@dataclass(frozen=True, slots=True)
class Filter:
    """Which records a product counts, by FLAG bits
    (:class:`~omicsweave.types.SamFlag`), MAPQ and the ``NH`` tag."""

    require: int = _D["filter"]["require"]
    """FLAG bits a record must all have."""
    exclude: int = _D["filter"]["exclude"]
    """FLAG bits a record must have none of (default: unmapped, secondary,
    QC-failed and duplicate, as ``samtools depth``)."""
    min_mapq: int = _D["filter"]["min_mapq"]
    unique: bool = _D["filter"]["unique"]
    """Only records whose ``NH`` tag is 1 (none counts as 1)."""


def _values(preset: Mapping[str, Any]) -> dict[str, Any]:
    """A preset as dataclass keyword arguments."""
    out = dict(preset)
    if "filter" in out:
        out["filter"] = Filter(**out["filter"])
    if "features" in out:
        out["features"] = from_native(out["features"])
    if out.get("gene_activity") is not None:
        out["gene_activity"] = from_native(out["gene_activity"])
    if "contexts" in out:
        out["contexts"] = tuple(out["contexts"])
    if "resolutions" in out:
        out["resolutions"] = tuple(out["resolutions"])
    return out


_COVERAGE = _values(_D["coverage"])


@dataclass(frozen=True, slots=True)
class Coverage:
    """Coverage as bigWig: ``<name>.bw``, or ``<name>.forward.bw`` and
    ``<name>.reverse.bw`` with ``strands="split"``.

    - ``unit``: ``"bases"`` (aligned ``M`` / ``=`` / ``X`` bases, as
      ``samtools depth``; with ``deletions`` also ``D``, ``samtools depth
      -J``), ``"span"``, ``"five-prime"``, ``"three-prime"``, or
      ``"fragment"`` (a pair once, from its leftmost mate, ``[pos, pos +
      |TLEN|)``, up to ``max_fragment`` bases when not 0).
    - ``strand``: with split strands, the transcript strand a read stands
      for (:data:`Strandedness`); ``"unstranded"`` takes the read's own.
    - ``scale``: ``"none"``, ``"cpm"`` (per million counted reads or
      fragments) or a factor.
    """

    name: str = _COVERAGE["name"]
    _: KW_ONLY
    unit: Literal["bases", "span", "five-prime", "three-prime", "fragment"] = _COVERAGE["unit"]
    deletions: bool = _COVERAGE["deletions"]
    max_fragment: int = _COVERAGE["max_fragment"]
    strands: Literal["merged", "split"] = _COVERAGE["strands"]
    strand: Strandedness = _COVERAGE["strand"]
    filter: Filter = _COVERAGE["filter"]
    scale: Literal["none", "cpm"] | float = _COVERAGE["scale"]
    sidecar: bool = _COVERAGE["sidecar"]
    kind: ClassVar[str] = "coverage"


_LABELED = _values(_D["labeled"])


@dataclass(frozen=True, slots=True)
class Labeled:
    """A labeled layer of a gene matrix: units (reads, pairs, UMIs) with at
    least ``min_conversions`` conversions.

    ``conversion`` is in transcript sense (``"T>C"`` for 4sU: on the genome
    T>C for forward transcripts, A>G for reverse ones); bases below
    ``min_quality`` and positions of the task's ``variants`` are not
    counted. The reference is the task's ``reference``, else each read's
    ``MD`` tag. A pair's conversions are those of both mates; a molecule is
    labeled when any of its reads is. Under an unstranded library, the
    strand is the assigned gene's.
    """

    conversion: str = _LABELED["conversion"]
    min_quality: int = _LABELED["min_quality"]
    min_conversions: int = _LABELED["min_conversions"]


_COUNTS = _values(_D["counts"])


@dataclass(frozen=True, slots=True)
class Counts:
    """A count matrix, rows × ``features``, written as a 10x directory
    ``<name>/`` (``barcodes.tsv.gz``, ``features.tsv.gz`` and the matrix).

    - ``features``: :class:`GeneExons`, :class:`GeneSpans`,
      :class:`Promoters`, :class:`GeneActivity`, :class:`Peaks`,
      :class:`Bins` or :class:`Regions`; ``subset`` (an
      :meth:`Annotation.subset <omicsweave.annotation.Annotation.subset>` of
      the task's annotation) keeps its genes, and for gene features their
      selected transcripts.
    - ``unit``: ``"read"``, ``"pair"`` (both mates together, featureCounts
      ``-p --countReadPairs``), ``"fragment"`` / ``"tn5-fragment"`` (a
      pair's span from its leftmost mate, up to ``max_fragment`` bases when
      not 0), ``"cut-site"`` (Tn5 insertions).
    - ``dedup``: ``"none"``, ``"umi"`` (molecules: units of a row with the
      same UMI, by ``umi_tag``, counted once) or ``"position"``
      (fragments of a row with the same ends counted once).
    - ``rows``: a tag (``"CB"``), ``"@RG"`` (read groups) or ``"@file"``
      (one row per input BAM); ``whitelist`` restricts and orders them.
    - A unit overlapping more than one feature is ambiguous unless
      ``multi_overlap``.
    - ``table``: also write a featureCounts-style ``counts.tsv``.
    - ``velocity``: add the ``spliced``, ``unspliced`` and ``ambiguous``
      layers (gene features; velocyto / STARsolo ``Velocyto`` rules),
      written as ``spliced.cpk`` etc. beside ``matrix.cpk``.
    - ``labeled``: add a ``labeled`` layer (:class:`Labeled`).

    The presets :meth:`rna_single_cell`, :meth:`rna_bulk`, :meth:`atac` and
    :meth:`chip` fill in what a library type usually needs. UMI and
    molecule identity follows the rules documented in
    :mod:`omicsweave.quant`.
    """

    features: Features
    _: KW_ONLY
    name: str = _COUNTS["name"]
    unit: Literal["read", "pair", "fragment", "tn5-fragment", "cut-site"] = _COUNTS["unit"]
    dedup: Literal["none", "umi", "position"] = _COUNTS["dedup"]
    rows: str = _COUNTS["rows"]
    whitelist: PathLike | None = None
    umi_tag: str = _COUNTS["umi_tag"]
    strand: Strandedness = _COUNTS["strand"]
    multi_overlap: bool = _COUNTS["multi_overlap"]
    filter: Filter = _COUNTS["filter"]
    max_fragment: int = _COUNTS["max_fragment"]
    table: bool = _COUNTS["table"]
    subset: Selection | None = None
    velocity: bool = _COUNTS["velocity"]
    labeled: Labeled | None = None
    kind: ClassVar[str] = "counts"

    @classmethod
    def _preset(cls, key: str, features: Features | None, options: dict[str, Any]) -> Self:
        values = _values(_D[key])
        if features is not None:
            values["features"] = features
        return cls(**(values | options))

    @classmethod
    def rna_single_cell(cls, features: Features | None = None, **options: Any) -> Self:
        """Single-cell RNA (STARsolo ``Gene``, Cell Ranger): genes × UMIs ×
        ``CB``, the raw matrix over every barcode, unique alignments,
        forward-stranded. ``options`` override any field."""
        return cls._preset("counts.single_cell", features, options)

    @classmethod
    def rna_bulk(cls, features: Features | None = None, **options: Any) -> Self:
        """Bulk RNA (featureCounts defaults): genes × reads (``unit="pair"``:
        read pairs) × input file, unstranded, every mapped record with
        ``NH`` 1 counted, with a ``counts.tsv`` table."""
        return cls._preset("counts.bulk", features, options)

    @classmethod
    def atac(cls, features: Features, **options: Any) -> Self:
        """Single-cell ATAC: ``features`` × Tn5 cut sites × ``CB``, duplicate
        fragments counted once, a cut site counted in every feature."""
        return cls._preset("counts.atac", features, options)

    @classmethod
    def chip(cls, features: Features, **options: Any) -> Self:
        """ChIP-seq / CUT&Tag: ``features`` (bins or peaks) × fragments ×
        input file, duplicates excluded by flag, with a ``counts.tsv``
        table."""
        return cls._preset("counts.chip", features, options)


_ATAC = _values(_D["atac"])


@dataclass(frozen=True, slots=True)
class Atac:
    """Single-cell ATAC under ``<name>/``: ``fragments.tsv.gz`` (+ ``.tbi``),
    ``barcode_metrics.csv``, and cut-site matrices over ``peaks`` (bigBed)
    and ``tiles``, and fragments over ``gene_activity`` windows when the
    task has an annotation.

    A fragment is a proper pair from its leftmost mate, Tn5-shifted (``pos +
    4`` to ``pos + |TLEN| - 5``); duplicates of a barcode collapse, the
    fifth column counting them.
    """

    name: str = _ATAC["name"]
    _: KW_ONLY
    barcode_tag: str = _ATAC["barcode_tag"]
    whitelist: PathLike | None = None
    filter: Filter = _ATAC["filter"]
    max_fragment: int = _ATAC["max_fragment"]
    fragments: bool = _ATAC["fragments"]
    peaks: PathLike | None = None
    tiles: int | None = _ATAC["tiles"]
    gene_activity: GeneActivity | None = _ATAC["gene_activity"]
    kind: ClassVar[str] = "atac"


@dataclass(frozen=True, slots=True)
class JunctionMatrix:
    """Junctions × rows, beside the junction table: one column per junction
    of the table. ``rows`` is a tag (``"CB"``), ``"@RG"`` or ``"@file"``;
    ``dedup`` ``"umi"`` counts UMIs (STARsolo ``SJ``), ``"none"`` reads or
    pairs; ``filter`` applies on top of the table's. The defaults are
    :meth:`single_cell`'s."""

    _: KW_ONLY
    rows: str = _D["junction_matrix.single_cell"]["rows"]
    whitelist: PathLike | None = None
    umi_tag: str = _D["junction_matrix.single_cell"]["umi_tag"]
    dedup: Literal["none", "umi"] = _D["junction_matrix.single_cell"]["dedup"]
    filter: Filter = _values(_D["junction_matrix.single_cell"])["filter"]

    @classmethod
    def single_cell(cls, **options: Any) -> Self:
        """UMIs per ``CB``, unique alignments (STARsolo ``SJ``)."""
        return cls(**(_values(_D["junction_matrix.single_cell"]) | options))

    @classmethod
    def bulk(cls, **options: Any) -> Self:
        """Reads (pairs) per input file, unique alignments."""
        return cls(**(_values(_D["junction_matrix.bulk"]) | options))


_JUNCTIONS = _values(_D["junctions"])


@dataclass(frozen=True, slots=True)
class Junctions:
    """Splice junctions as ``<name>/junctions.tab``, the nine columns of
    STAR's ``SJ.out.tab``: chromosome, first and last intron base (1-based),
    strand (0, 1 ``+``, 2 ``-``), motif (0 non-canonical, 1 GT/AG, 2 CT/AC,
    3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT; from the task's ``reference``),
    annotated, unique reads, multi-mapping reads, largest overhang.

    A pair crossing a junction counts once (``pairs``); every alignment of a
    multi-mapper counts. With ``matrix``, ``<name>/`` is also a 10x
    directory of junctions × rows.
    """

    name: str = _JUNCTIONS["name"]
    _: KW_ONLY
    filter: Filter = _JUNCTIONS["filter"]
    pairs: bool = _JUNCTIONS["pairs"]
    matrix: JunctionMatrix | None = None
    kind: ClassVar[str] = "junctions"


_LABELING = _values(_D["labeling"])


@dataclass(frozen=True, slots=True)
class Labeling:
    """Metabolic labeling statistics under ``<name>/``:
    ``conversions.tsv.gz`` (reads by coverage of the conversion's source
    base and number of conversions) and ``mismatches.tsv`` (per mate,
    sequencing cycle and transcript-sense reference base, the read bases),
    with overall substitution rates in ``summary.json``.

    A read's transcript strand follows ``strand``; ``"unstranded"`` reads
    as ``"forward"``. Reference bases come from the task's ``reference``,
    else the ``MD`` tag; bases below ``min_quality`` and positions of the
    task's ``variants`` are skipped.
    """

    name: str = _LABELING["name"]
    _: KW_ONLY
    filter: Filter = _LABELING["filter"]
    strand: Strandedness = _LABELING["strand"]
    conversion: str = _LABELING["conversion"]
    min_quality: int = _LABELING["min_quality"]
    kind: ClassVar[str] = "labeling"


_METHYL_MATRIX = _values(_D["methyl_matrix"])


@dataclass(frozen=True, slots=True)
class MethylMatrix:
    """Methylation calls of ``contexts`` per feature × row, written as a 10x
    directory ``<methylation>/<name>/`` with the layers ``methylated`` and
    ``total`` (``methylated.cpk``, ``total.cpk``).

    ``features`` and ``subset`` are those of :class:`Counts`; ``rows`` is
    a tag (``"CB"``), ``"@RG"`` or ``"@file"`` (one row per BAM:
    single-cell bisulfite with a BAM per cell).
    """

    name: str
    features: Features
    _: KW_ONLY
    rows: str = _METHYL_MATRIX["rows"]
    whitelist: PathLike | None = None
    contexts: Sequence[Context] = _METHYL_MATRIX["contexts"]
    subset: Selection | None = None


_METHYLATION = _values(_D["methylation"])


@dataclass(frozen=True, slots=True)
class Methylation:
    """DNA methylation (WGBS, EM-seq, single-cell bisulfite, long reads)
    under ``<name>/``: ``cytosines.tsv.gz`` (+ ``.tbi``; chrom, start, end,
    context, strand, methylated, unmethylated), ``beta.bw``,
    ``coverage.bw``, ``mbias.tsv`` (also in ``summary.json``) and the
    feature matrices of ``matrices``.

    ``source``: ``"bismark"`` (``XM`` / ``XG`` tags), ``"modifications"``
    (``MM`` / ``ML``, 5mC ``C+m``, called at ``min_probability``),
    ``"reference"`` (the read against the task's ``reference`` or its
    ``MD`` tag: C>T on the original top strand, G>A on the bottom; the
    strand from ``YD``, ``XG`` or read 1's orientation), or ``"auto"`` (per
    record, in that order of tags). Bases below ``min_quality`` and
    positions of the task's ``variants`` are not called. With ``overlap``,
    read 2's calls inside read 1's span are dropped. ``merge_cpg`` writes
    one line per CpG (both strands). The per-cytosine outputs need a single
    BAM; turn ``table`` and ``bigwig`` off to count matrices over several.
    """

    name: str = _METHYLATION["name"]
    _: KW_ONLY
    source: Literal["auto", "bismark", "reference", "modifications"] = _METHYLATION["source"]
    contexts: Sequence[Context] = _METHYLATION["contexts"]
    merge_cpg: bool = _METHYLATION["merge_cpg"]
    filter: Filter = _METHYLATION["filter"]
    min_quality: int = _METHYLATION["min_quality"]
    min_probability: float = _METHYLATION["min_probability"]
    overlap: bool = _METHYLATION["overlap"]
    table: bool = _METHYLATION["table"]
    bigwig: bool = _METHYLATION["bigwig"]
    matrices: Sequence[MethylMatrix] = ()
    kind: ClassVar[str] = "methylation"


_CONTACTS = _values(_D["contacts"])


@dataclass(frozen=True, slots=True)
class Contacts:
    """Hi-C / Micro-C contacts under ``<name>/``: ``contacts.pairs.gz`` (4DN
    pairs, header with ``#chromsize`` lines) and ``contacts.hic`` (Juicer,
    raw counts) at ``resolutions`` (empty: no ``.hic``; the default is
    Juicer ``pre``'s), with per-barcode ``cell_metrics.csv`` (and
    ``summary.json`` entries) by ``barcode_tag``.

    Needs a coordinate-sorted BAM with mate tags (``MC``, ``MQ``:
    ``samtools fixmate -m``): each pair is taken from its leftmost mate,
    its ends the two 5' positions. Both mates need MAPQ
    ``filter.min_mapq``; within a chromosome, contacts closer than
    ``min_distance`` are dropped, and with ``self_ligation`` inward and
    outward ones closer than it; ``dedup`` removes contacts with the same
    ends, strands and barcode. ``sidecar`` also writes the ``.hic``'s
    ``.sidecar``.
    """

    name: str = _CONTACTS["name"]
    _: KW_ONLY
    filter: Filter = _CONTACTS["filter"]
    min_distance: int = _CONTACTS["min_distance"]
    self_ligation: int = _CONTACTS["self_ligation"]
    dedup: bool = _CONTACTS["dedup"]
    pairs: bool = _CONTACTS["pairs"]
    resolutions: Sequence[int] = _CONTACTS["resolutions"]
    genome: str = _CONTACTS["genome"]
    barcode_tag: str | None = _CONTACTS["barcode_tag"]
    whitelist: PathLike | None = None
    sidecar: bool = _CONTACTS["sidecar"]
    kind: ClassVar[str] = "contacts"


@dataclass(frozen=True, slots=True)
class AlleleGenes:
    """Gene-level allele-specific expression: units counted per gene as
    ``hap1``, ``hap2`` or ``both`` by the alleles they show at the phased
    heterozygous sites of ``sample`` (default: the panel's first). ``strand``
    relates reads to the genes' strand: ``"forward"`` for 10x libraries,
    usually ``"unstranded"`` for bulk ones."""

    sample: str | None = None
    strand: Strandedness = "forward"


_ALLELES = _values(_D["alleles.single_cell"])


@dataclass(frozen=True, slots=True)
class Alleles:
    """Allele counts at the biallelic SNVs of ``sites`` (BCF with its CSI;
    cellsnp-lite / vartrix): a 10x directory ``<name>/`` of sites × rows
    with the layers ``ref``, ``alt`` and ``other`` (``ref.cpk``, ...),
    ``variants.tsv.gz`` (per site: REF / ALT / other units, ALT fraction),
    and with ``genes`` (needs an annotation) ``<name>/genes/``.

    ``unit``: ``"molecule"`` (UMI-based, by ``umi_tag``) or ``"fragment"``
    (read-name-based: a pair once), with the identity encoding documented
    in :mod:`omicsweave.quant`. Bases below ``min_quality`` are not
    counted. The defaults are :meth:`single_cell`'s.
    """

    sites: PathLike
    _: KW_ONLY
    name: str = _ALLELES["name"]
    rows: str = _ALLELES["rows"]
    whitelist: PathLike | None = None
    unit: Literal["molecule", "fragment"] = _ALLELES["unit"]
    umi_tag: str = _ALLELES["umi_tag"]
    filter: Filter = _ALLELES["filter"]
    min_quality: int = _ALLELES["min_quality"]
    table: bool = _ALLELES["table"]
    genes: AlleleGenes | None = None
    kind: ClassVar[str] = "alleles"

    @classmethod
    def single_cell(cls, sites: PathLike, **options: Any) -> Self:
        """Sites × ``CB``, UMIs, primary mapped records that are not QC
        failures, MAPQ 20 or more (cellsnp-lite)."""
        return cls(sites, **(_values(_D["alleles.single_cell"]) | options))

    @classmethod
    def bulk(cls, sites: PathLike, **options: Any) -> Self:
        """Sites × input file, read pairs by name, duplicates excluded too."""
        return cls(sites, **(_values(_D["alleles.bulk"]) | options))


Product = Coverage | Counts | Atac | Junctions | Labeling | Methylation | Contacts | Alleles
