"""``ow quant BAM... -o DIR``: coverage, count matrices, ATAC fragments, splice
junctions, metabolic-labeling statistics, DNA methylation, Hi-C contacts and
allele counts in one pass over coordinate-sorted BAMs.

Only the command line lives here: arguments, progress and messages. The
quantification, its presets included, is :mod:`omicsweave.quant`. An option
left unset keeps the preset of each product; the read-selection options
(``--min-mapq``, ``--require-flags``, ...), ``--library``, ``--umi-tag``,
``--barcode-tag`` and ``--whitelist`` apply to every product that has them.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from typing import Any

from .. import quant
from ..annotation import Annotation, read_aliases
from ._common import add_budget, budget_of, human, print_budget, progress_for

__all__ = ["register"]

_LIBRARY: dict[int, quant.Strandedness] = {0: "unstranded", 1: "forward", 2: "reverse"}


def _flags(text: str) -> int:
    try:
        return int(text, 0)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a FLAG value: {text!r} (e.g. 0x704)") from None


def _scale(text: str) -> str | float:
    if text in ("none", "cpm"):
        return text
    try:
        return float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a scale: {text!r} (none, cpm or a number)") from None


def _flank(text: str) -> tuple[int, int] | tuple[()]:
    if text == "none":
        return ()
    try:
        up, down = (int(x) for x in text.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError(f"not UP,DOWN: {text!r}") from None
    return (up, down)


def _ints(text: str) -> tuple[int, ...]:
    if text == "none":
        return ()
    try:
        return tuple(int(x) for x in text.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError(f"not N,N,...: {text!r}") from None


def _names(text: str) -> tuple[str, ...]:
    return tuple(x.strip() for x in text.split(",") if x.strip())


def register(commands: Any) -> None:
    p = commands.add_parser(
        "quant",
        help="coverage, count matrices, ATAC, junctions, labeling, methylation, Hi-C and alleles",
        description=(
            "Quantify coordinate-sorted BAMs in one pass, in parallel within a memory budget.\n\n"
            "Every product asked for is computed from the same scan and written under\n"
            "DIR: coverage as bigWig, count matrices as 10x directories (barcodes.tsv.gz,\n"
            "features.tsv.gz, matrix.cpk | matrix.mtx.gz | matrix.csv.gz, and layers such\n"
            "as spliced.cpk or labeled.cpk), ATAC fragments as fragments.tsv.gz with a\n"
            "tabix index, splice junctions as junctions.tab (STAR SJ.out.tab columns),\n"
            "labeling statistics, methylation (cytosines.tsv.gz, beta/coverage bigWigs,\n"
            "region matrices), Hi-C (contacts.pairs.gz, contacts.hic), allele counts\n"
            "(ref/alt/other matrices at SNVs), and summary.json. Files appear only\n"
            "once the whole run succeeded. Unsorted BAMs are refused: sort them first\n"
            "with `ow convert in.bam sorted.bam --sort`.\n\n"
            "Options left unset keep each product's preset; the shared options below\n"
            "apply to every product that has them."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("inputs", metavar="BAM", nargs="+", help="coordinate-sorted BAM files")
    p.add_argument("-o", "--out", required=True, metavar="DIR", help="the output directory")
    p.add_argument("--gtf", metavar="FILE", help="gene annotation (GTF / GFF3) for gene products")
    p.add_argument(
        "--aliases",
        metavar="FILE",
        help="chromosome aliases, two names per line (beyond chr1/1 and chrM/MT)",
    )
    p.add_argument(
        "--matrix",
        choices=("cpk", "mtx", "csv"),
        default="cpk",
        help="how count matrices are stored (default: %(default)s, CellPack)",
    )
    p.add_argument(
        "--reference",
        metavar="FILE",
        help="the genome the BAMs were aligned to (bgzip FASTA with .fai/.gzi, or 2bit): "
        "labeling and methylation reference bases (else MD tags) and junction motifs",
    )
    p.add_argument(
        "--variants",
        metavar="FILE.bcf",
        help="known variant sites (BCF with its .csi), masked when counting conversions "
        "and methylation calls",
    )

    shared = p.add_argument_group("shared (every product that has them; default: its preset)")
    shared.add_argument(
        "--require-flags", type=_flags, metavar="FLAGS", help="records must have these FLAG bits"
    )
    shared.add_argument(
        "--exclude-flags", type=_flags, metavar="FLAGS", help="skip records with any of these bits"
    )
    shared.add_argument("--min-mapq", type=int, metavar="Q", help="skip records below MAPQ Q")
    shared.add_argument(
        "--unique", action=argparse.BooleanOptionalAction, help="only records with NH 1"
    )
    shared.add_argument(
        "--library",
        type=int,
        choices=(0, 1, 2),
        help="strandedness as featureCounts -s: 0 unstranded, 1 forward, 2 reverse",
    )
    shared.add_argument("--barcode-tag", metavar="TAG", help="the cell barcode tag (CB)")
    shared.add_argument("--umi-tag", metavar="TAG", help="the UMI tag (UB)")
    shared.add_argument("--whitelist", metavar="FILE", help="barcodes to keep, one per line")
    shared.add_argument(
        "--max-fragment", type=int, metavar="N", help="longest fragment counted (0: any)"
    )
    shared.add_argument(
        "--min-base-quality",
        type=int,
        metavar="Q",
        help="labeling, labeled layers and alleles: bases below Q are not counted",
    )
    shared.add_argument(
        "--sidecar", action="store_true", help="also write the .sidecar of bigWigs and .hic"
    )

    cov = p.add_argument_group("coverage (bigWig)")
    cov.add_argument("--coverage", action="store_true", help="write coverage.bw")
    cov.add_argument(
        "--coverage-unit",
        choices=("bases", "bases+deletions", "span", "five-prime", "three-prime", "fragment"),
        help="what a read adds (default: bases, as samtools depth)",
    )
    cov.add_argument(
        "--strands", choices=("merged", "split"), help="one track, or .forward.bw and .reverse.bw"
    )
    cov.add_argument("--scale", type=_scale, help="none, cpm, or a factor (default: none)")

    genes = p.add_argument_group("gene counts (needs --gtf)")
    genes.add_argument(
        "--genes",
        choices=("single-cell", "bulk"),
        help="single-cell: genes x UMIs x barcodes (STARsolo Gene); "
        "bulk: genes x reads x BAM (featureCounts)",
    )
    genes.add_argument(
        "--gene-span",
        choices=("exons", "full"),
        default="exons",
        help="count exons, or whole genes with introns (GeneFull)",
    )
    genes.add_argument("--pairs", action="store_true", help="bulk: count read pairs")
    genes.add_argument(
        "--dedup", choices=("none", "umi", "position"), help="how units of a row collapse"
    )
    genes.add_argument("--biotypes", type=_names, metavar="T,T,...", help="only these biotypes")
    genes.add_argument(
        "--multi-overlap", action="store_true", help="count in every feature a unit overlaps"
    )
    genes.add_argument(
        "--velocity",
        action="store_true",
        help="add spliced / unspliced / ambiguous layers (velocyto rules)",
    )
    genes.add_argument(
        "--labeled",
        type=int,
        metavar="N",
        help="add a labeled layer: reads / pairs / UMIs with at least N conversions",
    )

    splice = p.add_argument_group("splice junctions")
    splice.add_argument(
        "--junctions", action="store_true", help="junctions/junctions.tab (STAR SJ.out.tab)"
    )
    splice.add_argument(
        "--junction-matrix",
        choices=("single-cell", "bulk"),
        help="also junctions x barcodes (UMIs, STARsolo SJ) or x BAM (reads)",
    )
    splice.add_argument(
        "--single-mates",
        action="store_true",
        help="count each mate crossing a junction (default: a pair once)",
    )

    label = p.add_argument_group("metabolic labeling (SLAM-seq, 4sU)")
    label.add_argument(
        "--labeling", action="store_true", help="labeling/conversions.tsv.gz and mismatches.tsv"
    )
    label.add_argument("--conversion", metavar="X>Y", help="in transcript sense (default: T>C)")

    atac = p.add_argument_group("single-cell ATAC")
    atac.add_argument("--atac", action="store_true", help="fragments, metrics and matrices")
    atac.add_argument("--peaks", metavar="FILE.bb", help="peaks (bigBed): a peak matrix")
    atac.add_argument("--tiles", type=int, metavar="N", help="tile size (default: 5000; 0: none)")
    atac.add_argument(
        "--gene-activity",
        type=_flank,
        metavar="UP,DOWN",
        help="gene activity flanks (default: 2000,0; none to skip)",
    )
    atac.add_argument("--no-fragments", action="store_true", help="skip fragments.tsv.gz")

    chip = p.add_argument_group("bulk ChIP / CUT&Tag (fragments x BAM)")
    chip.add_argument("--chip-bins", type=int, metavar="N", help="genome bins of N bases")
    chip.add_argument("--chip-peaks", metavar="FILE.bb", help="peaks (bigBed)")

    meth = p.add_argument_group("DNA methylation (bisulfite, EM-seq, MM/ML)")
    meth.add_argument(
        "--methylation", action="store_true", help="methylation/cytosines.tsv.gz, tracks, M-bias"
    )
    meth.add_argument(
        "--methylation-source",
        choices=("auto", "bismark", "reference", "modifications"),
        help="XM tags, MM/ML tags, or reads against --reference / MD (default: per record)",
    )
    meth.add_argument("--contexts", type=_names, metavar="C,C", help="CpG, CHG, CHH (default: CpG)")
    meth.add_argument("--merge-cpg", action="store_true", help="one line per CpG, both strands")
    meth.add_argument(
        "--methylation-min-quality", type=int, metavar="Q", help="call bases of quality Q or more"
    )
    meth.add_argument(
        "--min-probability",
        type=float,
        metavar="P",
        help="MM/ML: call 5mC at P, unmethylated at 1-P (default: 0.5)",
    )
    meth.add_argument("--keep-overlaps", action="store_true", help="count overlapping mates twice")
    meth.add_argument("--methylation-bins", type=int, metavar="N", help="a bins x rows matrix")
    meth.add_argument(
        "--methylation-genes", action="store_true", help="a gene-body x rows matrix (--gtf)"
    )
    meth.add_argument(
        "--methylation-promoters",
        type=_flank,
        metavar="UP,DOWN",
        help="a promoter x rows matrix (--gtf)",
    )
    meth.add_argument("--methylation-peaks", metavar="FILE.bb", help="a peak x rows matrix")
    meth.add_argument(
        "--methylation-rows", metavar="ROWS", help="matrix rows: a tag (CB), @RG or @file"
    )
    meth.add_argument(
        "--no-cytosines", action="store_true", help="matrices only (needed for several BAMs)"
    )

    hic = p.add_argument_group("Hi-C / Micro-C (needs MC and MQ tags: samtools fixmate -m)")
    hic.add_argument("--hic", action="store_true", help="hic/contacts.pairs.gz and contacts.hic")
    hic.add_argument(
        "--hic-resolutions",
        type=_ints,
        metavar="R,R",
        help="bin sizes of the .hic (default: the Juicer tools'; none: no .hic)",
    )
    hic.add_argument(
        "--min-distance", type=int, metavar="N", help="drop cis contacts closer than N"
    )
    hic.add_argument(
        "--self-ligation", type=int, metavar="N", help="drop inward/outward cis contacts below N"
    )
    hic.add_argument("--keep-duplicates", action="store_true", help="no contact deduplication")
    hic.add_argument("--no-pairs", action="store_true", help="skip contacts.pairs.gz")
    hic.add_argument("--genome", help="assembly name for pairs and .hic")
    hic.add_argument(
        "--hic-cells", action="store_true", help="per-barcode statistics by --barcode-tag (CB)"
    )

    al = p.add_argument_group("allele counts (cellsnp-lite / vartrix)")
    al.add_argument("--alleles", metavar="SITES.bcf", help="ref/alt/other matrices at these SNVs")
    al.add_argument(
        "--allele-mode",
        choices=("single-cell", "bulk"),
        default="single-cell",
        help="UMIs x barcodes, or read pairs x BAM (default: %(default)s)",
    )
    al.add_argument(
        "--ase",
        nargs="?",
        const="",
        metavar="SAMPLE",
        help="gene-level haplotype counts by SAMPLE's phased genotypes (default: the first)",
    )

    add_budget(p, "the output directory")
    p.add_argument(
        "-n", "--dry-run", action="store_true", help="check the task against the budget only"
    )
    p.add_argument("-q", "--quiet", action="store_true", help="print only the files written")
    p.set_defaults(run=run)


class _Shared:
    """The shared options of a command line, applied to a product."""

    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.filter = {
            k: v
            for k, v in (
                ("require", args.require_flags),
                ("exclude", args.exclude_flags),
                ("min_mapq", args.min_mapq),
                ("unique", args.unique),
            )
            if v is not None
        }
        self.strand = None if args.library is None else _LIBRARY[args.library]
        self.annotation = Annotation.read(args.gtf, attributes=None) if args.gtf else None
        self.subset = None
        if args.biotypes:
            if self.annotation is None:
                raise SystemExit("ow quant: error: --biotypes selects genes of --gtf")
            self.subset = self.annotation.subset(biotypes=args.biotypes)

    def apply(self, product: Any, **options: Any) -> Any:
        """``product`` with the shared read selection and ``options``; an
        option given as ``None`` keeps the product's value."""
        changes = {k: v for k, v in options.items() if v is not None}
        if self.filter and hasattr(product, "filter"):
            changes["filter"] = dataclasses.replace(product.filter, **self.filter)
        return dataclasses.replace(product, **changes)


def _labeled(args: argparse.Namespace) -> quant.Labeled | None:
    if args.labeled is None:
        return None
    options = {"min_conversions": args.labeled}
    if args.conversion is not None:
        options["conversion"] = args.conversion
    if args.min_base_quality is not None:
        options["min_quality"] = args.min_base_quality
    return quant.Labeled(**options)


def _coverage(args: argparse.Namespace, s: _Shared) -> quant.Coverage:
    unit = args.coverage_unit
    deletions = None
    if unit == "bases+deletions":
        unit, deletions = "bases", True
    return s.apply(
        quant.Coverage(),
        unit=unit,
        deletions=deletions,
        max_fragment=args.max_fragment if unit == "fragment" else None,
        strands=args.strands,
        strand=s.strand if args.strands == "split" else None,
        scale=args.scale,
        sidecar=args.sidecar or None,
    )


def _genes(args: argparse.Namespace, s: _Shared) -> quant.Counts:
    features = quant.GeneExons() if args.gene_span == "exons" else quant.GeneSpans()
    single = args.genes == "single-cell"
    base = quant.Counts.rna_single_cell if single else quant.Counts.rna_bulk
    return s.apply(
        base(features),
        rows=args.barcode_tag if single else None,
        whitelist=args.whitelist if single else None,
        umi_tag=args.umi_tag,
        strand=s.strand,
        unit="pair" if args.pairs else None,
        dedup=args.dedup,
        subset=s.subset,
        multi_overlap=args.multi_overlap or None,
        velocity=args.velocity or None,
        labeled=_labeled(args),
    )


def _junctions(args: argparse.Namespace, s: _Shared) -> quant.Junctions:
    matrix = None
    if args.junction_matrix == "single-cell":
        matrix = s.apply(
            quant.JunctionMatrix.single_cell(),
            rows=args.barcode_tag,
            whitelist=args.whitelist,
            umi_tag=args.umi_tag,
        )
    elif args.junction_matrix == "bulk":
        matrix = s.apply(quant.JunctionMatrix.bulk())
    return s.apply(quant.Junctions(), pairs=False if args.single_mates else None, matrix=matrix)


def _methylation(args: argparse.Namespace, s: _Shared) -> quant.Methylation:
    features: list[tuple[str, quant.Features]] = []
    if args.methylation_bins is not None:
        features.append(("bins", quant.Bins(args.methylation_bins)))
    if args.methylation_genes:
        features.append(("genes", quant.GeneSpans()))
    if args.methylation_promoters:
        up, down = args.methylation_promoters
        features.append(("promoters", quant.Promoters(up, down)))
    if args.methylation_peaks is not None:
        features.append(("peaks", quant.Peaks(args.methylation_peaks)))
    matrices = [
        s.apply(
            quant.MethylMatrix(name, f),
            rows=args.methylation_rows,
            whitelist=args.whitelist,
            contexts=args.contexts,
            subset=s.subset,
        )
        for name, f in features
    ]
    return s.apply(
        quant.Methylation(),
        source=args.methylation_source,
        contexts=args.contexts,
        merge_cpg=args.merge_cpg or None,
        min_quality=args.methylation_min_quality,
        min_probability=args.min_probability,
        overlap=False if args.keep_overlaps else None,
        table=False if args.no_cytosines else None,
        bigwig=False if args.no_cytosines else None,
        matrices=matrices,
    )


def _contacts(args: argparse.Namespace, s: _Shared) -> quant.Contacts:
    cells = args.hic_cells
    return s.apply(
        quant.Contacts(),
        min_distance=args.min_distance,
        self_ligation=args.self_ligation,
        dedup=False if args.keep_duplicates else None,
        pairs=False if args.no_pairs else None,
        resolutions=args.hic_resolutions,
        genome=args.genome,
        barcode_tag=(args.barcode_tag or "CB") if cells else None,
        whitelist=args.whitelist if cells else None,
        sidecar=args.sidecar or None,
    )


def _alleles(args: argparse.Namespace, s: _Shared) -> quant.Alleles:
    single = args.allele_mode == "single-cell"
    genes = None
    if args.ase is not None:
        strand = s.strand or ("forward" if single else "unstranded")
        genes = quant.AlleleGenes(sample=args.ase or None, strand=strand)
    if single:
        return s.apply(
            quant.Alleles.single_cell(args.alleles),
            rows=args.barcode_tag,
            whitelist=args.whitelist,
            umi_tag=args.umi_tag,
            min_quality=args.min_base_quality,
            genes=genes,
        )
    return s.apply(quant.Alleles.bulk(args.alleles), min_quality=args.min_base_quality, genes=genes)


def _products(args: argparse.Namespace, s: _Shared) -> list[quant.Product]:
    out: list[quant.Product] = []
    if args.coverage:
        out.append(_coverage(args, s))
    if args.genes:
        out.append(_genes(args, s))
    elif args.velocity or args.labeled is not None:
        raise SystemExit("ow quant: error: --velocity and --labeled are layers of --genes")
    if args.junctions or args.junction_matrix:
        out.append(_junctions(args, s))
    if args.labeling:
        out.append(
            s.apply(
                quant.Labeling(),
                strand=s.strand,
                conversion=args.conversion,
                min_quality=args.min_base_quality,
            )
        )
    if args.atac:
        atac = s.apply(
            quant.Atac(),
            barcode_tag=args.barcode_tag,
            whitelist=args.whitelist,
            max_fragment=args.max_fragment,
            fragments=False if args.no_fragments else None,
            peaks=args.peaks,
            tiles=args.tiles or None,
            gene_activity=quant.GeneActivity(*args.gene_activity) if args.gene_activity else None,
        )
        # 0 tiles and "none" switch those matrices off.
        if args.tiles == 0:
            atac = dataclasses.replace(atac, tiles=None)
        if args.gene_activity == ():
            atac = dataclasses.replace(atac, gene_activity=None)
        out.append(atac)
    if args.methylation:
        out.append(_methylation(args, s))
    if args.hic:
        out.append(_contacts(args, s))
    if args.alleles:
        out.append(_alleles(args, s))
    chip: list[tuple[str, quant.Features]] = []
    if args.chip_bins is not None:
        chip.append(("chip_bins", quant.Bins(args.chip_bins)))
    if args.chip_peaks is not None:
        chip.append(("chip_peaks", quant.Peaks(args.chip_peaks)))
    for name, features in chip:
        out.append(s.apply(quant.Counts.chip(features, name=name), max_fragment=args.max_fragment))
    return out


def run(args: argparse.Namespace) -> int:
    shared = _Shared(args)
    products = _products(args, shared)
    if not products:
        print(
            "ow quant: error: nothing to compute; ask for --coverage, --genes, --atac, "
            "--junctions, --labeling, --methylation, --hic, --alleles, --chip-bins or "
            "--chip-peaks",
            file=sys.stderr,
        )
        return 2
    task = quant.Task(
        args.inputs,
        args.out,
        products,
        annotation=shared.annotation,
        aliases=dict(read_aliases(args.aliases)) if args.aliases else None,
        matrix=args.matrix,
        reference=args.reference,
        variants=args.variants,
    )
    plan = quant.plan(task, **budget_of(args))
    if args.dry_run:
        print_budget(plan.budget, plan.memory, "buffers")
        for name, share in plan.stages:
            print(f"product {name}  buffer {human(share)}")
        print(f"reads   {', '.join(plan.needs)}")
        return 0
    progress = progress_for(args)
    try:
        done = plan.run(progress=progress)
    finally:
        if progress is not None:
            progress.clear()
    for file in done.files:
        print(file.path)
    if not args.quiet:
        print(
            f"quant: {done.records} records, {human(done.bytes_read)} read, "
            f"{len(products)} products, {done.seconds:.2f} s",
            file=sys.stderr,
        )
    return 0
