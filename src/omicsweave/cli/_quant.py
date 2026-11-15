"""``ow quant BAM... -o DIR``: coverage, count matrices, ATAC fragments, splice
junctions, metabolic-labeling statistics, DNA methylation, Hi-C contacts and
allele counts in one pass over coordinate-sorted BAMs.

Only the command line lives here: arguments, progress and messages. The
quantification is :mod:`omicsweave.quant`.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any

from .. import quant
from ..annotation import Annotation
from ..quant._products import Features
from ._convert import _Progress, _size, _threads, human

__all__ = ["register"]


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


def _flank(text: str) -> tuple[int, int] | None:
    if text == "none":
        return None
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


def _contexts(text: str) -> tuple[str, ...]:
    names = {c.lower(): c for c in ("CpG", "CHG", "CHH")}
    try:
        return tuple(names[x.strip().lower()] for x in text.split(","))
    except KeyError:
        raise argparse.ArgumentTypeError(f"not contexts: {text!r} (CpG,CHG,CHH)") from None


def _aliases(path: str) -> list[tuple[str, str]]:
    out = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                a, b = line.split()[:2]
                out.append((a, b))
    return out


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
            "with `ow convert in.bam sorted.bam --sort`."
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

    cov = p.add_argument_group("coverage (bigWig)")
    cov.add_argument("--coverage", action="store_true", help="write coverage.bw")
    cov.add_argument(
        "--coverage-unit",
        default="bases",
        choices=("bases", "bases+deletions", "span", "five-prime", "three-prime", "fragment"),
        help="what a read adds (default: %(default)s, as samtools depth)",
    )
    cov.add_argument(
        "--max-fragment", type=int, default=0, metavar="N", help="longest fragment counted (0: any)"
    )
    cov.add_argument(
        "--strands",
        choices=("merged", "split"),
        default="merged",
        help="one track, or .forward.bw and .reverse.bw",
    )
    cov.add_argument(
        "--scale", type=_scale, default="none", help="none, cpm, or a factor (default: none)"
    )
    cov.add_argument("--sidecar", action="store_true", help="also write each bigWig's .sidecar")
    cov.add_argument(
        "--exclude-flags",
        type=_flags,
        default=0x704,
        metavar="FLAGS",
        help="skip records with any of these FLAG bits (default: 0x704)",
    )
    cov.add_argument("--require-flags", type=_flags, default=0, metavar="FLAGS")
    cov.add_argument("--min-mapq", type=int, default=0, metavar="Q")
    cov.add_argument("--unique", action="store_true", help="only records with NH 1")

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
    genes.add_argument("--barcode-tag", default="CB", metavar="TAG", help="default: %(default)s")
    genes.add_argument("--whitelist", metavar="FILE", help="barcodes to keep, one per line")
    genes.add_argument("--biotypes", metavar="T,T,...", help="only genes of these biotypes")
    genes.add_argument("--multi-overlap", action="store_true", help="count in every feature hit")
    genes.add_argument(
        "--library",
        type=int,
        choices=(0, 1, 2),
        help="strandedness as featureCounts -s (default: 1 single-cell, 0 bulk; "
        "also the strand of --strands split)",
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
        "--labeling",
        action="store_true",
        help="labeling/conversions.tsv.gz and mismatches.tsv",
    )
    label.add_argument(
        "--conversion", default="T>C", metavar="X>Y", help="in transcript sense (default: T>C)"
    )
    label.add_argument(
        "--min-base-quality",
        type=int,
        default=20,
        metavar="Q",
        help="bases below Q are not counted (default: %(default)s)",
    )

    atac = p.add_argument_group("single-cell ATAC")
    atac.add_argument("--atac", action="store_true", help="fragments, metrics and matrices")
    atac.add_argument("--peaks", metavar="FILE.bb", help="peaks (bigBed): a peak matrix")
    atac.add_argument(
        "--tiles", type=int, default=5000, metavar="N", help="tile matrix bin size (0: none)"
    )
    atac.add_argument(
        "--gene-activity",
        type=_flank,
        default=(2000, 0),
        metavar="UP,DOWN",
        help="gene activity flanks (default: 2000,0; none to skip)",
    )
    atac.add_argument("--atac-min-mapq", type=int, default=30, metavar="Q")
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
        default="auto",
        help="XM tags, MM/ML tags, or reads against --reference / MD (default: per record)",
    )
    meth.add_argument(
        "--contexts", type=_contexts, default=("CpG",), metavar="C,C", help="default: CpG"
    )
    meth.add_argument("--merge-cpg", action="store_true", help="one line per CpG, both strands")
    meth.add_argument(
        "--methylation-min-quality", type=int, default=5, metavar="Q", help="default: %(default)s"
    )
    meth.add_argument(
        "--min-probability",
        type=float,
        default=0.5,
        metavar="P",
        help="MM/ML: call 5mC at P, unmethylated at 1-P (default: %(default)s)",
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
        "--methylation-rows",
        default="@file",
        metavar="ROWS",
        help="matrix rows: a tag (CB), @RG or @file (default: @file)",
    )
    meth.add_argument(
        "--no-cytosines",
        action="store_true",
        help="matrices only (needed for several BAMs)",
    )

    hic = p.add_argument_group("Hi-C / Micro-C (needs MC and MQ tags: samtools fixmate -m)")
    hic.add_argument("--hic", action="store_true", help="hic/contacts.pairs.gz and contacts.hic")
    hic.add_argument(
        "--hic-resolutions",
        type=_ints,
        default=None,
        metavar="R,R",
        help="bin sizes of the .hic (default: the Juicer tools'; none: no .hic)",
    )
    hic.add_argument(
        "--hic-min-mapq", type=int, default=30, metavar="Q", help="both mates (default: 30)"
    )
    hic.add_argument(
        "--min-distance", type=int, default=0, metavar="N", help="drop cis contacts closer than N"
    )
    hic.add_argument(
        "--self-ligation",
        type=int,
        default=0,
        metavar="N",
        help="drop inward/outward cis contacts closer than N",
    )
    hic.add_argument("--keep-duplicates", action="store_true", help="no contact deduplication")
    hic.add_argument("--no-pairs", action="store_true", help="skip contacts.pairs.gz")
    hic.add_argument("--genome", default="unknown", help="assembly name for pairs and .hic")
    hic.add_argument(
        "--hic-cells", action="store_true", help="per-barcode statistics by --barcode-tag"
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

    budget = p.add_argument_group("budget")
    budget.add_argument(
        "-t", "--threads", type=_threads, metavar="N", help="threads (default: every core)"
    )
    budget.add_argument(
        "-m",
        "--memory",
        type=_size,
        default="16G",
        metavar="SIZE",
        help="memory for buffers (default: %(default)s)",
    )
    budget.add_argument(
        "--tmp-dir", metavar="DIR", help="where buffers spill (default: the output directory)"
    )
    p.add_argument(
        "-n", "--dry-run", action="store_true", help="check the task against the budget only"
    )
    p.add_argument("-q", "--quiet", action="store_true", help="print only the files written")
    p.set_defaults(run=run)


def _products(args: argparse.Namespace) -> list[quant.Product]:
    out: list[quant.Product] = []
    labeled = (
        None
        if args.labeled is None
        else quant.Labeled(args.conversion, args.min_base_quality, args.labeled)
    )
    if args.coverage:
        out.append(
            quant.Coverage(
                unit=args.coverage_unit,
                max_fragment=args.max_fragment,
                strands=args.strands,
                library=args.library or 0,
                filter=quant.Filter(
                    args.require_flags, args.exclude_flags, args.min_mapq, args.unique
                ),
                scale=args.scale,
                sidecar=args.sidecar,
            )
        )
    if args.genes:
        out.append(
            quant.Genes(
                args.genes,
                span=args.gene_span,
                pairs=args.pairs,
                library=args.library,
                barcode_tag=args.barcode_tag,
                whitelist=args.whitelist,
                biotypes=args.biotypes.split(",") if args.biotypes else None,
                multi_overlap=args.multi_overlap,
                velocity=args.velocity,
                labeled=labeled,
            )
        )
    elif args.velocity or labeled is not None:
        raise SystemExit("ow quant: error: --velocity and --labeled are layers of --genes")
    if args.junctions or args.junction_matrix:
        matrix = None
        if args.junction_matrix == "single-cell":
            matrix = quant.JunctionMatrix(rows=args.barcode_tag, whitelist=args.whitelist)
        elif args.junction_matrix == "bulk":
            matrix = quant.JunctionMatrix.bulk()
        out.append(quant.Junctions(pairs=not args.single_mates, matrix=matrix))
    if args.labeling:
        out.append(
            quant.Labeling(
                library=1 if args.library is None else args.library,
                conversion=args.conversion,
                min_quality=args.min_base_quality,
            )
        )
    if args.atac:
        base = quant.Atac()
        out.append(
            quant.Atac(
                barcode_tag=args.barcode_tag,
                whitelist=args.whitelist,
                filter=quant.Filter(
                    base.filter.require, base.filter.exclude, args.atac_min_mapq, False
                ),
                max_fragment=args.max_fragment,
                fragments=not args.no_fragments,
                peaks=args.peaks,
                tiles=args.tiles or None,
                gene_activity=args.gene_activity,
            )
        )
    if args.methylation:
        matrices = []
        rows = args.methylation_rows
        specs: list[tuple[str, Features, dict[str, Any] | None]] = [
            ("bins", "bins", {"bin_size": args.methylation_bins}),
            ("genes", "gene-full", {} if args.methylation_genes else None),
            (
                "promoters",
                "promoters",
                None
                if args.methylation_promoters is None
                else dict(
                    upstream=args.methylation_promoters[0],
                    downstream=args.methylation_promoters[1],
                ),
            ),
            ("peaks", "peaks", {"peaks": args.methylation_peaks}),
        ]
        for name, features, extra in specs:
            if extra is None or (extra and next(iter(extra.values())) is None):
                continue
            matrices.append(
                quant.MethylMatrix(name, features, rows=rows, contexts=args.contexts, **extra)
            )
        out.append(
            quant.Methylation(
                source=args.methylation_source,
                contexts=args.contexts,
                merge_cpg=args.merge_cpg,
                min_quality=args.methylation_min_quality,
                min_probability=args.min_probability,
                overlap=not args.keep_overlaps,
                table=not args.no_cytosines,
                bigwig=not args.no_cytosines,
                matrices=matrices,
            )
        )
    if args.hic:
        base = quant.Contacts()
        out.append(
            quant.Contacts(
                filter=quant.Filter(
                    base.filter.require, base.filter.exclude, args.hic_min_mapq, False
                ),
                min_distance=args.min_distance,
                self_ligation=args.self_ligation,
                dedup=not args.keep_duplicates,
                pairs=not args.no_pairs,
                resolutions=base.resolutions
                if args.hic_resolutions is None
                else args.hic_resolutions,
                genome=args.genome,
                barcode_tag=args.barcode_tag if args.hic_cells else None,
                whitelist=args.whitelist if args.hic_cells else None,
            )
        )
    if args.alleles:
        genes = None
        if args.ase is not None:
            genes = quant.AlleleGenes(
                sample=args.ase or None,
                library=(1 if args.allele_mode == "single-cell" else 0)
                if args.library is None
                else args.library,
            )
        if args.allele_mode == "single-cell":
            out.append(
                quant.Alleles(
                    args.alleles, rows=args.barcode_tag, whitelist=args.whitelist, genes=genes
                )
            )
        else:
            out.append(quant.Alleles.bulk(args.alleles, genes=genes))
    chip: list[tuple[str, Features, dict[str, Any]]] = [
        ("chip_bins", "bins", {"bin_size": args.chip_bins}),
        ("chip_peaks", "peaks", {"peaks": args.chip_peaks}),
    ]
    for name, features, extra in chip:
        if next(iter(extra.values())) is not None:
            out.append(
                quant.Counts(
                    name=name,
                    features=features,
                    unit="fragment",
                    rows="@file",
                    multi_overlap=True,
                    max_fragment=args.max_fragment,
                    table=True,
                    **extra,
                )
            )
    return out


def run(args: argparse.Namespace) -> int:
    products = _products(args)
    if not products:
        print(
            "ow quant: error: nothing to compute; ask for --coverage, --genes, --atac, "
            "--junctions, --labeling, --methylation, --hic, --alleles, --chip-bins or "
            "--chip-peaks",
            file=sys.stderr,
        )
        return 2
    annotation = Annotation.read(args.gtf, attributes=None) if args.gtf else None
    task = quant.Task(
        inputs=args.inputs,
        out_dir=args.out,
        products=products,
        annotation=annotation,
        aliases=_aliases(args.aliases) if args.aliases else None,
        matrix=args.matrix,
        reference=args.reference,
        variants=args.variants,
    )
    threads = args.threads or os.cpu_count() or 1
    if args.dry_run:
        plan = quant.plan(task, threads=threads, memory=args.memory)
        m = plan.memory
        print(
            f"memory  fixed {human(m.fixed)}, pipeline {human(m.pipeline)}, "
            f"buffers {human(m.buffer)} (budget {human(plan.memory_budget)}, "
            f"{plan.threads} threads)"
        )
        for name, share in plan.stages:
            print(f"product {name}  buffer {human(share)}")
        print(f"reads   {', '.join(plan.needs)}")
        return 0
    progress = None if args.quiet or not sys.stderr.isatty() else _Progress(sys.stderr)
    try:
        done = quant.run(
            task, threads=threads, memory=args.memory, tmp_dir=args.tmp_dir, progress=progress
        )
    finally:
        if progress is not None:
            progress.clear()
    for _, path in done.files:
        print(path)
    if not args.quiet:
        print(
            f"quant: {done.records} records, {human(done.bytes_read)} read, "
            f"{len(products)} products, {done.seconds:.2f} s",
            file=sys.stderr,
        )
    return 0
