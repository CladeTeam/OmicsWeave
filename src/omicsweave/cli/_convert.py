"""``ow convert SOURCE TARGET``: one file into another format.

Only the command line lives here: arguments, progress and messages. The
conversion itself, formats included, is :mod:`omicsweave.convert`.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from collections import defaultdict
from typing import Any

from .. import convert
from ..formats import FORMATS, BbiOptions, File
from ._common import add_budget, budget_of, human, print_budget, progress_for

__all__ = ["register"]


def _resolutions(text: str) -> list[int]:
    try:
        return [int(r) for r in text.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a list of bin sizes: {text!r}") from None


def _routes() -> str:
    """The conversions there are, one target per line."""
    sources: dict[str, list[str]] = defaultdict(list)
    for source, target in convert.routes():
        sources[target.format].append(source.format)
    width = max(len(", ".join(s)) for s in sources.values())
    lines = [f"  {', '.join(s):<{width}}  ->  {target}" for target, s in sources.items()]
    return "conversions (source formats -> target format):\n" + "\n".join(lines)


def register(commands: Any) -> None:
    routes = convert.routes()
    p = commands.add_parser(
        "convert",
        help="convert a file into another format",
        description=(
            "Convert SOURCE into TARGET, streaming and in parallel within a memory budget.\n\n"
            "SOURCE's format is told by its content, TARGET's by its name (.fa.gz,\n"
            ".2bit, .bcf, .bam, .bb, .bw, .hic); --from and --to state them instead.\n"
            "What a conversion passes over (rows of unknown chromosomes) is printed\n"
            "as a warning.\n"
            "TARGET and its companion files (.fai, .gzi, .csi, .bai, .sidecar) appear\n"
            "only once the whole conversion succeeded."
        ),
        epilog=_routes(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("source", metavar="SOURCE", help="the file to convert")
    p.add_argument("target", metavar="TARGET", help="the file to write")
    p.add_argument(
        "--from",
        dest="source_format",
        metavar="FORMAT",
        choices=sorted({s.format for s, _ in routes}),
        help="the format of SOURCE: %(choices)s",
    )
    p.add_argument(
        "--to",
        dest="target_format",
        metavar="FORMAT",
        choices=sorted({t.format for _, t in routes}),
        help="the format of TARGET: %(choices)s",
    )

    add_budget(p, "beside TARGET")

    output = p.add_argument_group("output")
    output.add_argument(
        "-l",
        "--level",
        type=int,
        metavar="L",
        help="compression level: bgzip, BCF, BAM 0-12; bigBed, bigWig, hic 1-9",
    )
    output.add_argument(
        "--index",
        nargs="?",
        const="auto",
        choices=("auto", "bai", "csi"),
        help="also index a sorted BCF (csi) or BAM (auto: bai when every reference fits, else csi)",
    )
    output.add_argument(
        "--min-shift",
        type=int,
        metavar="N",
        help="BAM --index csi: bins of 2^N bases (default: 14)",
    )
    output.add_argument(
        "--chrom-sizes",
        metavar="FILE",
        help="bigBed / bigWig chromosome lengths (chrom.sizes, .fai or .2bit; "
        "default: each chromosome's largest interval end); pairs to hic: the "
        "chromosomes in order (default: the #chromsize lines; short text needs it)",
    )
    output.add_argument(
        "--resolutions",
        type=_resolutions,
        metavar="R,R,...",
        help=".mcool resolutions to convert (default: all); pairs to hic: bin sizes, "
        "each a multiple of the finest (default: Juicer pre's, 2.5M down to 5k)",
    )
    methyl = p.add_argument_group("bedMethyl to bigWig")
    methyl.add_argument(
        "--value",
        metavar="COLUMN",
        help="the column stored: percent (default), coverage, score, nmod, ncanonical, "
        "nother_mod, ndelete, nfail, ndiff, nnocall",
    )
    methyl.add_argument(
        "--mod",
        metavar="CODE",
        help="the items of this modification code (column 4: m, h, a, ...); "
        "required when there are several",
    )
    methyl.add_argument(
        "--strand",
        choices=("+", "-", "."),
        help="the items of this strand (default: all; overlapping strands are an error)",
    )
    output.add_argument(
        "--sidecar",
        action="store_true",
        help="also write the bigBed / bigWig / hic .sidecar (for data loaders)",
    )
    output.add_argument(
        "--block-size",
        type=int,
        metavar="N",
        help="bigBed / bigWig: children per index node (default: 256)",
    )
    output.add_argument(
        "--items-per-slot",
        type=int,
        metavar="N",
        help="bigBed / bigWig: items per data block (default: 512 bigBed, 1024 bigWig)",
    )
    output.add_argument(
        "--zoom-levels",
        type=int,
        metavar="N",
        help="bigBed / bigWig: most zoom levels derived (default: 10)",
    )
    output.add_argument(
        "--add-contigs",
        action="store_true",
        help="BCF: declare contigs the VCF header lacks (one more pass over SOURCE)",
    )
    output.add_argument(
        "--sort",
        action="store_true",
        help="BAM to BAM: sort the records by coordinate, as samtools sort does",
    )
    output.add_argument(
        "--reference",
        metavar="FASTA",
        help="CRAM to BAM: the reference the CRAM was compressed against "
        "(FASTA + .fai, bgzip FASTA + .fai/.gzi, or 2bit)",
    )
    output.add_argument(
        "--no-verify",
        dest="verify_checksums",
        action="store_false",
        help="CRAM to BAM: skip the block CRC32 and slice reference MD5 checks",
    )

    p.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="check the conversion against the budget and show the plan; write nothing",
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="print only the files written",
    )
    p.set_defaults(run=run)


def _file(path: str, name: str | None) -> File | str:
    """``path`` as a file of format ``name``; a bare path when none is given."""
    return path if name is None else FORMATS[name](path)


def _bbi(args: argparse.Namespace) -> BbiOptions | None:
    given = {
        k: v
        for k, v in (
            ("block_size", args.block_size),
            ("items_per_slot", args.items_per_slot),
            ("zoom_levels", args.zoom_levels),
        )
        if v is not None
    }
    return BbiOptions(**given) if given else None


def run(args: argparse.Namespace) -> int:
    plan = convert.plan(
        _file(args.source, args.source_format),
        _file(args.target, args.target_format),
        **budget_of(args),
        level=args.level,
        index=args.index or False,
        min_shift=args.min_shift,
        chrom_sizes=args.chrom_sizes,
        resolutions=args.resolutions,
        value=args.value,
        mod=args.mod,
        strand=args.strand,
        sidecar=args.sidecar,
        add_contigs=args.add_contigs,
        sort=args.sort,
        reference=args.reference,
        verify_checksums=args.verify_checksums,
        bbi=_bbi(args),
    )
    if args.dry_run:
        print(f"source  {plan.source.path}  ({plan.source.format})")
        print(f"target  {plan.target.path}  ({plan.target.format})")
        print_budget(plan.budget, plan.memory, "sort buffer")
        return 0
    progress = progress_for(args)
    try:
        # Warnings are printed below, after the progress line is gone.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", convert.ConvertWarning)
            done = plan.run(progress=progress)
    finally:
        if progress is not None:
            progress.clear()
    for note in done.warnings:
        print(f"warning: {note}", file=sys.stderr)
    for file in done.files:
        print(file.path)
    if not args.quiet:
        print(
            f"{done.source.format} -> {done.target.format}: {done.records} records, "
            f"{human(done.bytes_read)} read, {human(done.spilled)} spilled, "
            f"{done.seconds:.2f} s",
            file=sys.stderr,
        )
    return 0
