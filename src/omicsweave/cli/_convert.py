"""``ow convert SOURCE TARGET``: one file into another format.

Only the command line lives here: arguments, progress and messages. The
conversion itself, formats included, is :mod:`omicsweave.convert`.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from collections import defaultdict
from typing import Any, TextIO

from .. import _native, convert
from ..errors import OmicsWeaveError
from ..formats import FORMATS, File

__all__ = ["register"]

human = _native.human_size


def _size(text: str) -> str:
    try:
        size = _native.parse_size(text)
    except OmicsWeaveError:
        size = 0
    if size == 0:
        raise argparse.ArgumentTypeError(f"not a memory budget: {text!r} (e.g. 512M, 8G)")
    return text


def _threads(text: str) -> int:
    try:
        n = int(text)
    except ValueError:
        n = 0
    if not 1 <= n <= sys.maxsize:
        raise argparse.ArgumentTypeError(f"not a thread count: {text!r} (at least 1)")
    return n


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

    budget = p.add_argument_group("budget")
    budget.add_argument(
        "-t",
        "--threads",
        type=_threads,
        metavar="N",
        help="threads for the parallel stages (default: every core)",
    )
    budget.add_argument(
        "-m",
        "--memory",
        type=_size,
        default="4G",
        metavar="SIZE",
        help="memory for buffers, e.g. 512M, 8G (default: %(default)s)",
    )
    budget.add_argument(
        "--tmp-dir",
        metavar="DIR",
        help="where sorts spill when they outgrow the budget (default: beside TARGET)",
    )

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
        "--no-sidecar",
        dest="sidecar",
        action="store_false",
        help="do not write the bigBed / bigWig / hic .sidecar",
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


class _Progress:
    """One line on ``stream``, rewritten as the conversion goes."""

    def __init__(self, stream: TextIO) -> None:
        self.stream = stream
        self.shown = False

    def __call__(self, phase: str, read: int, total: int) -> None:
        share = f" {100 * read / total:5.1f}%" if total else ""
        self.stream.write(f"\r{phase:<8} {human(read)} / {human(total)}{share}\x1b[K")
        self.stream.flush()
        self.shown = True

    def clear(self) -> None:
        if self.shown:
            self.stream.write("\r\x1b[K")
            self.stream.flush()


def run(args: argparse.Namespace) -> int:
    source = _file(args.source, args.source_format)
    target = _file(args.target, args.target_format)
    threads = args.threads or os.cpu_count() or 1
    if args.dry_run:
        plan = convert.plan(
            source,
            target,
            threads=threads,
            memory=args.memory,
            chrom_sizes=args.chrom_sizes,
            resolutions=args.resolutions,
            value=args.value,
            mod=args.mod,
            strand=args.strand,
            sort=args.sort,
        )
        b, m = plan.budget, plan.memory
        print(f"source  {plan.source.path}  ({plan.source.format})")
        print(f"target  {plan.target.path}  ({plan.target.format})")
        print(
            f"memory  fixed {human(m.fixed)}, pipeline {human(m.pipeline)}, "
            f"sort buffer {human(m.buffer)} (budget {human(b.memory)}, {b.threads} threads)"
        )
        return 0
    progress = None if args.quiet or not sys.stderr.isatty() else _Progress(sys.stderr)
    try:
        # Warnings are printed below, after the progress line is gone.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", convert.ConvertWarning)
            done = convert.run(
                source,
                target,
                threads=threads,
                memory=args.memory,
                tmp_dir=args.tmp_dir,
                level=args.level,
                index=args.index or False,
                chrom_sizes=args.chrom_sizes,
                resolutions=args.resolutions,
                value=args.value,
                mod=args.mod,
                strand=args.strand,
                sidecar=args.sidecar,
                add_contigs=args.add_contigs,
                sort=args.sort,
                reference=args.reference,
                progress=progress,
            )
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
