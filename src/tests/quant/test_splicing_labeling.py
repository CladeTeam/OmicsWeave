import gzip
import json
import random

import pytest

import omicsweave as ow
from omicsweave import quant
from omicsweave.annotation import Annotation
from omicsweave.cli.app import main
from omicsweave.errors import OmicsWeaveError

from ._support import CHROMS, GTF, _layer
from ._support import data as data


def test_splicing_and_labeling(tmp_path, capsys):
    """Gene A (+): exons [1000, 1500) and [2500, 3000). Reads spliced across
    its intron, inside an exon, or inside the intron; some with T>C."""
    rng = random.Random(7)
    genome = {c: [rng.choice("ACGT") for _ in range(n)] for c, n in CHROMS}
    genome["chr1"][1500:1502] = "GT"
    genome["chr1"][2498:2500] = "AG"
    fa = tmp_path / "genome.fa"
    fa.write_text("".join(f">{c}\n{''.join(s)}\n" for c, s in genome.items()))
    ref = tmp_path / "genome.fa.gz"
    ow.convert.run(fa, ref, threads=1)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)

    kinds = {"spliced": (1470, "30M1000N20M"), "exonic": (1100, "50M"), "intronic": (1800, "50M")}
    lines = ["@HD\tVN:1.6\tSO:coordinate\n"] + [f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in CHROMS]
    recs, want = [], {"hist": {}, "velocity": {}, "labeled": {}}
    for i in range(300):
        kind = rng.choice(sorted(kinds))
        pos, cigar = kinds[kind]
        pos += rng.randrange(0, 20) if kind != "spliced" else 0
        span = [pos + k for k in range(30)] + [2500 + k for k in range(20)]
        if kind != "spliced":
            span = list(range(pos, pos + 50))
        labeled = rng.random() < 0.5
        seq, cov, conv = [], 0, 0
        for p in span:
            b = genome["chr1"][p]
            if b == "T":
                cov += 1
                if labeled and rng.random() < 0.2:
                    b, conv = "C", conv + 1
            seq.append(b)
        cb = f"BC{rng.randrange(3)}-1"
        recs.append(
            (
                pos,
                f"r{i}\t0\tchr1\t{pos + 1}\t60\t{cigar}\t*\t0\t0\t{''.join(seq)}\t*"
                f"\tCB:Z:{cb}\tUB:Z:ACGT{i:05d}\tNH:i:1\n",
            )
        )
        want["hist"][(cov, conv)] = want["hist"].get((cov, conv), 0) + 1
        layer = "unspliced" if kind == "intronic" else "spliced"
        want["velocity"][(layer, cb)] = want["velocity"].get((layer, cb), 0) + 1
        if kind != "intronic" and conv >= 1:
            want["labeled"][(cb, "A")] = want["labeled"].get((cb, "A"), 0) + 1
    recs.sort(key=lambda r: r[0])
    sam = tmp_path / "sl.sam"
    sam.write_text("".join(lines + [r for _, r in recs]))
    bam = tmp_path / "sl.bam"
    ow.convert.run(sam, bam, threads=1)

    out = tmp_path / "sl"
    args = ["quant", str(bam), "-o", str(out), "--gtf", str(gtf), "--reference", str(ref)]
    args += ["--genes", "single-cell", "--velocity", "--labeled", "1", "--matrix", "mtx"]
    args += ["--junctions", "--junction-matrix", "single-cell", "--labeling", "-t", "2"]
    assert main(args) == 0
    capsys.readouterr()

    tab = (out / "junctions/junctions.tab").read_text().splitlines()
    spliced_reads = sum(1 for _, r in recs if "1000N" in r)
    assert tab == [f"chr1\t1501\t2500\t1\t1\t1\t{spliced_reads}\t0\t20"]
    sj = _layer(out / "junctions", "matrix")
    assert sum(sj.values()) == spliced_reads
    assert set(f for _, f in sj) == {"chr1:1500-2500"}

    genes = out / "genes"
    for layer in ("spliced", "unspliced"):
        got = {cb: v for (cb, g), v in _layer(genes, layer).items() if g == "A"}
        assert got == {cb: v for (k, cb), v in want["velocity"].items() if k == layer}, layer
    assert _layer(genes, "ambiguous") == {}
    assert _layer(genes, "labeled") == want["labeled"]

    text = gzip.decompress((out / "labeling/conversions.tsv.gz").read_bytes()).decode()
    rows = [line.split("\t") for line in text.splitlines()[1:]]
    assert {(int(a), int(b)): int(n) for a, b, n in rows} == want["hist"]
    summary = json.loads((out / "summary.json").read_text())
    kinds_seen = [p["product"] for p in summary["products"]]
    assert kinds_seen == ["counts", "junctions", "labeling"]
    assert summary["products"][2]["mismatch_rates"]["T>C"] > 0

    task = quant.Task(
        [bam],
        tmp_path / "p",
        [quant.Counts.rna_single_cell(velocity=True, labeled=quant.Labeled()), quant.Labeling()],
        annotation=Annotation.read(gtf),
        reference=ref,
    )
    assert "reference" in quant.plan(task, threads=1, memory="256M").needs
    with pytest.raises(OmicsWeaveError, match="velocity and labeled layers"):
        quant.plan(
            quant.Task(
                [bam], tmp_path / "p", [quant.Counts(quant.Bins(100), name="x", velocity=True)]
            )
        )
