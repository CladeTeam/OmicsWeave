"""omicsweave.quant and ``ow quant``: products from one pass over a small
sorted BAM, checked against counts computed here."""

import gzip
import json
import pickle
import random

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import quant
from omicsweave.annotation import Annotation
from omicsweave.cli.app import main
from omicsweave.errors import OmicsWeaveError
from omicsweave.query import TrackReader

CHROMS = [("chr1", 20_000), ("chr2", 8_000)]
GTF = (
    'chr1\tt\tgene\t1001\t3000\t.\t+\t.\tgene_id "A"; gene_name "a";\n'
    'chr1\tt\texon\t1001\t1500\t.\t+\t.\tgene_id "A"; transcript_id "A1";\n'
    'chr1\tt\texon\t2501\t3000\t.\t+\t.\tgene_id "A"; transcript_id "A1";\n'
    'chr1\tt\tgene\t8001\t9000\t.\t-\t.\tgene_id "B"; gene_name "b";\n'
    'chr1\tt\texon\t8001\t9000\t.\t-\t.\tgene_id "B"; transcript_id "B1";\n'
    'chr2\tt\tgene\t101\t2000\t.\t+\t.\tgene_id "C"; gene_name "c";\n'
    'chr2\tt\texon\t101\t2000\t.\t+\t.\tgene_id "C"; transcript_id "C1";\n'
)
EXONS = {
    "A": ("chr1", [(1000, 1500), (2500, 3000)], "+"),
    "B": ("chr1", [(8000, 9000)], "-"),
    "C": ("chr2", [(100, 2000)], "+"),
}


def reads(seed=0, n=800):
    """Paired reads with CB / UB tags: (chrom, pos, flag, mate pos, tlen,
    cb, ub), sorted."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        chrom, length = rng.choice(CHROMS)
        s = rng.randrange(0, length - 1000)
        frag = rng.randrange(150, 500)
        p2 = s + frag - 50
        cb = f"BC{rng.randrange(6)}-1"
        ub = "".join(rng.choice("ACGT") for _ in range(4))
        r1_left = rng.random() < 0.5
        f1, f2 = (0x40, 0x80) if r1_left else (0x80, 0x40)
        out.append((chrom, s, 0x1 | 0x2 | 0x20 | f1, p2, frag, cb, ub, f"r{i}"))
        out.append((chrom, p2, 0x1 | 0x2 | 0x10 | f2, s, -frag, cb, ub, f"r{i}"))
    out.sort(key=lambda r: ([c for c, _ in CHROMS].index(r[0]), r[1]))
    return out


def write_bam(tmp_path, recs, name="in", header="@HD\tVN:1.6\tSO:coordinate\n"):
    lines = [header] + [f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in CHROMS]
    for chrom, pos, flag, mpos, tlen, cb, ub, qname in recs:
        lines.append(
            f"{qname}\t{flag}\t{chrom}\t{pos + 1}\t60\t50M\t=\t{mpos + 1}\t{tlen}\t{'A' * 50}\t*"
            f"\tCB:Z:{cb}\tUB:Z:{ub}\tNH:i:1\n"
        )
    sam = tmp_path / f"{name}.sam"
    sam.write_text("".join(lines))
    bam = tmp_path / f"{name}.bam"
    ow.convert.run(sam, bam, threads=1)
    return bam


@pytest.fixture
def data(tmp_path):
    recs = reads()
    bam = write_bam(tmp_path, recs)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)
    return tmp_path, recs, bam, Annotation.read(gtf)


def read_10x(d):
    rows = gzip.decompress((d / "barcodes.tsv.gz").read_bytes()).decode().split()
    feats = [
        line.split("\t")[0]
        for line in gzip.decompress((d / "features.tsv.gz").read_bytes()).decode().splitlines()
    ]
    text = gzip.decompress((d / "matrix.mtx.gz").read_bytes()).decode().splitlines()
    body = [line for line in text if not line.startswith("%")]
    n_f, n_r, nnz = map(int, body[0].split())
    assert (n_f, n_r, nnz) == (len(feats), len(rows), len(body) - 1)
    return {
        (rows[int(r) - 1], feats[int(f) - 1]): int(v)
        for f, r, v in (line.split() for line in body[1:])
    }


def test_coverage_and_single_cell_genes(data):
    tmp, recs, bam, ann = data
    out = tmp / "out"
    task = quant.Task(
        [bam], out, [quant.Coverage(), quant.Counts.rna_single_cell()], annotation=ann, matrix="mtx"
    )
    p = quant.plan(task, threads=2, memory="256M")
    assert [s[0] for s in p.stages] == ["coverage", "genes"]
    assert p.memory.total == p.budget.memory == 256 << 20
    assert p.budget.threads == 2
    assert [i.path for i in p.inputs] == [bam]
    report = p.run()
    assert report.records == len(recs)
    kinds = {k for k, _ in report.files}
    assert {"bigwig", "mtx", "barcodes", "features", "json"} <= kinds
    assert json.loads((out / "summary.json").read_text()) == report.summary

    depth = {c: np.zeros(n, np.float32) for c, n in CHROMS}
    for chrom, pos, *_ in recs:
        depth[chrom][pos : pos + 50] += 1
    bw = TrackReader(out / "coverage.bw")
    for c, n in CHROMS:
        got = np.nan_to_num(bw.signal((c, 0, n)).values)
        np.testing.assert_array_equal(got, depth[c])

    want = {}
    seen = set()
    for chrom, pos, flag, _, _, cb, ub, _ in recs:
        # Forward-stranded: the fragment strand is read 1's.
        rev = bool(flag & 0x10) ^ bool(flag & 0x80)
        hits = [
            g
            for g, (c, ex, strand) in EXONS.items()
            if c == chrom
            and (strand == "-") == rev
            and any(a < pos + 50 and pos < b for a, b in ex)
        ]
        if len(hits) == 1 and (cb, hits[0], ub) not in seen:
            seen.add((cb, hits[0], ub))
            want[(cb, hits[0])] = want.get((cb, hits[0]), 0) + 1
    assert read_10x(out / "genes") == want


def test_atac_and_cli(data, capsys):
    tmp, recs, bam, ann = data
    out = tmp / "cli"
    gtf = tmp / "genes.gtf"
    args = ["quant", str(bam), "-o", str(out), "--gtf", str(gtf), "--atac", "--tiles", "1000"]
    args += ["--coverage", "--strands", "split", "--library", "2", "-t", "3", "-m", "256M"]
    assert main(args + ["-n"]) == 0
    assert "buffer" in capsys.readouterr().out
    assert not out.exists()
    assert main(args) == 0
    printed = capsys.readouterr().out.split()
    for f in [
        "coverage.forward.bw",
        "coverage.reverse.bw",
        "atac/fragments.tsv.gz",
        "atac/fragments.tsv.gz.tbi",
        "atac/barcode_metrics.csv",
        "atac/tiles/matrix.cpk",
        "atac/gene_activity/matrix.cpk",
        "summary.json",
    ]:
        assert str(out / f) in printed, f
    frags = {}
    for chrom, pos, _, _, tlen, cb, _, _ in recs:
        if tlen > 0:
            key = (chrom, pos + 4, pos + tlen - 5, cb)
            frags[key] = frags.get(key, 0) + 1
    lines = gzip.decompress((out / "atac/fragments.tsv.gz").read_bytes()).decode().splitlines()
    got = {(c, int(s), int(e), b): int(n) for c, s, e, b, n in (x.split("\t") for x in lines)}
    assert got == frags


def test_refusals(tmp_path, data):
    _, recs, _, _ = data
    swapped = list(recs)
    swapped[10], swapped[500] = swapped[500], swapped[10]
    bam = write_bam(tmp_path, swapped, "swapped")
    task = quant.Task([bam], tmp_path / "o", [quant.Coverage()])
    with pytest.raises(OmicsWeaveError, match="not sorted by coordinate"):
        quant.run(task, threads=1)
    unsorted = write_bam(tmp_path, recs, "unsorted", "@HD\tVN:1.6\tSO:unsorted\n")
    with pytest.raises(OmicsWeaveError, match="--sort"):
        quant.plan(quant.Task([unsorted], tmp_path / "o", [quant.Coverage()]))
    with pytest.raises(OmicsWeaveError, match="needs at least"):
        quant.plan(task, threads=4, memory="4M")
    with pytest.raises(OmicsWeaveError, match="annotation"):
        quant.plan(quant.Task([bam], tmp_path / "o", [quant.Counts.rna_single_cell()]))
    assert main(["quant", str(bam), "-o", str(tmp_path / "x")]) == 2


def _layer(d, layer):
    rows = gzip.decompress((d / "barcodes.tsv.gz").read_bytes()).decode().split()
    feats = [
        line.split("\t")[0]
        for line in gzip.decompress((d / "features.tsv.gz").read_bytes()).decode().splitlines()
    ]
    text = gzip.decompress((d / f"{layer}.mtx.gz").read_bytes()).decode().splitlines()
    body = [line for line in text if not line.startswith("%")]
    return {
        (rows[int(r) - 1], feats[int(f) - 1]): int(v)
        for f, r, v in (line.split() for line in body[1:])
    }


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


def _genome(rng):
    return {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CHROMS}


def _fasta(tmp_path, genome):
    fa = tmp_path / "genome.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    ref = tmp_path / "genome.fa.gz"
    ow.convert.run(fa, ref, threads=1)
    return ref


def _sam_to_bam(tmp_path, name, lines):
    head = ["@HD\tVN:1.6\tSO:coordinate\n"] + [f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in CHROMS]
    sam = tmp_path / f"{name}.sam"
    sam.write_text("".join(head + lines))
    bam = tmp_path / f"{name}.bam"
    ow.convert.run(sam, bam, threads=1)
    return bam


def test_methylation_and_cli(tmp_path, capsys):
    """Single-end top-strand bisulfite reads: a reference C stays C when its
    position is a multiple of 3 (methylated), else reads T."""
    rng = random.Random(11)
    genome = _genome(rng)
    ref = _fasta(tmp_path, genome)
    recs, want = [], {}
    for i in range(400):
        chrom = rng.choice([c for c, _ in CHROMS])
        g = genome[chrom]
        pos = rng.randrange(0, len(g) - 60)
        seq = []
        for p in range(pos, pos + 50):
            b = g[p]
            if b == "C":
                meth = p % 3 == 0
                n1, n2 = (g[p + 1], g[p + 2]) if p + 2 < len(g) else ("N", "N")
                ctx = "CpG" if n1 == "G" else "CHG" if n2 == "G" else "CHH"
                m, u = want.get((chrom, p, ctx), (0, 0))
                want[(chrom, p, ctx)] = (m + meth, u + (not meth))
                b = "C" if meth else "T"
            seq.append(b)
        cb = f"BC{rng.randrange(3)}"
        recs.append(
            (
                ([c for c, _ in CHROMS].index(chrom), pos),
                f"r{i}\t0\t{chrom}\t{pos + 1}\t60\t50M\t*\t0\t0\t{''.join(seq)}\t*\tCB:Z:{cb}\n",
            )
        )
    recs.sort(key=lambda r: r[0])
    bam = _sam_to_bam(tmp_path, "bs", [r for _, r in recs])
    out = tmp_path / "m"
    args = ["quant", str(bam), "-o", str(out), "--reference", str(ref), "--methylation"]
    args += ["--contexts", "CpG,CHG,CHH", "--methylation-bins", "1000"]
    args += ["--methylation-rows", "CB", "--matrix", "mtx", "-t", "2"]
    assert main(args) == 0
    capsys.readouterr()
    order = {c: i for i, (c, _) in enumerate(CHROMS)}
    lines = gzip.decompress((out / "methylation/cytosines.tsv.gz").read_bytes()).decode()
    got = {}
    for line in lines.splitlines():
        c, s, e, ctx, strand, m, u = line.split("\t")
        assert (int(e) - int(s), strand) == (1, "+")
        got[(c, int(s), ctx)] = (int(m), int(u))
    assert got == want
    keys = [(order[c], s) for c, s, _ in got]
    assert keys == sorted(keys)
    assert (out / "methylation/cytosines.tsv.gz.tbi").exists()
    total = _layer(out / "methylation/bins", "total")
    assert sum(total.values()) == sum(m + u for m, u in want.values())
    meth = _layer(out / "methylation/bins", "methylated")
    assert sum(meth.values()) == sum(m for m, _ in want.values())
    bw = TrackReader(out / "methylation/coverage.bw")
    c, s, _ = next(iter(want))
    cov = np.nan_to_num(bw.signal((c, s, s + 1)).values)
    assert cov[0] == sum(m + u for (cc, ss, _), (m, u) in want.items() if (cc, ss) == (c, s))
    summary = json.loads((out / "summary.json").read_text())
    p = summary["products"][0]
    assert p["product"] == "methylation" and p["reads_called"] == 400
    assert p["mbias"]["CpG"]["read1"]["methylated"]
    with pytest.raises(OmicsWeaveError, match="one BAM at a time"):
        quant.plan(quant.Task([bam, bam], tmp_path / "x", [quant.Methylation()]))
    assert (
        quant.plan(
            quant.Task([bam, bam], tmp_path / "x", [quant.Methylation(table=False, bigwig=False)])
        ).stages[0][0]
        == "methylation"
    )


def test_hic_and_alleles(tmp_path):
    rng = random.Random(12)
    genome = _genome(rng)
    recs, contacts = [], []
    names = [c for c, _ in CHROMS]
    for i in range(300):
        a = (rng.randrange(2), rng.randrange(0, 7000), rng.random() < 0.5)
        b = (rng.randrange(2), rng.randrange(0, 7000), rng.random() < 0.5)
        cb = f"BC{rng.randrange(3)}"
        for k, (me, mate) in enumerate([(a, b), (b, a)]):
            flag = (
                0x1 | (0x40 if k == 0 else 0x80) | (0x10 if me[2] else 0) | (0x20 if mate[2] else 0)
            )
            rnext = "=" if me[0] == mate[0] else names[mate[0]]
            recs.append(
                (
                    (me[0], me[1]),
                    f"p{i}\t{flag}\t{names[me[0]]}\t{me[1] + 1}\t60\t40M\t{rnext}\t{mate[1] + 1}\t0"
                    f"\t{'A' * 40}\t*\tMC:Z:40M\tMQ:i:60\tCB:Z:{cb}\n",
                )
            )
        ends = sorted((e[0], e[1] + (40 if e[2] else 1), "-" if e[2] else "+") for e in (a, b))
        contacts.append((ends[0], ends[1], cb))
    recs.sort(key=lambda r: r[0])
    bam = _sam_to_bam(tmp_path, "hic", [r for _, r in recs])
    out = tmp_path / "h"
    task = quant.Task(
        [bam], out, [quant.Contacts(resolutions=[1000, 5000], barcode_tag="CB", dedup=False)]
    )
    report = quant.run(task, threads=2, memory="256M")
    p = report.summary["products"][0]
    assert p["contacts"] == len(contacts)
    assert {c["barcode"] for c in p["cells"]} == {c for *_, c in contacts}
    text = gzip.decompress((out / "hic/contacts.pairs.gz").read_bytes()).decode()
    body = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
    got = sorted(
        ((names.index(f[1]), int(f[2]), f[5]), (names.index(f[3]), int(f[4]), f[6]), f[10])
        for f in body
    )
    assert got == sorted(contacts)
    assert "#chromsize: chr1 20000" in text
    assert (out / "hic/contacts.hic").stat().st_size > 0

    # Alleles: reads over three SNVs with REF or ALT bases.
    sites = [("chr1", 1000, "A"), ("chr1", 5000, "G"), ("chr2", 300, "T")]
    vcf = tmp_path / "s.vcf"
    lines = ["##fileformat=VCFv4.2\n"] + [f"##contig=<ID={c},length={n}>\n" for c, n in CHROMS]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
    for c, p0, alt in sites:
        r = genome[c][p0]
        alt = alt if alt != r else "C" if r != "C" else "G"
        lines.append(f"{c}\t{p0 + 1}\t.\t{r}\t{alt}\t.\t.\t.\tGT\t0|1\n")
    vcf.write_text("".join(lines))
    bcf = tmp_path / "s.bcf"
    ow.convert.run(vcf, bcf, index=True, threads=1)
    alts = {}
    for line in lines[-3:]:
        f = line.split("\t")
        alts[(f[0], int(f[1]) - 1)] = (f[3], f[4])
    recs, want = [], {}
    for i in range(200):
        (c, p0), (r, a) = rng.choice(sorted(alts.items()))
        pos = p0 - rng.randrange(0, 40)
        allele = rng.randrange(2)
        seq = list(genome[c][pos : pos + 50])
        seq[p0 - pos] = (r, a)[allele]
        cb, ub = f"BC{rng.randrange(3)}", f"U{i}"
        recs.append(
            (
                (names.index(c), pos),
                f"a{i}\t0\t{c}\t{pos + 1}\t60\t50M\t*\t0\t0\t{''.join(seq)}\t*"
                f"\tCB:Z:{cb}\tUB:Z:{ub}\n",
            )
        )
        key = (cb, f"{c}:{p0 + 1}:{r}>{a}")
        layer = ("ref", "alt")[allele]
        want.setdefault(layer, {})
        want[layer][key] = want[layer].get(key, 0) + 1
    recs.sort(key=lambda r: r[0])
    bam = _sam_to_bam(tmp_path, "al", [r for _, r in recs])
    out = tmp_path / "a"
    quant.run(quant.Task([bam], out, [quant.Alleles(bcf)], matrix="mtx"), threads=2, memory="256M")
    for layer in ("ref", "alt"):
        assert _layer(out / "alleles", layer) == want[layer], layer
    table = gzip.decompress((out / "alleles/variants.tsv.gz").read_bytes()).decode().splitlines()
    assert len(table) == 4 and table[0].startswith("chrom\tpos\tref\talt")


def test_regions_count_as_their_own_features(data):
    tmp, _, bam, _ = data
    size = 2_000
    chroms = [c for c, _ in CHROMS]
    seqid, start, end = [], [], []
    for i, (_, n) in enumerate(CHROMS):
        for s in range(0, n, size):
            seqid.append(i)
            start.append(s)
            end.append(min(s + size, n))
    tiles = ow.annotation.IntervalSet.from_arrays(
        chroms, np.array(seqid), np.array(start), np.array(end)
    )
    task = quant.Task(
        [bam],
        tmp / "regions",
        [
            quant.Counts.chip(quant.Bins(size), name="bins"),
            quant.Counts.chip(quant.Regions(tiles), name="regions"),
        ],
        matrix="mtx",
    )
    quant.run(task, threads=2, memory="256M")
    bins = read_10x(tmp / "regions" / "bins")
    assert bins and read_10x(tmp / "regions" / "regions") == bins


def test_presets_come_from_the_native_layer():
    single, bulk = quant.Counts.rna_single_cell(), quant.Counts.rna_bulk(unit="pair")
    assert (single.rows, single.dedup, single.strand) == ("CB", "umi", "forward")
    assert single.filter.unique and isinstance(single.features, quant.GeneExons)
    assert (bulk.rows, bulk.unit, bulk.table, bulk.strand) == ("@file", "pair", True, "unstranded")
    chip = quant.Counts.chip(quant.Peaks("p.bb"))
    assert chip.features == quant.Peaks("p.bb") and chip.multi_overlap
    assert quant.Atac().gene_activity == quant.GeneActivity(2000, 0)
    assert quant.Alleles.bulk("s.bcf").unit == "fragment"
    assert quant.JunctionMatrix.bulk().rows == "@file"
    assert quant.Contacts().resolutions[-1] == 5_000


def test_products_are_checked_by_the_native_layer(data):
    tmp, _, bam, _ = data

    def plan(*products):
        return quant.plan(quant.Task([bam], tmp / "x", products), memory="256M")

    with pytest.raises(ValueError, match="split strands"):
        plan(quant.Coverage(strand="forward"))
    with pytest.raises(ValueError, match="fragment unit"):
        plan(quant.Coverage(max_fragment=500))
    with pytest.raises(ValueError, match="strand must be"):
        plan(quant.Coverage(strands="split", strand="sense"))
    with pytest.raises(ValueError, match="two characters"):
        plan(quant.Counts(quant.Bins(1000), umi_tag="UMI"))
    with pytest.raises(ValueError, match="none or umi"):
        plan(quant.Junctions(matrix=quant.JunctionMatrix(dedup="position")))
    with pytest.raises(AttributeError):
        quant.plan(quant.Task([bam], tmp / "x", [object()]))


def test_a_plan_runs_without_planning_again(data):
    tmp, recs, bam, _ = data
    p = quant.plan(quant.Task([bam], tmp / "once", [quant.Coverage()]), threads=1, memory="256M")
    ((path, refs, groups, size),) = [
        (i.path, i.references, i.read_groups, i.size) for i in p.inputs
    ]
    assert path == bam and refs == tuple(CHROMS) and groups == () and size > 0
    report = p.run(progress=lambda *a: None)
    assert report.records == len(recs) and 0 <= report.scan_seconds <= report.seconds
    assert [f.kind for f in report.files if f.path.name == "coverage.bw"] == ["bigwig"]
    with pytest.raises(TypeError, match="process-local"):
        pickle.dumps(p)


def test_a_subset_of_the_annotation_picks_the_genes(data):
    tmp, _, bam, ann = data
    only_a = ann.subset(genes=["A"])
    task = quant.Task(
        [bam],
        tmp / "subset",
        [
            quant.Counts.rna_single_cell(name="all"),
            quant.Counts.rna_single_cell(name="a", subset=only_a, velocity=True),
        ],
        annotation=ann,
        matrix="mtx",
    )
    quant.run(task, threads=2, memory="256M")
    every, a = read_10x(tmp / "subset" / "all"), read_10x(tmp / "subset" / "a")
    assert a and {f for _, f in a} == {"A"}
    assert a == {k: v for k, v in every.items() if k[1] == "A"}
    other = Annotation.read(tmp / "genes.gtf").subset(genes=["C"])
    rows = type(other)(other.genes + 10, other.transcripts)
    with pytest.raises(OmicsWeaveError, match="rows of the annotation"):
        quant.plan(
            quant.Task([bam], tmp / "x", [quant.Counts.rna_bulk(subset=rows)], annotation=ann)
        )
