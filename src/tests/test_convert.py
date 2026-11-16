import gzip
import struct

import numpy as np
import pytest

import omicsweave as ow
from omicsweave.formats import (
    Bai,
    Bam,
    Bcf,
    BedGraph,
    BgzipFasta,
    BigWig,
    Csi,
    Fai,
    Fasta,
    Gzi,
    GzipBedGraph,
    GzipFasta,
    Sidecar,
    TwoBit,
    Vcf,
    Wig,
    detect,
)
from omicsweave.query import IntervalReader, SequenceReader, TrackReader


def random_fasta(rng, n=3):
    records = {}
    for i in range(n):
        length = int(rng.integers(1, 50_000))
        records[f"seq{i}"] = bytes(rng.choice(list(b"ACGTNacgt"), length).astype(np.uint8))
    lines = []
    for name, seq in records.items():
        lines.append(b">" + name.encode())
        lines.extend(seq[i : i + 60] for i in range(0, len(seq), 60))
    return records, b"\n".join(lines) + b"\n"


@pytest.mark.parametrize("threads", [1, 4])
def test_fasta_to_bgzip_and_2bit(tmp_path, threads):
    rng = np.random.default_rng(threads)
    records, text = random_fasta(rng)
    src = tmp_path / "a.fa.gz"
    src.write_bytes(gzip.compress(text))

    report = ow.convert.run(src, tmp_path / "a.bgz.fa.gz", threads=threads, memory="64M")
    assert report.records == len(records)
    assert report.source == GzipFasta(src)
    assert report.target == BgzipFasta(tmp_path / "a.bgz.fa.gz")
    assert report.companions == (
        Fai(tmp_path / "a.bgz.fa.gz.fai"),
        Gzi(tmp_path / "a.bgz.fa.gz.gzi"),
    )
    with SequenceReader(tmp_path / "a.bgz.fa.gz") as f:
        for name, seq in records.items():
            assert f.fetch((name, 0, len(seq))).bases.tobytes() == seq

    ow.convert.run(src, tmp_path / "a.2bit", threads=threads, memory="64M")
    with SequenceReader(tmp_path / "a.2bit", soft_mask=True) as f:
        for name, seq in records.items():
            want = bytes(c if chr(c) in "ACGTacgt" else ord("N") for c in seq)
            assert f.fetch((name, 0, len(seq))).bases.tobytes() == want


def test_bed_to_bigbed_reads_back(tmp_path):
    rng = np.random.default_rng(3)
    rows = []
    for chrom, length in (("chr2", 90_000), ("chr1", 120_000)):
        for i in range(500):
            s = int(rng.integers(0, length - 100))
            rows.append((chrom, s, s + int(rng.integers(1, 100)), f"p{i}"))
    rng.shuffle(rows)
    (tmp_path / "a.bed").write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in rows))
    (tmp_path / "sizes").write_text("chr1\t120000\nchr2\t90000\n")

    report = ow.convert.run(
        tmp_path / "a.bed",
        tmp_path / "a.bb",
        threads=2,
        memory="32M",
        chrom_sizes=tmp_path / "sizes",
        sidecar=True,
    )
    assert report.records == len(rows)
    assert (tmp_path / "a.bb.sidecar").exists()
    with IntervalReader(tmp_path / "a.bb") as f:
        assert dict(f.chroms) == {"chr2": 90_000, "chr1": 120_000}
        got = f.fetch(("chr1", 0, 120_000))
        want = sorted((s, e, n) for c, s, e, n in rows if c == "chr1")
        assert (
            sorted(zip(got.start.tolist(), got.end.tolist(), got.annotation.to_list(), strict=True))
            == want
        )


@pytest.mark.parametrize("compressed", [False, True])
def test_bedgraph_and_wig_to_bigwig_read_back(tmp_path, compressed):
    rng = np.random.default_rng(5)
    rows = []
    for chrom, length in (("chr2", 90_000), ("chr1", 120_000)):
        ends = np.sort(rng.choice(np.arange(1, length), 800, replace=False))
        starts = np.concatenate([[0], ends[:-1]]) + rng.integers(0, 2, len(ends))
        for s, e in zip(starts.tolist(), ends.tolist(), strict=True):
            if s < e:
                rows.append((chrom, s, e, float(np.float32(rng.normal()))))
    order = rng.permutation(len(rows))
    text = "".join(f"{rows[i][0]}\t{rows[i][1]}\t{rows[i][2]}\t{rows[i][3]}\n" for i in order)
    src = tmp_path / ("a.bg.gz" if compressed else "a.bg")
    src.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
    (tmp_path / "sizes").write_text("chr1\t120000\nchr2\t90000\n")
    assert detect(src) == (GzipBedGraph(src) if compressed else BedGraph(src))

    report = ow.convert.run(
        src,
        tmp_path / "a.bw",
        threads=2,
        memory="32M",
        chrom_sizes=tmp_path / "sizes",
        sidecar=True,
    )
    assert report.records == len(rows)
    assert report.target == BigWig(tmp_path / "a.bw")
    assert report.companions == (Sidecar(tmp_path / "a.bw.sidecar"),)
    with TrackReader(report.target) as f:
        assert dict(f.chroms) == {"chr2": 90_000, "chr1": 120_000}
        for chrom in ("chr1", "chr2"):
            got = f.intervals((chrom, 0, f.chroms[chrom]))
            want = sorted((s, e, v) for c, s, e, v in rows if c == chrom)
            assert (
                list(zip(got.start.tolist(), got.end.tolist(), got.value.tolist(), strict=True))
                == want
            )
    with pytest.raises(ow.errors.InvalidError, match="convert it to bigWig"):
        ow.open(src)

    wig = tmp_path / "a.wig"
    wig.write_text(
        "track type=wiggle_0\nfixedStep chrom=chr1 start=11 step=10 span=5\n1\n2.5\n"
        "variableStep chrom=chr2 span=2\n101 -1\n51 3\nchr1 1000 1010 4\n"
    )
    assert detect(wig) == Wig(wig)
    report = ow.convert.run(wig, tmp_path / "w.bw", threads=1)
    assert report.companions == ()
    with TrackReader(report.target) as f:
        assert dict(f.chroms) == {"chr1": 1010, "chr2": 102}
        values = f.signal(("chr1", 8, 26)).values
        assert np.array_equal(
            values,
            np.array([np.nan] * 2 + [1] * 5 + [np.nan] * 5 + [2.5] * 5 + [np.nan], "f4"),
            equal_nan=True,
        )
        got = f.intervals(("chr2", 0, 102))
        assert got.start.tolist() == [50, 100] and got.value.tolist() == [3, -1]
    (tmp_path / "bad.wig").write_text("variableStep chrom=chr1 span=5\n1 1\n3 2\n")
    with pytest.raises(ow.errors.FormatError, match="line 3"):
        ow.convert.run(tmp_path / "bad.wig", tmp_path / "bad.bw")
    assert not (tmp_path / "bad.bw").exists()


def test_plan_and_errors(tmp_path):
    (tmp_path / "a.fa").write_bytes(b">a\nACGT\n")
    plan = ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.2bit", threads=2, memory="1G")
    assert plan.source == Fasta(tmp_path / "a.fa") and plan.target == TwoBit(tmp_path / "a.2bit")
    assert plan.budget.threads == 2 and plan.budget.memory == 1 << 30
    assert plan.memory.fixed + plan.memory.pipeline <= 1 << 30
    with pytest.raises(ow.errors.InvalidError):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.2bit", threads=2, memory=1000)
    with pytest.raises(ow.errors.InvalidError):
        ow.convert.run(tmp_path / "a.fa", tmp_path / "a.xyz")
    with pytest.raises(ow.errors.InvalidError):
        ow.convert.run(tmp_path / "a.fa", tmp_path / "a.2bit", index=True)
    assert not (tmp_path / "a.2bit").exists()


def test_progress_is_reported(tmp_path):
    rng = np.random.default_rng(9)
    _, text = random_fasta(rng, n=2)
    (tmp_path / "a.fa").write_bytes(text)
    seen = []
    ow.convert.run(
        tmp_path / "a.fa", tmp_path / "a.fa.gz", threads=2, progress=lambda *p: seen.append(p)
    )
    assert all(isinstance(phase, str) and read <= total for phase, read, total in seen)


def test_formats_given_by_file_objects(tmp_path):
    (tmp_path / "a.txt").write_bytes(b">a\nACGTN\n")
    with pytest.raises(ow.errors.InvalidError, match="from its name"):
        ow.convert.run(tmp_path / "a.txt", tmp_path / "out.bin")
    report = ow.convert.run(tmp_path / "a.txt", TwoBit(tmp_path / "out.bin"), threads=1)
    assert report.source == Fasta(tmp_path / "a.txt") and report.companions == ()
    with SequenceReader(report.target) as f:
        assert f.fetch(("a", 0, 5)).text == "ACGTN"
    with pytest.raises(ow.errors.InvalidError, match="no conversion from 2bit to bcf"):
        ow.convert.run(report.target, tmp_path / "x.bcf")
    pairs = {(s.format, t.format) for s, t in ow.convert.routes()}
    assert ("gzip-fasta", "2bit") in pairs and ("mcool", "hic") in pairs


def test_options_and_companions_come_typed(tmp_path):
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c,length=100>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "c\t5\t.\tA\tG\t50\tPASS\t.\n"
    )
    report = ow.convert.run(vcf, tmp_path / "v.bcf", index=True, threads=1)
    assert report.source == Vcf(vcf) and report.target == Bcf(tmp_path / "v.bcf")
    assert report.companions == (Csi(tmp_path / "v.bcf.csi"),)
    with pytest.raises(ow.errors.InvalidError, match="does not apply to VCF to BCF"):
        ow.convert.run(vcf, tmp_path / "w.bcf", chrom_sizes=vcf)
    with pytest.raises(ow.errors.InvalidError, match="not .bai"):
        ow.convert.run(vcf, tmp_path / "w.bcf", index="bai")
    assert not (tmp_path / "w.bcf").exists()


def unsorted_sam(rng, n=600):
    """SAM text of `n` reads in no order over two references, some unmapped."""
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:a\tLN:5000", "@SQ\tSN:b\tLN:5000"]
    for i in range(n):
        ref = rng.choice(["a", "b", "*"], p=[0.45, 0.45, 0.1])
        if ref == "*":
            lines.append(f"u{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*")
        else:
            pos = int(rng.integers(1, 4000))
            flag = int(rng.choice([0, 16]))
            lines.append(f"r{i}\t{flag}\t{ref}\t{pos}\t30\t4M\t*\t0\t0\tACGT\t*")
    return "\n".join(lines) + "\n"


def bam_records(path):
    """The header text and the (tid, pos, name) of every record of a BAM."""
    raw = gzip.decompress(path.read_bytes())
    assert raw[:4] == b"BAM\x01"
    (l_text,) = struct.unpack_from("<I", raw, 4)
    at = 8 + l_text
    (n_ref,) = struct.unpack_from("<I", raw, at)
    at += 4
    for _ in range(n_ref):
        (l_name,) = struct.unpack_from("<I", raw, at)
        at += 8 + l_name
    records = []
    while at < len(raw):
        size, tid, pos, l_name = struct.unpack_from("<IiiB", raw, at)
        records.append((tid, pos, raw[at + 36 : at + 35 + l_name].decode()))
        at += 4 + size
    return raw[8 : 8 + l_text], records


def test_bam_sorts_by_coordinate(tmp_path):
    rng = np.random.default_rng(5)
    (tmp_path / "reads.sam").write_text(unsorted_sam(rng))
    unsorted, out = tmp_path / "unsorted.bam", tmp_path / "sorted.bam"
    ow.convert.run(tmp_path / "reads.sam", unsorted, threads=1)
    with pytest.raises(ow.errors.InvalidError, match="pass sort"):
        ow.convert.run(unsorted, out)
    with pytest.raises(ow.errors.InvalidError, match="index needs a coordinate-sorted"):
        ow.convert.run(unsorted, out, index=True)
    with pytest.raises(ow.errors.InvalidError, match="pass sort"):
        ow.convert.plan(unsorted, out)
    with pytest.raises(ow.errors.InvalidError, match="does not apply to SAM to BAM"):
        ow.convert.run(tmp_path / "reads.sam", out, sort=True)
    assert not out.exists()
    plan = ow.convert.plan(unsorted, out, sort=True, threads=2, memory="64M")
    assert plan.target == Bam(out) and plan.memory.buffer > 0
    report = ow.convert.run(unsorted, out, sort=True, index=True, threads=2, memory="64M")
    assert report.source == Bam(unsorted) and report.target == Bam(out)
    assert report.companions == (Bai(tmp_path / "sorted.bam.bai"),)
    assert report.records == 600
    text, records = bam_records(out)
    assert text.startswith(b"@HD\tVN:1.6\tSO:coordinate\n")
    keys = [(tid % (1 << 32), pos) for tid, pos, _ in records]
    assert keys == sorted(keys) and keys[-1][0] == (1 << 32) - 1
    assert sorted(r[2] for r in records) == sorted(r[2] for r in bam_records(unsorted)[1])


def test_pairs_to_hic_reads_back_and_warns(tmp_path):
    from omicsweave.formats import GzipPairs, Hic, JuicerShort, Pairs
    from omicsweave.query import ContactReader
    from omicsweave.types import ContactWindow

    rng = np.random.default_rng(9)
    sizes = {"chr1": 300_000, "chr2": 200_000}
    rows = []
    for i in range(3_000):
        c1, c2 = rng.choice(["chr1", "chr2", "chrM"], 2, p=[0.6, 0.38, 0.02])
        p1 = int(rng.integers(1, sizes.get(c1, 16_000)))
        p2 = int(rng.integers(1, sizes.get(c2, 16_000)))
        rows.append((f"r{i}", c1, p1, c2, p2))
    head = "## pairs format v1.0\n" + "".join(f"#chromsize: {c} {n}\n" for c, n in sizes.items())
    head += "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    text = head + "".join(f"{r}\t{a}\t{p}\t{b}\t{q}\t+\t-\n" for r, a, p, b, q in rows)
    src = tmp_path / "a.pairs.gz"
    src.write_bytes(gzip.compress(text.encode()))
    assert detect(src) == GzipPairs(src)

    with pytest.warns(ow.convert.ConvertWarning, match=r"chrM \(\d+\)"):
        report = ow.convert.run(src, tmp_path / "a.hic", resolutions=[10_000, 50_000], threads=2)
    known = [r for r in rows if r[1] != "chrM" and r[3] != "chrM"]
    assert report.records == len(known) and len(report.warnings) == 1
    assert report.target == Hic(tmp_path / "a.hic")
    want: dict[tuple[int, int], float] = {}
    for _, a, p, b, q in known:
        if a == b == "chr1":
            x, y = sorted((p // 10_000, q // 10_000))
            want[(x * 10_000, y * 10_000)] = want.get((x * 10_000, y * 10_000), 0) + 1
    with ContactReader(report.target) as f:
        region = ("chr1", 0, sizes["chr1"])
        got = f.pixels(ContactWindow(region, region, 10_000))
        pixels = {
            (r, c): v
            for r, c, v in zip(got.row.tolist(), got.col.tolist(), got.value.tolist(), strict=True)
            if r <= c
        }
    assert pixels == want

    short = tmp_path / "a.txt"
    short.write_text("".join(f"0 {a} {p} 0 16 {b} {q} 1\n" for _, a, p, b, q in known))
    assert detect(short) == JuicerShort(short)
    with pytest.raises(ow.errors.InvalidError, match="chromosome sizes"):
        ow.convert.run(short, tmp_path / "b.hic")
    (tmp_path / "sizes").write_text("".join(f"{c}\t{n}\n" for c, n in sizes.items()))
    ow.convert.run(
        short, tmp_path / "b.hic", chrom_sizes=tmp_path / "sizes", resolutions=[10_000, 50_000]
    )
    assert (tmp_path / "b.hic").read_bytes() == (tmp_path / "a.hic").read_bytes()
    with pytest.raises(ow.errors.InvalidError, match="not a multiple of the finest"):
        ow.convert.run(src, tmp_path / "c.hic", resolutions=[10_000, 25_000])
    with pytest.raises(ow.errors.InvalidError, match="not a multiple of the finest"):
        ow.convert.plan(src, tmp_path / "c.hic", resolutions=[10_000, 25_000])
    assert ow.convert.plan(src, tmp_path / "c.hic", memory="256M").memory.buffer > 0
    assert Pairs(tmp_path / "x.pairs").format == "pairs"


def test_bedmethyl_to_bigbed_and_bigwig(tmp_path):
    from omicsweave.formats import BedMethyl, BigBed

    rng = np.random.default_rng(4)
    lines = []
    for chrom, length in (("chr1", 50_000), ("chr2", 30_000)):
        for pos in sorted(rng.choice(np.arange(0, length - 2, 2), 300, replace=False).tolist()):
            for strand, at in (("+", pos), ("-", pos + 1)):
                for code in ("m", "h"):
                    valid = int(rng.integers(1, 40))
                    nmod = int(rng.integers(0, valid + 1))
                    counts = f"{valid} {100 * nmod / valid:.2f} {nmod} {valid - nmod} 0 0 1 0 0"
                    lines.append(
                        (chrom, at, code, strand, valid, nmod, f"{100 * nmod / valid:.2f}", counts)
                    )
    rng.shuffle(lines)
    # Older modkit: the nine counts separated by spaces.
    text = "".join(
        f"{c}\t{a}\t{a + 1}\t{m}\t{v}\t{s}\t{a}\t{a + 1}\t255,0,0\t{counts}\n"
        for c, a, m, s, v, _, _, counts in lines
    )
    src = tmp_path / "calls.bed"
    src.write_text(text)
    assert detect(src) == BedMethyl(src)

    report = ow.convert.run(src, tmp_path / "calls.bb", threads=2, memory="64M")
    assert report.target == BigBed(tmp_path / "calls.bb") and report.records == len(lines)
    with IntervalReader(report.target) as f:
        got = f.fetch(("chr2", 0, f.chroms["chr2"]))
        want = sorted(
            (a, f"{m}\t{v}\t{s}\t{a}\t{a + 1}\t255,0,0\t" + counts.replace(" ", "\t"))
            for c, a, m, s, v, _, _, counts in lines
            if c == "chr2"
        )
        assert sorted(zip(got.start.tolist(), got.annotation.to_list(), strict=True)) == want

    with pytest.raises(ow.errors.InvalidError, match=r"modification codes \(h, m\)"):
        ow.convert.run(src, tmp_path / "calls.bw")
    with pytest.raises(ow.errors.InvalidError, match="does not apply"):
        ow.convert.run(src, tmp_path / "calls.bb", mod="m")
    report = ow.convert.run(src, tmp_path / "calls.bw", mod="m", strand="+", value="nmod")
    with TrackReader(report.target) as f:
        got = f.intervals(("chr1", 0, f.chroms["chr1"]))
        want = sorted(
            (a, a + 1, float(n))
            for c, a, m, s, _, n, _, _ in lines
            if c == "chr1" and m == "m" and s == "+"
        )
        assert (
            list(zip(got.start.tolist(), got.end.tolist(), got.value.tolist(), strict=True)) == want
        )
    report = ow.convert.run(src, tmp_path / "pct.bw", mod="h")
    with TrackReader(report.target) as f:
        got = f.intervals(("chr2", 0, f.chroms["chr2"]))
        want = sorted(
            (a, float(np.float32(float(p))))
            for c, a, m, _, _, _, p, _ in lines
            if c == "chr2" and m == "h"
        )
        assert list(zip(got.start.tolist(), got.value.tolist(), strict=True)) == want
