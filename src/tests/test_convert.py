import gzip

import numpy as np
import pytest

import omicsweave as ow
from omicsweave.formats import Bcf, BgzipFasta, Csi, Fai, Fasta, Gzi, GzipFasta, TwoBit, Vcf
from omicsweave.query import IntervalReader, SequenceReader


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
