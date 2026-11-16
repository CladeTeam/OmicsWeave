import gzip

from omicsweave.cli.app import main

from ._support import data as data


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
