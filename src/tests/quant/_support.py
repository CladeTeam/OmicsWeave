import gzip
import random

import pytest

import omicsweave as ow
from omicsweave.annotation import Annotation

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


def _genome(rng):
    return {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CHROMS}


def _sam_to_bam(tmp_path, name, lines):
    head = ["@HD\tVN:1.6\tSO:coordinate\n"] + [f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in CHROMS]
    sam = tmp_path / f"{name}.sam"
    sam.write_text("".join(head + lines))
    bam = tmp_path / f"{name}.bam"
    ow.convert.run(sam, bam, threads=1)
    return bam
