import gzip
import json
import random

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import quant
from omicsweave.cli.app import main
from omicsweave.errors import OmicsWeaveError
from omicsweave.query import TrackReader

from ._support import CHROMS, _genome, _layer, _sam_to_bam
from ._support import data as data


def _fasta(tmp_path, genome):
    fa = tmp_path / "genome.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    ref = tmp_path / "genome.fa.gz"
    ow.convert.run(fa, ref, threads=1)
    return ref


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
