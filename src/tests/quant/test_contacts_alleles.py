import gzip
import random

import omicsweave as ow
from omicsweave import quant

from ._support import CHROMS, _genome, _layer, _sam_to_bam
from ._support import data as data


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
