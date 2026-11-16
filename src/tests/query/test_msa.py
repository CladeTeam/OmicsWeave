import gzip
import pickle

import numpy as np
import pytest

import omicsweave as ow
from omicsweave.formats import BigBed, BigBedSchema, BigBedWriter, GzipMaf, Maf, detect
from omicsweave.query import IntervalReader, MsaReader
from omicsweave.types import Region

SPECIES = ("hg38", "panTro4", "mm10", "canFam3")
CHROMS = {"chr1": 6_000, "chr2": 3_000}
FILL = ord("-")
ABSENT = 0


def expected(rows, fill=FILL):
    """``rows`` of the truth with absent bases as ``fill``."""
    return np.where(rows == ABSENT, fill, rows).astype(np.uint8)


def alignment(rng):
    """MAF text over both chromosomes, and the expected matrix of each:
    one row per species, one column per reference base, 0 where no block
    gives a species a base."""
    truth = {c: np.full((len(SPECIES), n), ABSENT, dtype=np.uint8) for c, n in CHROMS.items()}
    blocks = []
    for chrom, length in CHROMS.items():
        at = int(rng.integers(0, 50))
        while True:
            size = int(rng.integers(1, 60))
            if at + size > length:
                break
            ref = rng.choice(list(b"ACGTacgt"), size).astype(np.uint8)
            rows = {0: ref}
            for s in range(1, len(SPECIES)):
                if rng.random() < 0.7:
                    rows[s] = np.where(
                        rng.random(size) < 0.15, FILL, rng.choice(list(b"ACGTN"), size)
                    ).astype(np.uint8)
            for s, bases in rows.items():
                truth[chrom][s, at : at + size] = bases
            # Columns the reference lacks: insertions in the other species.
            cols = {s: list(bases) for s, bases in rows.items()}
            for _ in range(int(rng.integers(0, 4))):
                where = int(rng.integers(0, size + 1))
                for s, col in cols.items():
                    col.insert(where, FILL if s == 0 else int(rng.choice(list(b"ACGT"))))
            lines = [f"a score={rng.random() * 1000:.3f}"]
            for s, col in cols.items():
                text = bytes(col).decode()
                n = sum(c != "-" for c in text)
                strand = "+" if s == 0 else str(rng.choice(["+", "-"]))
                src_size = length if s == 0 else 10**7
                start = at if s == 0 else int(rng.integers(0, 10**6))
                name = chrom if s == 0 else "chrUn"
                lines.append(f"s {SPECIES[s]}.{name} {start} {n} {strand} {src_size} {text}")
            if len(cols) < len(SPECIES):
                lines.append(f"e {SPECIES[3]}.chr9 5 3 + 900 I")
            blocks.append("\n".join(lines))
            at += size + int(rng.integers(0, 25))
    order = rng.permutation(len(blocks))
    text = "##maf version=1\n\n" + "\n\n".join(blocks[i] for i in order) + "\n"
    return text, truth


@pytest.fixture(scope="module")
def bigmaf(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("maf")
    text, truth = alignment(np.random.default_rng(11))
    src = tmp / "multiz.maf.gz"
    src.write_bytes(gzip.compress(text.encode()))
    assert detect(src) == GzipMaf(src)
    report = ow.convert.run(src, tmp / "multiz.bigmaf", threads=2, memory="64M")
    assert report.target == BigBed(tmp / "multiz.bigmaf")
    return report.target, truth


def test_rows_follow_the_reference(bigmaf):
    path, truth = bigmaf
    rng = np.random.default_rng(3)
    with MsaReader(path, SPECIES) as f:
        assert dict(f.chroms) == CHROMS
        assert f.species == SPECIES
        for _ in range(60):
            chrom = str(rng.choice(list(CHROMS)))
            start = int(rng.integers(0, CHROMS[chrom] - 1))
            end = int(rng.integers(start + 1, min(start + 700, CHROMS[chrom]) + 1))
            want = expected(truth[chrom][:, start:end])
            got = f.fetch((chrom, start, end))
            assert got.species == SPECIES
            np.testing.assert_array_equal(got.bases, want)
            np.testing.assert_array_equal(got.row("mm10"), want[2])
            rev = f.fetch(Region(chrom, start, end, "-")).bases
            comp = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
            for r, w in zip(rev, want, strict=True):
                assert r.tobytes() == w.tobytes().translate(comp)[::-1]


def test_species_fill_out_and_pickling(bigmaf):
    path, truth = bigmaf
    with MsaReader(path, ["mm10", "hg38", "dog"], fill=".") as f:
        got = f.fetch(("chr2", 0, 3_000))
        assert got.species == ("mm10", "hg38", "dog")
        np.testing.assert_array_equal(got.bases[:2], expected(truth["chr2"][[2, 0]], ord(".")))
        assert set(got.bases[2].tolist()) == {ord(".")}
        out = np.zeros((3, 100), dtype=np.uint8)
        same = f.fetch(("chr2", 500, 600), out=out)
        assert same.bases is out
        np.testing.assert_array_equal(out, got.bases[:, 500:600])
        with pytest.raises(ValueError, match="shape"):
            f.fetch(("chr2", 500, 600), out=np.zeros((3, 99), dtype=np.uint8))
        out.flags.writeable = False
        with pytest.raises(ValueError, match="writable"):
            f.fetch(("chr2", 500, 600), out=out)
        again = pickle.loads(pickle.dumps(f))
        np.testing.assert_array_equal(again.fetch(("chr2", 0, 3_000)).bases, got.bases)
        again.close()
        clone = f.clone()
        np.testing.assert_array_equal(clone.fetch(("chr2", 0, 3_000)).bases, got.bases)


def test_rows_by_appearance_without_species(bigmaf):
    path, truth = bigmaf
    with MsaReader(path) as f:
        assert f.species is None
        got = f.fetch(("chr1", 0, 6_000))
        assert got.species[0] == "hg38" and set(got.species) == set(SPECIES)
        for name in SPECIES:
            np.testing.assert_array_equal(
                got.row(name), expected(truth["chr1"][SPECIES.index(name)])
            )
        with pytest.raises(ValueError, match="species fixed"):
            f.fetch(("chr1", 0, 10), out=np.zeros((4, 10), dtype=np.uint8))
        empty = f.fetch(("chr1", 5_999, 6_000))
        assert empty.bases.shape[1] == 1


def test_maf_text_and_other_bigbeds_are_refused(tmp_path, bigmaf):
    src = tmp_path / "a.maf"
    src.write_text("##maf version=1\na score=1\ns hg38.chr1 0 2 + 10 AC\n")
    assert detect(src) == Maf(src)
    with pytest.raises(ow.errors.InvalidError, match="ow convert"):
        MsaReader(src)
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t0\t5\tx\n")
    bb = ow.convert.run(bed, tmp_path / "a.bb").target
    with pytest.raises(ow.errors.InvalidError, match="not a bigMaf"):
        MsaReader(bb)
    with pytest.raises(ValueError, match="one ASCII character"):
        MsaReader(bigmaf[0], fill="--")


def test_reference_species_option(tmp_path):
    src = tmp_path / "x.maf"
    src.write_text(
        "##maf version=1\n\na score=2\ns mm10.chr3 5 2 + 50 AC\ns hg38.chr1 0 2 + 9 AC\n"
    )
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t0\t5\n")
    with pytest.raises(ow.errors.InvalidError, match="does not apply"):
        ow.convert.plan(bed, tmp_path / "a.bb", ref_species="mm10")
    report = ow.convert.run(src, tmp_path / "x.bb", ref_species="mm10")
    with IntervalReader(report.target) as f:
        assert dict(f.chroms) == {"chr3": 50}
        assert f.schema.fields[3].name == "mafBlock"
        got = f.fetch(("chr3", 0, 50))
        assert got.annotation.to_list() == [
            "a score=2.000000;s mm10.chr3 5 2 + 50 AC;s hg38.chr1 0 2 +  9 AC;"
        ]
    with pytest.raises(ow.errors.FormatError, match="line 3: .*reference hg38"):
        ow.convert.run(src, tmp_path / "y.bb", ref_species="hg38")


def test_bigmaf_checks_each_alignment_row(tmp_path):
    path = tmp_path / "bad.bigmaf"
    schema = BigBedSchema(
        4,
        auto_sql='table bigMaf\n"alignment"\n(\nstring chrom;\nuint chromStart;\n'
        "uint chromEnd;\nstring mafBlock;\n)\n",
    )
    block = "a;s hg38.c 0 3 + 10 AAA;s mm10.c 0 2 + 10 AA;s dog.c 0 4 + 10 AAAA;"
    with BigBedWriter(path, {"c": 10}, schema) as writer:
        writer.add("c", [0], [3], [block])
    with MsaReader(path) as reader, pytest.raises(ow.errors.InvalidError, match="2 columns, not 3"):
        reader.fetch(("c", 0, 3))
