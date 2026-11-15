import pickle
import re

import numpy as np
import pytest

import omicsweave as ow
from omicsweave.errors import InvalidError, IoError, UnknownChromError
from omicsweave.formats import Bai, Bam, Csi
from omicsweave.query import AlignmentReader
from omicsweave.types import AlignmentField

CHROMS = {"chr1": 5000, "chr2": 3000}

# name, flag, chrom, 1-based pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags
RECORDS = [
    ("a", 99, "chr1", 101, 60, "10M", "=", 251, 160, "ACGTACGTAC", "IIIIIIIIII", "NM:i:0"),
    ("b", 0, "chr1", 105, 30, "3S5M2D4M", "*", 0, 0, "GGGACGTAAAAA", "*", ""),
    ("c", 16, "chr1", 150, 10, "4M100N4M", "*", 0, 0, "ACGTACGT", "!!!!####", "XS:A:-"),
    ("d", 1024, "chr1", 200, 60, "8M", "*", 0, 0, "ACGTACGT", "*", ""),
    ("u", 69, "chr1", 220, 0, "*", "=", 101, 0, "ACGT", "*", ""),
    ("a", 147, "chr1", 251, 60, "10M", "=", 101, -160, "ACGTACGTAC", "IIIIIIIIII", ""),
    ("e", 0, "chr1", 4995, 60, "6M", "*", 0, 0, "ACGTAC", "*", ""),
    ("f", 0, "chr2", 10, 60, "5M", "*", 0, 0, "ACGTA", "*", ""),
    ("x", 4, "*", 0, 0, "*", "*", 0, 0, "ACGT", "*", ""),
]


def sam_text() -> str:
    head = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
        f"@SQ\tSN:{n}\tLN:{ln}\n" for n, ln in CHROMS.items()
    )
    body = ""
    for r in RECORDS:
        fields = [str(x) for x in r[:11]]
        if r[11]:
            fields.append(r[11])
        body += "\t".join(fields) + "\n"
    return head + body


def ref_len(cigar: str) -> int:
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X")


def span(r) -> tuple[int, int]:
    start = r[3] - 1
    rlen = 0 if r[1] & 4 or r[5] == "*" else ref_len(r[5])
    return start, start + (rlen or 1)


@pytest.fixture
def bam(tmp_path):
    sam = tmp_path / "in.sam"
    sam.write_text(sam_text())
    out = tmp_path / "in.bam"
    report = ow.convert.run(sam, out, index="bai", threads=2, memory="64M")
    assert isinstance(report.companions[0], Bai)
    return Bam(out)


def expected(chrom: str, start: int, end: int):
    return [r for r in RECORDS if r[2] == chrom and span(r)[1] > start and span(r)[0] < end]


@pytest.mark.parametrize(
    "region",
    [
        ("chr1", 0, 5000),
        ("chr1", 104, 106),
        ("chr1", 160, 250),
        ("chr1", 4999, 6000),
        ("chr2", 0, 9),
    ],
)
def test_fetch_matches_the_records(bam, region):
    with AlignmentReader(bam) as reader:
        got = reader.fetch(region, "all")
    want = expected(*region)
    assert len(got) == len(want)
    assert got.names.to_list() == [r[0] for r in want]
    assert got.flag.tolist() == [r[1] for r in want]
    assert got.start.tolist() == [span(r)[0] for r in want]
    assert got.end.tolist() == [span(r)[1] for r in want]
    assert got.mapq.tolist() == [r[4] for r in want]
    assert [got.cigars.string(i) for i in range(len(got))] == [r[5] for r in want]
    assert got.seqs.to_list() == [r[9] for r in want]
    for i, r in enumerate(want):
        qual = got.quals.bytes(i)
        if r[10] == "*":
            assert qual == b"\xff" * len(r[9])
        else:
            assert qual == bytes(ord(c) - 33 for c in r[10])
    assert got.start.dtype == np.int64 and got.flag.dtype == np.uint16
    assert got.mate_tid.dtype == np.int32 and got.cigars.ops.dtype == np.uint32


def test_fields_are_optional(bam):
    reader = AlignmentReader(bam)
    bare = reader.fetch(("chr1", 0, 5000))
    assert bare.names is None and bare.seqs is None and bare.cigars is None
    some = reader.fetch(("chr1", 0, 5000), AlignmentField.NAME | AlignmentField.AUX)
    assert some.names is not None and some.aux is not None and some.quals is None
    assert some.aux.bytes(0) == b"NMC\x00"  # htslib stores 0 in one byte
    assert reader.fetch(("chr1", 0, 5000), ("cigar",)).cigars is not None
    with pytest.raises(ValueError, match="unknown alignment field"):
        reader.fetch(("chr1", 0, 10), ("names",))


def test_mates_and_header(bam):
    reader = AlignmentReader(bam)
    assert reader.chroms == CHROMS
    assert reader.references == ("chr1", "chr2")
    assert reader.header_text.startswith(b"@HD\tVN:1.6")
    got = reader.fetch(("chr1", 100, 101))
    assert reader.references[got.mate_tid[0]] == "chr1"
    assert got.mate_start.tolist() == [250] and got.tlen.tolist() == [160]


def test_errors(bam, tmp_path):
    reader = AlignmentReader(bam)
    with pytest.raises(UnknownChromError):
        reader.fetch(("chrX", 0, 10))
    assert len(reader.fetch(("chr1", 10, 10))) == 0
    bare = tmp_path / "bare.bam"
    bare.write_bytes(bam.path.read_bytes())
    with pytest.raises(IoError, match="build_index"):
        AlignmentReader(bare)
    reader.close()
    with pytest.raises(ValueError, match="closed"):
        reader.fetch(("chr1", 0, 10))


def test_built_indexes(bam, tmp_path):
    converted = Bai.beside(bam).path.read_bytes()
    copy = Bam(tmp_path / "copy.bam")
    copy.path.write_bytes(bam.path.read_bytes())
    bai = copy.build_index("bai", threads=2)
    assert isinstance(bai, Bai) and bai.path.read_bytes() == converted
    csi = copy.build_index("csi", min_shift=12)
    assert isinstance(csi, Csi) and csi.path == tmp_path / "copy.bam.csi"
    # htslib looks for the CSI first.
    reader = AlignmentReader(copy)
    assert isinstance(reader.index, Csi)
    with_bai = AlignmentReader(copy, index=bai)
    for region in [("chr1", 0, 5000), ("chr1", 140, 160)]:
        assert reader.fetch(region, "name").names.to_list() == (
            with_bai.fetch(region, "name").names.to_list()
        )
    with pytest.raises(InvalidError):
        copy.build_index("tbi")  # type: ignore[arg-type]


def test_open_and_pickle(bam):
    reader = ow.open(bam)
    assert isinstance(reader, AlignmentReader)
    again = pickle.loads(pickle.dumps(reader))
    assert again.fetch(("chr1", 0, 5000), "name").names.to_list() == [
        r[0] for r in expected("chr1", 0, 5000)
    ]
