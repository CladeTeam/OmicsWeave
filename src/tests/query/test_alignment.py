import pickle

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.errors import InvalidError, IoError, UnknownChromError
from omicsweave.formats import Bai, Bam, Csi
from omicsweave.query import (
    AlignmentReader,
)
from omicsweave.types import AlignmentField

from ..support.alignment import CHROMS, expected, sam_text, span

# name, flag, chrom, 1-based pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags


@pytest.fixture
def bam(tmp_path):
    sam = tmp_path / "in.sam"
    sam.write_text(sam_text())
    out = tmp_path / "in.bam"
    report = ow.convert.run(sam, out, index="bai", threads=2, memory="64M")
    assert isinstance(report.companions[0], Bai)
    return Bam(out)


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


def test_cram_decoding_options_belong_to_cram(tmp_path):
    sam = tmp_path / "r.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100\nr\t0\tc\t5\t9\t4M\t*\t0\t0\tACGT\t*\n"
    )
    bam = ow.convert.run(sam, tmp_path / "r.bam", index=True, threads=1).target
    with AlignmentReader(bam) as reads:
        assert len(reads.fetch(("c", 0, 100))) == 1
    with pytest.raises(errors.InvalidError, match="CRAM decoding"):
        AlignmentReader(bam, verify_checksums=False)
