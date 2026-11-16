import numpy as np
import pytest

import omicsweave as ow
from omicsweave.errors import InvalidError
from omicsweave.formats import Bam
from omicsweave.query import AlignmentReader
from omicsweave.types import AlignmentField, ModificationFlag, Modifications

M, H = ord("m"), ord("h")

# name, flag, 1-based pos, cigar, seq, tags
RECORDS = [
    # C at 1, 3, 6, 9; the first is soft-clipped.
    ("fwd", 0, 300, "2S8M", "ACGCGTCAGC", "MM:Z:C+mh,0,1;\tML:B:C,200,10,50,60"),
    # Sequenced as TAGC: its C is the stored G at 0.
    ("rev", 16, 400, "4M", "GCTA", "MM:Z:C+m,0;\tML:B:C,99"),
    ("explicit", 0, 500, "2M1I3M", "CACCAA", "MM:Z:C+m?,1;A+17596,2;"),
    ("plain", 0, 600, "4M", "ACGT", ""),
]


@pytest.fixture
def bam(tmp_path):
    head = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5000\n"
    body = "".join(
        "\t".join([name, str(flag), "chr1", str(pos), "60", cigar, "*", "0", "0", seq, "*"])
        + (f"\t{tags}" if tags else "")
        + "\n"
        for name, flag, pos, cigar, seq, tags in RECORDS
    )
    sam = tmp_path / "mods.sam"
    sam.write_text(head + body)
    out = tmp_path / "mods.bam"
    ow.convert.run(sam, out, index="bai", threads=1, memory="64M")
    return Bam(out)


def rows(mods: Modifications, i: int) -> list[tuple]:
    s = mods.rows(i)
    return list(
        zip(
            mods.qpos[s].tolist(),
            mods.ref_pos[s].tolist(),
            [Modifications.code_name(c) for c in mods.code[s].tolist()],
            bytes(mods.canonical[s]).decode(),
            mods.probability[s].tolist(),
            mods.flags[s].tolist(),
            strict=True,
        )
    )


def test_modifications_decode_mm_and_ml(bam):
    with AlignmentReader(bam) as reader:
        got = reader.fetch(("chr1", 0, 5000), "modifications")
    mods = got.modifications
    assert got.seqs is None and len(mods) == 4
    assert mods.offsets.tolist() == [0, 4, 5, 7, 7]
    assert rows(mods, 0) == [
        (1, -1, "m", "C", 200, 0),
        (1, -1, "h", "C", 10, 0),
        (6, 303, "m", "C", 50, 0),
        (6, 303, "h", "C", 60, 0),
    ]
    assert rows(mods, 1) == [(0, 399, "m", "C", 99, 0)]
    unscored = int(ModificationFlag.EXPLICIT | ModificationFlag.UNSCORED)
    # The inserted base (query 2) is not on the reference.
    assert rows(mods, 2) == [
        (2, -1, "m", "C", 0, unscored),
        (5, 503, "17596", "A", 0, int(ModificationFlag.UNSCORED)),
    ]
    assert mods.record.tolist() == [0, 0, 0, 0, 1, 2, 2]
    assert mods.ref_pos.dtype == np.int64 and mods.qpos.dtype == np.int32
    assert mods.code.dtype == np.int32 and mods.flags.dtype == np.uint8


def test_modification_codes_filter_rows(bam):
    reader = AlignmentReader(bam)
    region = ("chr1", 0, 5000)
    only_h = reader.fetch(region, AlignmentField.MODIFICATIONS, modification_codes="h")
    assert only_h.modifications.code.tolist() == [H, H]
    assert only_h.modifications.offsets.tolist() == [0, 2, 2, 2, 2]
    chebi = reader.fetch(region, "modifications", modification_codes=[17596, "m"])
    assert chebi.modifications.code.tolist() == [M, M, M, M, -17596]
    assert reader.fetch(region, "all").modifications.offsets[-1] == 7
    with pytest.raises(ValueError, match="'modifications' field"):
        reader.fetch(region, "seq", modification_codes="m")
    with pytest.raises(ValueError, match="invalid modification code"):
        reader.fetch(region, "modifications", modification_codes="mh")


def test_large_modification_columns_preserve_each_base(tmp_path):
    n = 20_000
    sam = tmp_path / "dense.sam"
    sam.write_text(
        f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:{n}\n"
        f"dense\t0\tc\t1\t60\t{n}M\t*\t0\t0\t{'C' * n}\t*\t"
        f"MM:Z:C+m,{','.join(['0'] * n)};\tML:B:C,{','.join(str(i % 256) for i in range(n))}\n"
    )
    path = tmp_path / "dense.bam"
    ow.convert.run(sam, path, index=True, threads=1, memory="64M")
    with AlignmentReader(path) as reader:
        mods = reader.fetch(("c", 0, n), "modifications").modifications
    np.testing.assert_array_equal(mods.offsets, [0, n])
    np.testing.assert_array_equal(mods.ref_pos, np.arange(n, dtype=np.int64))
    np.testing.assert_array_equal(mods.qpos, np.arange(n, dtype=np.int32))
    np.testing.assert_array_equal(mods.probability, np.arange(n).astype(np.uint8))
    assert np.all(mods.code == M) and np.all(mods.canonical == ord("C"))
    assert np.all(mods.flags == 0)


def test_malformed_tags_fail_the_query(tmp_path):
    sam = tmp_path / "bad.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
        "bad\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\t*\tMM:Z:C+m,3;\n"
    )
    out = tmp_path / "bad.bam"
    ow.convert.run(sam, out, index="bai", threads=1, memory="64M")
    reader = AlignmentReader(out)
    with pytest.raises(InvalidError, match="bad"):
        reader.fetch(("chr1", 0, 100), "modifications")
    assert len(reader.fetch(("chr1", 0, 100), "aux")) == 1
