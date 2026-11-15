import gzip
import pickle
import shutil
from pathlib import Path

import numpy as np
import pytest
from test_alignment import CHROMS, RECORDS, expected, span

import omicsweave as ow
from omicsweave.errors import InvalidError, IoError, UnknownChromError
from omicsweave.formats import Bam, Crai, Cram, Fasta, detect
from omicsweave.query import AlignmentReader, SequenceReader

# Written by htslib 1.25 from the records of test_alignment against ref.fa:
# ref.cram (CRAM 3.1, external reference) and embed.cram (CRAM 3.0, the
# reference embedded), each with the .crai samtools index writes.
DATA = Path(__file__).parent / "data" / "cram"


@pytest.fixture
def files(tmp_path):
    for name in ("ref.fa", "ref.cram", "ref.cram.crai", "embed.cram", "embed.cram.crai"):
        shutil.copy(DATA / name, tmp_path / name)
    fasta = Fasta(tmp_path / "ref.fa")
    fasta.build_index()
    return fasta, Cram(tmp_path / "ref.cram"), Cram(tmp_path / "embed.cram")


REGIONS = [("chr1", 0, 5000), ("chr1", 104, 106), ("chr1", 160, 250), ("chr1", 4999, 6000)]


@pytest.mark.parametrize("region", REGIONS)
def test_fetch_matches_the_records(files, region):
    fasta, cram, _ = files
    with AlignmentReader(cram, reference=SequenceReader(fasta)) as reader:
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
    assert got.start.dtype == np.int64 and got.cigars.ops.dtype == np.uint32


def test_same_as_the_bam_converted_from_it(files, tmp_path):
    fasta, cram, embed = files
    bam = tmp_path / "from.bam"
    report = ow.convert.run(cram, bam, reference=fasta, index="bai", threads=2, memory="256M")
    assert isinstance(report.target, Bam) and report.records == len(RECORDS)
    bam_reader = AlignmentReader(bam)
    for source, reference in [(cram, fasta), (cram, fasta.path), (embed, None)]:
        reader = AlignmentReader(source, reference=reference)
        assert reader.chroms == CHROMS
        assert reader.header_text == bam_reader.header_text
        for region in REGIONS + [("chr2", 0, 3000)]:
            a, b = reader.fetch(region, "all"), bam_reader.fetch(region, "all")
            for field in ("start", "end", "flag", "mapq", "mate_tid", "mate_start", "tlen"):
                assert getattr(a, field).tolist() == getattr(b, field).tolist(), field
            for field in ("names", "seqs", "quals", "aux"):
                # Byte columns: absent qualities are 0xff, binary aux is not text.
                x, y = getattr(a, field), getattr(b, field)
                assert [x.bytes(i) for i in range(len(x))] == [y.bytes(i) for i in range(len(y))], (
                    field
                )
            assert a.cigars.ops.tolist() == b.cigars.ops.tolist()


def test_md_and_nm_are_generated(files):
    fasta, cram, _ = files
    reader = AlignmentReader(cram, reference=fasta)
    got = reader.fetch(("chr1", 100, 101), "aux")
    aux = got.aux.bytes(0)
    assert b"MDZ" in aux and b"NM" in aux


def test_fixed_fields_need_no_reference(files):
    _, cram, _ = files
    reader = AlignmentReader(cram)
    got = reader.fetch(("chr1", 0, 5000), ("cigar", "name", "qual"))
    want = expected("chr1", 0, 5000)
    assert got.names.to_list() == [r[0] for r in want]
    assert got.end.tolist() == [span(r)[1] for r in want]
    with pytest.raises(InvalidError, match="reference"):
        reader.fetch(("chr1", 0, 5000), "seq")


def test_formats_and_index(files, tmp_path):
    _, cram, _ = files
    assert isinstance(detect(cram), Cram)
    assert isinstance(detect(Crai.beside(cram)), Crai)
    copy = Cram(tmp_path / "copy.cram")
    copy.path.write_bytes(cram.path.read_bytes())
    with pytest.raises(IoError, match="build_index"):
        AlignmentReader(copy)
    crai = copy.build_index()
    assert isinstance(crai, Crai) and crai.path == tmp_path / "copy.cram.crai"
    assert gzip.decompress(crai.path.read_bytes()) == gzip.decompress(
        Crai.beside(cram).path.read_bytes()
    )
    reader = AlignmentReader(copy)
    assert isinstance(reader.index, Crai)


def test_errors_open_and_pickle(files):
    fasta, cram, _ = files
    reader = ow.open(cram, reference=fasta)
    assert isinstance(reader, AlignmentReader)
    with pytest.raises(UnknownChromError):
        reader.fetch(("chrX", 0, 10))
    assert len(reader.fetch(("chr1", 10, 10))) == 0
    again = pickle.loads(pickle.dumps(reader))
    assert again.fetch(("chr1", 0, 5000), "seq").seqs.to_list() == [
        r[9] for r in expected("chr1", 0, 5000)
    ]
    clone = reader.clone()
    assert len(clone.fetch(("chr1", 0, 5000))) == len(expected("chr1", 0, 5000))
    with pytest.raises(ow.errors.UnsupportedError):
        reader.cache_namespace  # noqa: B018
    reader.close()
    with pytest.raises(ValueError, match="closed"):
        reader.fetch(("chr1", 0, 10))
