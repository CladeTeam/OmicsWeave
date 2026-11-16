import errno

import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.formats import (
    Bed,
    BgzipFasta,
    Fai,
    Fasta,
)
from omicsweave.query import SequenceReader


def test_error_hierarchy():
    assert issubclass(errors.UnknownChromError, LookupError)
    assert issubclass(errors.FormatError, ValueError)
    assert issubclass(errors.BigEndianError, NotImplementedError)
    e = errors.IoError("I/O (NotFound): gone", kind="NotFound", detail="gone")
    assert isinstance(e, OSError) and e.errno == errno.ENOENT
    for name in errors.__all__:
        assert issubclass(getattr(errors, name), errors.OmicsWeaveError)


def test_missing_file_is_os_error(tmp_path):
    with pytest.raises(OSError) as e:
        SequenceReader(tmp_path / "none.fa")
    assert e.value.errno == errno.ENOENT


def test_open_picks_the_reader(tmp_path):
    (tmp_path / "a.fa").write_text(">x\nACGT\n")
    with pytest.raises(errors.IoError):
        ow.open(tmp_path / "a.fa")
    assert Fasta(tmp_path / "a.fa").build_index() == Fai(tmp_path / "a.fa.fai")
    with ow.open(tmp_path / "a.fa") as f:
        assert isinstance(f, SequenceReader) and f.source == Fasta(tmp_path / "a.fa")
        assert f.fetch(("x", 1, 3)).text == "CG"
    assert f.closed
    (tmp_path / "v.vcf").write_text("##fileformat=VCFv4.2\n")
    with pytest.raises(errors.InvalidError, match="convert it to BCF"):
        ow.open(tmp_path / "v.vcf")


def test_file_objects_are_checked_against_content(tmp_path):
    (tmp_path / "a.fa").write_text(">x\nACGT\n")
    with pytest.raises(errors.InvalidError, match="holds fasta data, not BgzipFasta"):
        SequenceReader(BgzipFasta(tmp_path / "a.fa"))
    with pytest.raises(errors.InvalidError):
        BgzipFasta(tmp_path / "a.fa").build_index()
    with pytest.raises(TypeError, match="SequenceReader reads"):
        SequenceReader(Bed(tmp_path / "a.fa"))
