import dataclasses
import errno
import gzip
import pickle
import threading

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.formats import (
    Bcf,
    Bed,
    BgzipFasta,
    Fai,
    Fasta,
    FastaWriter,
    GzipBed,
    GzipVcf,
    Sidecar,
    Vcf,
    detect,
)
from omicsweave.query import SequenceReader
from omicsweave.types import BinUnit, ContactWindow, Region


def test_region_and_window_keys():
    r = Region("chr1", np.int64(5), 9)
    assert (r.start, len(r), str(r)) == (5, 4, "chr1:5-9")
    assert dataclasses.replace(r, end=20).end == 20
    assert hash(r) == hash(Region("chr1", 5, 9))
    with pytest.raises(ValueError):
        Region("chr1", -1, 3)
    with pytest.raises(TypeError):
        Region("chr1", 1.5, 3)
    w = ContactWindow(("a", 0, 2500), ("b", 999, 1001), 1000, "FRAG")
    assert w.unit is BinUnit.FRAG
    assert w.rows == Region("a", 0, 2500)
    assert w.shape == (3, 2)
    assert ContactWindow(("a", 7, 7), ("a", 0, 1), 5).shape == (0, 1)


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


def test_files_are_path_like_values(tmp_path):
    fa = Fasta(tmp_path / "a.fa")
    assert fa.path == tmp_path / "a.fa" and str(tmp_path / "a.fa") == fa.__fspath__()
    assert fa == Fasta(str(tmp_path / "a.fa")) and fa != BgzipFasta(fa)
    assert len({fa, Fasta(fa.path), Fai(fa.path)}) == 2
    assert pickle.loads(pickle.dumps(fa)) == fa
    assert Fai.beside(fa) == Fai(tmp_path / "a.fa.fai")
    assert Sidecar.beside(fa).path.name == "a.fa.sidecar"
    assert not fa.exists()


def test_detect_reads_content_not_names(tmp_path):
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    cases = {
        "a.fa": (b">x\nACGT\n", Fasta),
        "a.bed": (b"x\t0\t1\n", Bed),
        "noext": (b"x\t0\t1\n", Bed),
        "a.txt": (vcf.encode(), Vcf),
        "a.vcf.gz": (gzip.compress(vcf.encode()), GzipVcf),
        "a.gz": (gzip.compress(b"x\t0\t1\n"), GzipBed),
        "a.bin": (gzip.compress(b"BCF\x02\x02"), Bcf),
    }
    for name, (content, cls) in cases.items():
        (tmp_path / name).write_bytes(content)
        assert detect(tmp_path / name) == cls(tmp_path / name)
    (tmp_path / "bad.gz").write_bytes(b"\x1f\x8b\x08\x00" + b"\xff" * 20)
    with pytest.raises(errors.OmicsWeaveError):
        detect(tmp_path / "bad.gz")


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


def test_threads_share_a_reader(tmp_path):
    rng = np.random.default_rng(0)
    seq = np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, 200_000)].tobytes()
    with FastaWriter(tmp_path / "g.fa") as w:
        w.add("c", seq)
    f = SequenceReader(tmp_path / "g.fa")
    failures = []

    def work(seed):
        r = np.random.default_rng(seed)
        for _ in range(200):
            s = int(r.integers(0, len(seq)))
            e = int(r.integers(s, len(seq) + 1))
            if f.fetch(("c", s, e)).bases.tobytes() != seq[s:e]:
                failures.append((s, e))

    threads = [threading.Thread(target=work, args=(i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert failures == []
