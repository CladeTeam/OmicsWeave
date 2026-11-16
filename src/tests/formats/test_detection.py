import gzip
import pickle

import pytest

from omicsweave import errors
from omicsweave.formats import (
    Bcf,
    Bed,
    BgzipFasta,
    Fai,
    Fasta,
    GzipBed,
    GzipVcf,
    Sidecar,
    Vcf,
    detect,
)


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
