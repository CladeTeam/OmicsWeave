import pickle
import struct

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.formats import BgzipFasta, Fai, Fasta, FastaWriter, Gzi, TwoBit
from omicsweave.query import SequenceReader
from omicsweave.types import Region


def random_genome(rng: np.random.Generator, n: int) -> dict[str, bytes]:
    alphabet = np.frombuffer(b"ACGTacgtN", dtype=np.uint8)
    return {
        f"chr{i}": alphabet[rng.integers(0, len(alphabet), int(rng.integers(1, 3000)))].tobytes()
        for i in range(n)
    }


def write_fasta(path, genome, **kw):
    with FastaWriter(path, **kw) as w:
        for name, seq in genome.items():
            w.add(name, seq)
    return w


@pytest.mark.parametrize("suffix", [".fa", ".fa.gz"])
@pytest.mark.parametrize("line_width", [0, 1, 7, 60])
def test_written_fasta_reads_back(tmp_path, suffix, line_width):
    rng = np.random.default_rng(line_width)
    genome = random_genome(rng, 5)
    path = tmp_path / f"g{suffix}"
    w = FastaWriter(path, line_width=line_width)
    for name, seq in genome.items():
        w.add(name, seq)
    companions = w.close()
    assert w.close() == ()
    if suffix == ".fa.gz":
        assert w.file == BgzipFasta(path)
        assert companions == (Fai(f"{path}.fai"), Gzi(f"{path}.gzi"))
    else:
        assert w.file == Fasta(path)
        assert companions == (Fai(f"{path}.fai"),)
    assert all(c.exists() for c in companions)
    assert not Gzi.beside(path).exists() or suffix == ".fa.gz"
    f = SequenceReader(path)
    assert f.source == w.file
    assert dict(f.chroms) == {k: len(v) for k, v in genome.items()}
    for _ in range(50):
        name = list(genome)[int(rng.integers(len(genome)))]
        n = len(genome[name])
        start = int(rng.integers(0, n + 1))
        end = int(rng.integers(start, n + 1))
        s = f.fetch(Region(name, start, end))
        assert s.bases.dtype == np.uint8
        assert s.bases.tobytes() == genome[name][start:end]
        assert s.text == genome[name][start:end].decode()


def test_indexes_are_required_and_built_explicitly(tmp_path):
    genome = random_genome(np.random.default_rng(1), 3)
    path = tmp_path / "g.fa.gz"
    write_fasta(path, genome)
    fai, gzi = tmp_path / "g.fa.gz.fai", tmp_path / "g.fa.gz.gzi"
    written = fai.read_bytes(), gzi.read_bytes()
    gzi.unlink()
    with pytest.raises(OSError, match="build_index") as e:
        SequenceReader(path)
    assert isinstance(e.value, errors.IoError) and e.value.kind == "NotFound"
    fai.unlink()
    with pytest.raises(errors.IoError):
        SequenceReader(path)
    assert not fai.exists() and not gzi.exists()
    assert BgzipFasta(path).build_index() == (Fai(fai), Gzi(gzi))
    assert (fai.read_bytes(), gzi.read_bytes()) == written
    f = SequenceReader(path)
    assert f.fetch(("chr0", 0, 10)).bases.tobytes() == genome["chr0"][:10]


def test_build_index_matches_written_index(tmp_path):
    genome = random_genome(np.random.default_rng(2), 4)
    path = tmp_path / "g.fa"
    write_fasta(path, genome, line_width=13)
    written = (tmp_path / "g.fa.fai").read_bytes()
    (tmp_path / "g.fa.fai").unlink()
    assert Fasta(path).build_index() == Fai(tmp_path / "g.fa.fai")
    assert (tmp_path / "g.fa.fai").read_bytes() == written


def test_explicit_index_paths(tmp_path):
    genome = random_genome(np.random.default_rng(3), 2)
    path = tmp_path / "g.fa"
    write_fasta(path, genome)
    moved = tmp_path / "elsewhere.fai"
    (tmp_path / "g.fa.fai").rename(moved)
    f = SequenceReader(path, Fai(moved))
    assert f.fetch(("chr1", 0, 5)).bases.tobytes() == genome["chr1"][:5]
    assert not (tmp_path / "g.fa.fai").exists()
    with pytest.raises(OSError) as e:
        SequenceReader(path, fai=tmp_path / "missing.fai")
    assert isinstance(e.value, errors.IoError)
    assert e.value.kind == "NotFound"


def test_out_argument_is_filled_in_place(tmp_path):
    genome = {"a": b"ACGTACGTAC"}
    path = tmp_path / "g.fa"
    write_fasta(path, genome)
    f = SequenceReader(path)
    out = np.zeros(4, dtype=np.uint8)
    s = f.fetch(("a", 2, 6), out=out)
    assert s.bases is out
    assert out.tobytes() == b"GTAC"
    with pytest.raises(errors.OutputSizeError) as e:
        f.fetch(("a", 2, 6), out=np.zeros(3, dtype=np.uint8))
    assert (e.value.have, e.value.need) == (3, 4)


def test_errors(tmp_path):
    path = tmp_path / "g.fa"
    write_fasta(path, {"a": b"ACGT"})
    f = SequenceReader(path)
    with pytest.raises(LookupError) as e:
        f.fetch(("b", 0, 1))
    assert isinstance(e.value, errors.UnknownChromError)
    assert e.value.chrom == "b"
    with pytest.raises(IndexError) as e:
        f.fetch(("a", 2, 9))
    assert isinstance(e.value, errors.RegionOutOfRangeError)
    assert e.value.length == 4
    with pytest.raises(errors.InvalidRegionError):
        Region("a", 3, 2)
    with pytest.raises(errors.UnknownChromError):
        f.length("b")
    bad = tmp_path / "bad.fa"
    bad.write_bytes(b">a\nACGT\nAC\nACGT\n")
    with pytest.raises(errors.FormatError) as e:
        Fasta(bad).build_index()
    assert e.value.line >= 1
    f.close()
    with pytest.raises(ValueError, match="closed"):
        f.fetch(("a", 0, 1))


def test_writer_rejects_bad_records(tmp_path):
    with FastaWriter(tmp_path / "g.fa") as w:
        w.add("a", "ACGT")
        with pytest.raises(errors.InvalidError):
            w.add("a", "ACGT")
        with pytest.raises(errors.InvalidError):
            w.add("has space", "ACGT")
        with pytest.raises(errors.InvalidError):
            w.add("empty", "")


def test_pickle_reopens(tmp_path):
    path = tmp_path / "g.fa"
    write_fasta(path, {"a": b"ACGTACGT"})
    f = SequenceReader(path)
    g = pickle.loads(pickle.dumps(f))
    assert g is not f
    assert g.fetch(("a", 1, 4)).text == "CGT"


def write_2bit(path, genome: dict[str, str]) -> None:
    """A minimal 2bit file: N runs as N blocks, lowercase runs as mask blocks."""
    code = {"T": 0, "C": 1, "A": 2, "G": 3, "N": 0}

    def runs(seq, test):
        out, i = [], 0
        while i < len(seq):
            if test(seq[i]):
                j = i
                while j < len(seq) and test(seq[j]):
                    j += 1
                out.append((i, j - i))
                i = j
            else:
                i += 1
        return out

    records = []
    for seq in genome.values():
        n_runs = runs(seq, lambda c: c in "Nn")
        mask = runs(seq, str.islower)
        upper = seq.upper()
        packed = bytearray()
        for i in range(0, len(upper), 4):
            chunk = upper[i : i + 4].ljust(4, "T")
            b = 0
            for c in chunk:
                b = (b << 2) | code[c]
            packed.append(b)
        rec = struct.pack("<II", len(seq), len(n_runs))
        rec += b"".join(struct.pack("<I", s) for s, _ in n_runs)
        rec += b"".join(struct.pack("<I", n) for _, n in n_runs)
        rec += struct.pack("<I", len(mask))
        rec += b"".join(struct.pack("<I", s) for s, _ in mask)
        rec += b"".join(struct.pack("<I", n) for _, n in mask)
        rec += struct.pack("<I", 0) + bytes(packed)
        records.append(rec)
    header = struct.pack("<IIII", 0x1A412743, 0, len(genome), 0)
    index_len = sum(1 + len(name) + 4 for name in genome)
    offset = len(header) + index_len
    index = b""
    for name, rec in zip(genome, records, strict=True):
        index += bytes([len(name)]) + name.encode() + struct.pack("<I", offset)
        offset += len(rec)
    path.write_bytes(header + index + b"".join(records))


def test_two_bit(tmp_path):
    genome = {"x": "ACGTNNacgtAC", "y": "GGGTTTAAACCC"}
    path = tmp_path / "g.2bit"
    write_2bit(path, genome)
    f = ow.open(path)
    assert isinstance(f, SequenceReader)
    assert f.source == TwoBit(path)
    assert dict(f.chroms) == {"x": 12, "y": 12}
    assert f.fetch(("x", 0, 12)).text == "ACGTNNACGTAC"
    masked = SequenceReader(path, soft_mask=True)
    assert masked.fetch(("x", 0, 12)).text == "ACGTNNacgtAC"
