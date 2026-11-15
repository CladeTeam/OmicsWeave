import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from omicsweave.cli.app import main
from omicsweave.query import SequenceReader

FASTA = b">a\nACGTNacgt\n>b\nGGGG\n"


def test_convert_detects_source_and_names_target(tmp_path, capsys):
    src = tmp_path / "in.fa.gz"
    src.write_bytes(gzip.compress(FASTA))
    dst = tmp_path / "out.fa.gz"
    assert main(["convert", str(src), str(dst), "-t", "2", "-m", "64M"]) == 0
    out, err = capsys.readouterr()
    assert out.splitlines() == [str(dst), f"{dst}.fai", f"{dst}.gzi"]
    assert err.startswith("gzip-fasta -> bgzip-fasta: 2 records")
    with SequenceReader(dst) as f:
        assert f.fetch(("a", 0, 9)).text == "ACGTNacgt"


def test_from_and_to_state_formats(tmp_path, capsys):
    src, dst = tmp_path / "genome.txt", tmp_path / "genome.out"
    src.write_bytes(FASTA)
    assert main(["convert", str(src), str(dst), "--to", "2bit", "-q"]) == 0
    out, err = capsys.readouterr()
    assert out.splitlines() == [str(dst)] and err == ""
    with SequenceReader(dst) as f:
        assert f.fetch(("b", 0, 4)).text == "GGGG"
    (tmp_path / "peaks").write_text("c\t0\t5\n")
    assert main(["convert", str(tmp_path / "peaks"), str(tmp_path / "p.bb"), "--from", "bed"]) == 0


def test_dry_run_writes_nothing(tmp_path, capsys):
    src, dst = tmp_path / "a.fa", tmp_path / "a.2bit"
    src.write_bytes(FASTA)
    assert main(["convert", str(src), str(dst), "--dry-run", "-t", "3", "-m", "1G"]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == f"source  {src}  (fasta)" and out[1] == f"target  {dst}  (2bit)"
    assert out[2].startswith("memory  fixed ") and "(budget 1.0 GiB, 3 threads)" in out[2]
    assert not dst.exists()


def test_bam_sort(tmp_path, capsys):
    sam = tmp_path / "r.sam"
    sam.write_text(
        "@SQ\tSN:c\tLN:1000\n"
        "r1\t0\tc\t500\t9\t4M\t*\t0\t0\tACGT\t*\n"
        "r2\t16\tc\t10\t9\t4M\t*\t0\t0\tACGT\t*\n"
    )
    unsorted, out = tmp_path / "r.bam", tmp_path / "s.bam"
    assert main(["convert", str(sam), str(unsorted), "-q"]) == 0
    capsys.readouterr()
    assert main(["convert", str(unsorted), str(out), "--index"]) == 1
    assert "pass sort" in capsys.readouterr().err
    assert main(["convert", str(unsorted), str(out), "--dry-run"]) == 1
    assert main(["convert", str(unsorted), str(out), "--sort", "--dry-run"]) == 0
    assert not out.exists()
    capsys.readouterr()
    assert main(["convert", str(unsorted), str(out), "--sort", "--index", "-q"]) == 0
    assert capsys.readouterr().out.splitlines() == [str(out), f"{out}.bai"]


def test_failures_and_usage(tmp_path, capsys):
    src = tmp_path / "a.fa"
    src.write_bytes(FASTA)
    assert main(["convert", str(src), str(tmp_path / "a.xyz")]) == 1
    assert "ow: error:" in capsys.readouterr().err
    assert main(["convert", str(tmp_path / "none.fa"), str(tmp_path / "a.2bit")]) == 1
    assert "ow: error:" in capsys.readouterr().err
    for argv in (
        ["convert", str(src)],
        ["convert", "a", "b", "--to", "nope"],
        ["convert", "a", "b", "-m", "lots"],
        [],
    ):
        with pytest.raises(SystemExit) as e:
            main(argv)
        assert e.value.code == 2
    capsys.readouterr()


def test_help_lists_conversions(capsys):
    with pytest.raises(SystemExit):
        main(["convert", "--help"])
    out = capsys.readouterr().out
    assert "conversions (source formats -> target format):" in out
    assert "vcf, gzip-vcf" in out and "->  bcf" in out


def test_ow_executable(tmp_path):
    ow = shutil.which("ow", path=str(Path(sys.executable).parent))
    if ow is None:
        pytest.skip("the ow executable is not installed beside this interpreter")
    src, dst = tmp_path / "a.fa", tmp_path / "a.2bit"
    src.write_bytes(FASTA)
    done = subprocess.run(
        [ow, "convert", str(src), str(dst), "-q"], capture_output=True, text=True, check=True
    )
    assert done.stdout.splitlines() == [str(dst)]


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--threads", "0"),
        ("--threads", "-1"),
        ("--threads", str(1 << 64)),
        ("--memory", "0"),
        ("--memory", "-1G"),
        ("--memory", "-0.5"),
        ("--memory", "inf"),
        ("--memory", "-inf"),
        ("--memory", "nan"),
        ("--memory", "1e1000G"),
        ("--memory", str(1 << 64)),
        ("--memory", "lots"),
    ],
)
def test_invalid_budgets_are_usage_errors(tmp_path, capsys, option, value):
    target = tmp_path / "out.2bit"
    with pytest.raises(SystemExit) as e:
        main(["convert", str(tmp_path / "missing.fa"), str(target), f"{option}={value}"])
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert f"argument {option}" in err or f"/{option}" in err
    assert "Traceback" not in err
    assert not target.exists()


@pytest.mark.parametrize("cores", [None, 3])
def test_default_threads_use_available_cores(tmp_path, capsys, monkeypatch, cores):
    monkeypatch.setattr("omicsweave.cli._convert.os.cpu_count", lambda: cores)
    src, dst = tmp_path / "a.fa", tmp_path / "a.2bit"
    src.write_bytes(FASTA)
    assert main(["convert", str(src), str(dst), "--dry-run"]) == 0
    assert f"{cores or 1} threads)" in capsys.readouterr().out
    assert not dst.exists()
