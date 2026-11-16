import gzip

import pytest

import omicsweave as ow
from omicsweave import convert
from omicsweave.cli.app import main
from omicsweave.formats import GzipHicProPairs, HicProPairs, detect
from omicsweave.query import ContactReader
from omicsweave.types import ContactWindow

TEXT = (
    "r1\tchr1\t5000\t+\tchr2\t10000\t-\t0\tf_1\tf_1\t0\t0\t\n"
    "r2\tchr2\t10000\t-\tchr1\t5000\t+\t0\tf_1\tf_1\t0\t0\t1\n"
    "r3\tchr1\t4999\t+\tchr2\t5000\t-\t0\tf_1\tf_1\t60\t60\n"
)


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("explicit", [False, True])
def test_hicpro_convert_plan_and_contact_bins(tmp_path, compressed, explicit):
    source = tmp_path / ("contacts.validPairs.gz" if compressed else "contacts.allValidPairs")
    source.write_bytes(gzip.compress(TEXT.encode()) if compressed else TEXT.encode())
    cls = GzipHicProPairs if compressed else HicProPairs
    assert detect(source) == cls(source)
    with pytest.raises(ow.errors.InvalidError, match=r"ow convert .* --chrom-sizes"):
        ow.open(source)
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t600000000\nchr2\t700000000\n")
    target = tmp_path / "contacts.hic"
    plan = convert.plan(
        cls(source) if explicit else source,
        target,
        chrom_sizes=sizes,
        resolutions=[5000, 10000],
        threads=2,
        memory="32M",
    )
    assert plan.source == cls(source)
    report = plan.run()
    assert report.records == 3
    assert report.warnings == ()
    with ContactReader(target) as reader:
        pixels = reader.pixels(ContactWindow(("chr1", 0, 15000), ("chr2", 0, 15000), 5000))
        assert list(zip(pixels.row, pixels.col, pixels.value, strict=True)) == [
            (0, 5000, 1),
            (5000, 10000, 2),
        ]
        assert reader.resolutions("FRAG") == ()
    with pytest.raises(ValueError, match="chromosome sizes"):
        convert.run(source, tmp_path / "missing.hic", threads=1, memory="32M")


@pytest.mark.parametrize("explicit", [False, True])
def test_hicpro_cli_and_missing_chromosome_warning(tmp_path, capsys, explicit):
    source = tmp_path / "contacts.allValidPairs"
    source.write_text(TEXT + "r4\tunknown\t1\t+\tchr2\t2\t-\t0\tf1\tf2\t60\t60\n")
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t600000000\nchr2\t700000000\n")
    target = tmp_path / "contacts.hic"
    assert main(["convert", str(source), str(target), "-t", "1", "-m", "32M"]) == 1
    assert "chromosome sizes" in capsys.readouterr().err
    assert (
        main(
            [
                "convert",
                str(source),
                str(target),
                "--chrom-sizes",
                str(sizes),
                "--resolutions",
                "5000,10000",
                *(["--from", "hicpro-pairs"] if explicit else []),
                "-t",
                "2",
                "-m",
                "32M",
            ]
        )
        == 0
    )
    captured = capsys.readouterr()
    assert str(target) in captured.out
    assert "unknown" in captured.err
    assert "skipped 1 of 4 rows" in captured.err
