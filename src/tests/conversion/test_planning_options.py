import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.query import (
    TrackReader,
)


def test_conversion_options_are_checked_when_planned(tmp_path):
    (tmp_path / "a.fa").write_bytes(b">a\nACGT\n")
    with pytest.raises(errors.InvalidError, match="index"):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.2bit", index=True)
    with pytest.raises(errors.InvalidError, match="level"):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.2bit", level=3)
    with pytest.raises(errors.InvalidError, match="verify_checksums"):
        ow.convert.plan(tmp_path / "a.fa", tmp_path / "a.fa.gz", verify_checksums=False)
    sam = tmp_path / "r.sam"
    sam.write_text("@SQ\tSN:c\tLN:100\nr\t0\tc\t5\t9\t4M\t*\t0\t0\tACGT\t*\n")
    with pytest.raises(errors.InvalidError, match="min_shift"):
        ow.convert.plan(sam, tmp_path / "r.bam", index="bai", min_shift=12)
    report = ow.convert.run(sam, tmp_path / "r.bam", index="csi", min_shift=12, threads=1)
    assert [f.format for f in report.companions] == ["csi"]
    (tmp_path / "a.bedgraph").write_text("c\t0\t10\t1.5\nc\t10\t20\t2\n")
    small = ow.convert.run(
        tmp_path / "a.bedgraph", tmp_path / "s.bw", bbi=ow.formats.BbiOptions(zoom_levels=0)
    )
    with TrackReader(small.target) as t:
        assert t.zoom_reductions == ()
