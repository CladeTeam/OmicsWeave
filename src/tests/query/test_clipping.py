import numpy as np
import pytest

from omicsweave.formats import BigBedWriter, BigWigWriter
from omicsweave.query import IntervalReader, StrandedTrackReader, TrackReader

REGION = ("c", 5, 15)
STORED = [(0, 10), (12, 30)]
CLIPPED = [(5, 10), (12, 15)]


def bounds(result):
    return list(zip(result.start.tolist(), result.end.tolist(), strict=True))


def bigwig(path):
    with BigWigWriter(path, {"c": 100}) as writer:
        writer.add("c", np.array([0, 12]), np.array([10, 30]), np.array([1.0, 2.0]))
    return path


def test_bigbed_features_keep_stored_bounds_unless_clipped(tmp_path):
    path = tmp_path / "a.bb"
    with BigBedWriter(path, {"c": 100}) as writer:
        writer.add("c", np.array([0, 12]), np.array([10, 30]), ["a", "b"])
    with IntervalReader(path) as reader:
        whole = reader.fetch(REGION)
        cut = reader.fetch(REGION, clip=True)
    assert bounds(whole) == STORED
    assert bounds(cut) == CLIPPED
    assert whole.annotation.to_list() == cut.annotation.to_list() == ["a", "b"]


def test_bigwig_runs_are_clipped_unless_asked_for_stored_bounds(tmp_path):
    with TrackReader(bigwig(tmp_path / "a.bw")) as reader:
        cut = reader.intervals(REGION)
        whole = reader.intervals(REGION, clip=False)
    assert bounds(cut) == CLIPPED
    assert bounds(whole) == STORED
    assert cut.value.tolist() == whole.value.tolist() == [1.0, 2.0]


@pytest.mark.parametrize(("clip", "expected"), [(True, CLIPPED), (False, STORED)])
def test_stranded_runs_pass_clip_through(tmp_path, clip, expected):
    plus = bigwig(tmp_path / "plus.bw")
    minus = bigwig(tmp_path / "minus.bw")
    with StrandedTrackReader(plus, minus) as reader:
        assert bounds(reader.intervals(("c", 5, 15, "-"), clip=clip)) == expected
