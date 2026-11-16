from __future__ import annotations

import pickle

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import formats, query, types
from omicsweave.errors import InvalidError

from ._support import DNA
from ._support import collection as collection


def test_ffm_round_trip(collection, tmp_path):
    path, motifs = collection
    with query.MotifReader(path) as r:
        assert r.alphabet == DNA
        assert dict(r.schema) == {"name": "string", "rank": "integer", "score": "float"}
        assert r.ids == tuple(m.id for m in motifs)
        assert "D0" in r and r.index("D0") == 5
        for want, got in zip(motifs, r, strict=True):
            assert got.width == want.width and got.sites == want.sites
            for a, b in zip(want.factors, got.factors, strict=True):
                assert (a.kind, a.positions) == (b.kind, b.positions)
                np.testing.assert_array_equal(a.table, b.table)
        assert r["M1"].annotation == {"name": "one", "rank": 3, "score": 0.5}
        assert r["M0"].annotation == {}
        assert r[-1].matrix is None and r[0].matrix.shape == (r[0].width, 4)
        np.testing.assert_array_equal(r.widths(), [m.width for m in motifs])
        clone = pickle.loads(pickle.dumps(r))
        assert clone.ids == r.ids
        clone.close()
    assert isinstance(ow.open(path), query.MotifReader)

    gapped = types.Motif(
        "G",
        3,
        tuple(types.Factor("joint", (i,), np.array([3.0, 1.0, 0.0, 1.0])) for i in range(3)),
        sites=5.0,
        gaps=types.Gaps(np.array([0.0, 1.0, 0.0]), np.array([2.0, 0.0])),
    )
    with formats.FfmWriter(tmp_path / "g.ffm", DNA) as w:
        w.add(gapped)
    back = query.MotifReader(tmp_path / "g.ffm")[0]
    np.testing.assert_array_equal(back.gaps.insertions, [2.0, 0.0])
    with formats.FfmWriter(tmp_path / "bad.ffm", DNA) as w, pytest.raises(InvalidError):
        w.add(types.Motif("bad", 2, (types.Factor("joint", (0, 0), np.ones(16)),)))
