import pickle
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from omicsweave.cache import BlockCache
from omicsweave.formats import BbiOptions, BigBedSchema, BigBedWriter
from omicsweave.query import MsaReader


@pytest.fixture
def alignment_file(tmp_path):
    path = tmp_path / "alignment.bigmaf"
    schema = BigBedSchema(
        4,
        auto_sql='table bigMaf\n"alignment"\n(\nstring chrom;\nuint chromStart;\n'
        "uint chromEnd;\nstring mafBlock;\n)\n",
    )
    with BigBedWriter(path, {"chr1": 20}, schema, options=BbiOptions(zoom_levels=0)) as writer:
        writer.add(
            "chr1",
            [2, 12],
            [8, 16],
            [
                "a;s ref.chr1 2 6 + 20 Aa-RYun;s other.chr2 0 6 + 30 cGATNn-;",
                "a;s ref.chr1 12 4 + 20 ACGT;e other.chr2 0 4 + 30 I;",
            ],
        )
    return path


def test_cached_clones_have_independent_staging_and_owned_results(alignment_file):
    cache = BlockCache(max_bytes=1 << 20)
    reader = MsaReader(alignment_file, ["ref", "other"], cache=cache)
    result = reader.fetch(("chr1", 0, 24))
    expected = result.bases.copy()
    assert result.bases.flags.owndata
    assert result.bases[0].tobytes() == b"--AaRYun----ACGT--------"
    assert result.bases[1].tobytes() == b"--cGTNn-----------------"
    loads = cache.stats()["loads"]
    assert loads > 0
    clones = [reader.clone() for _ in range(3)]
    reader.close()
    try:
        with ThreadPoolExecutor(max_workers=3) as pool:
            values = list(pool.map(lambda r: r.fetch(("chr1", 0, 24)).bases, clones))
        for value in values:
            np.testing.assert_array_equal(value, expected)
        assert cache.stats()["loads"] == loads
        cache.clear()
        clones[0].fetch(("chr1", 12, 16))
        clones[0].set_cache(None)
        loads = cache.stats()["loads"]
        clones[0].fetch(("chr1", 2, 8))
        assert cache.stats()["loads"] == loads
    finally:
        for clone in clones:
            clone.close()
    np.testing.assert_array_equal(result.bases, expected)


def test_reverse_complements_bases_and_preserves_missing_fill(alignment_file):
    with MsaReader(alignment_file, ["ref", "absent"], fill="A") as reader:
        forward = reader.fetch(("chr1", 0, 10)).bases
        reverse = reader.fetch(("chr1", 0, 10, "-")).bases
        assert forward[0].tobytes() == b"AAAaRYunAA"
        assert reverse[0].tobytes() == b"AAnaRYtTAA"
        assert reverse[1].tobytes() == b"A" * 10


def test_output_validation_and_empty_queries(alignment_file):
    with MsaReader(alignment_file, ["ref", "other"]) as reader:
        for out in (
            np.empty((2, 6), dtype=np.int8),
            np.empty((2, 12), dtype=np.uint8)[:, ::2],
        ):
            with pytest.raises(ValueError, match="uint8"):
                reader.fetch(("chr1", 2, 8), out=out)
        output = np.empty((2, 6), dtype=np.uint8)
        result = reader.fetch(("chr1", 2, 8), out=output)
        assert result.bases is output
        assert reader.fetch(("chr1", 20, 20)).bases.shape == (2, 0)
        restored = pickle.loads(pickle.dumps(reader))
        try:
            np.testing.assert_array_equal(restored.fetch(("chr1", 2, 8)).bases, output)
        finally:
            restored.close()
    with MsaReader(alignment_file) as reader:
        result = reader.fetch(("chr1", 20, 24))
        assert result.species == ()
        assert result.bases.shape == (0, 4)


def test_native_output_shape_and_numpy_allocation_are_atomic(alignment_file, monkeypatch):
    with MsaReader(alignment_file, ["ref", "other"]) as reader:
        output = np.zeros((2, 6), dtype=np.uint8)
        for wrong in (output.reshape(3, 4), output.reshape(1, 12)):
            with pytest.raises(ValueError, match="shape"):
                reader._native.fetch_into("chr1", 2, 8, False, wrong)
        assert not output.any()
        expected = reader.fetch(("chr1", 2, 8)).bases.copy()

        def replaced_empty(*args, **kwargs):
            raise AssertionError("reader allocation called Python numpy.empty")

        with monkeypatch.context() as patch:
            patch.setattr(np, "empty", replaced_empty)
            result = reader.fetch(("chr1", 2, 8))
        assert result.bases.flags.owndata
        np.testing.assert_array_equal(result.bases, expected)


def test_native_readonly_output_is_rejected_without_changing_storage(alignment_file):
    with MsaReader(alignment_file, ["ref", "other"]) as reader:
        output = np.full((2, 6), 17, dtype=np.uint8)
        before = output.copy()
        output.flags.writeable = False
        with pytest.raises(ValueError, match="out cannot be written"):
            reader._native.fetch_into("chr1", 2, 8, False, output)
        np.testing.assert_array_equal(output, before)
        output.flags.writeable = True
        reader._native.fetch_into("chr1", 2, 8, False, output)
        assert output[0].tobytes() == b"AaRYun"
        assert output[1].tobytes() == b"cGTNn-"


def test_dynamic_queries_share_one_atomic_context(alignment_file):
    with MsaReader(alignment_file) as reader:
        regions = [("chr1", 2, 8), ("chr1", 12, 16)] * 16
        with ThreadPoolExecutor(max_workers=4) as pool:
            outputs = list(pool.map(reader.fetch, regions))
        for region, result in zip(regions, outputs, strict=True):
            assert result.bases.flags.owndata
            if region[1] == 2:
                assert result.species == ("ref", "other")
                assert result.bases.shape == (2, 6)
                assert result.bases[0].tobytes() == b"AaRYun"
                assert result.bases[1].tobytes() == b"cGTNn-"
            else:
                assert result.species == ("ref",)
                assert result.bases.shape == (1, 4)
                assert result.bases[0].tobytes() == b"ACGT"
