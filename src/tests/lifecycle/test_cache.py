import multiprocessing
import os
import pickle
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import omicsweave as ow
from omicsweave.cache import BlockCache
from omicsweave.formats import (
    FastaWriter,
)
from omicsweave.query import (
    StrandedTrackReader,
    TrackReader,
)

from ..support.tracks import track


def reference(tmp_path, suffix=".fa.gz"):
    path = tmp_path / ("ref" + suffix)
    bases = b"ACGTN" * 40000
    with FastaWriter(path, line_width=61) as writer:
        writer.add("a", bases)
    return path, bases


def test_cached_reads_clones_and_independent_readers(tmp_path):
    path, bases = reference(tmp_path)
    cache = BlockCache(max_bytes=1 << 20, max_entries=64)
    reader = ow.query.SequenceReader(path, cache=cache)
    region = ow.Region("a", 120, 300)
    result = reader.fetch(region)
    assert result.bases.tobytes() == bases[120:300]
    loads = cache.stats()["loads"]
    assert loads > 0 and cache.stats()["resident_entries"] > 0
    assert reader.fetch(region).bases.tobytes() == bases[120:300]
    assert cache.stats()["loads"] == loads
    clone = reader.clone()
    reader.close()
    assert clone.fetch(region).bases.tobytes() == bases[120:300]
    assert cache.stats()["loads"] == loads
    # A reader opened on its own is a separate source session.
    independent = ow.query.SequenceReader(path, cache=cache)
    independent.fetch(region)
    assert cache.stats()["loads"] > loads
    cache.clear()
    assert cache.stats()["live_bytes"] == cache.stats()["resident_entries"] == 0
    assert clone.fetch(region).bases.tobytes() == result.bases.tobytes()
    clone.close()
    independent.close()


def test_zero_budget_bypasses_cache_and_explicit_detach(tmp_path):
    path, bases = reference(tmp_path)
    cache = BlockCache(max_bytes=0, max_entries=0)
    with ow.query.SequenceReader(path, cache=cache) as reader:
        assert reader.fetch(("a", 20, 30)).bases.tobytes() == bases[20:30]
        assert cache.stats()["admission_rejections"] > 0
        assert cache.stats()["resident_entries"] == 0
        cache.resize(max_bytes=1 << 20, max_entries=64)
        reader.fetch(("a", 20, 30))
        assert cache.stats()["resident_entries"] > 0
        assert (cache.stats()["max_bytes"], cache.stats()["max_entries"]) == (1 << 20, 64)
        reader.set_cache(None)
        loads = cache.stats()["loads"]
        out = np.empty(10, np.uint8)
        assert reader.fetch(("a", 20, 30), out=out).bases is out
        assert out.tobytes() == bases[20:30]
        assert cache.stats()["loads"] == loads


def test_cloned_readers_share_one_concurrent_block_load(tmp_path):
    path, bases = reference(tmp_path)
    cache = BlockCache(max_bytes=1 << 20)
    original = ow.query.SequenceReader(path, cache=cache)
    readers = [original.clone() for _ in range(4)]
    original.close()
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            values = list(pool.map(lambda r: r.fetch(("a", 7, 127)).bases.tobytes(), readers))
        assert values == [bases[7:127]] * 4
        assert cache.stats()["loads"] == 1
    finally:
        for reader in readers:
            reader.close()


def test_cache_configuration_and_plain_fasta(tmp_path):
    path, bases = reference(tmp_path, ".fa")
    cache = BlockCache(max_bytes=4096)
    with pytest.raises(ow.errors.UnsupportedError):
        ow.query.SequenceReader(path, cache=cache)
    with ow.query.SequenceReader(path) as reader:
        with pytest.raises(ow.errors.UnsupportedError):
            reader.set_cache(cache)
        assert reader.fetch(("a", 2, 12)).bases.tobytes() == bases[2:12]
    with pytest.raises(TypeError):
        BlockCache(max_bytes=True)
    with pytest.raises(ValueError):
        cache.resize(max_bytes=-1, max_entries=4)
    with pytest.raises(TypeError):
        ow.query.SequenceReader(path, cache={})
    with pytest.raises(TypeError):
        pickle.dumps(cache)


def test_pickle_reopens_without_the_cache(tmp_path):
    path, bases = reference(tmp_path)
    cache = BlockCache(max_bytes=1 << 20)
    with ow.query.SequenceReader(path, cache=cache) as reader:
        reader.fetch(("a", 1, 120))
        restored = pickle.loads(pickle.dumps(reader))
        try:
            loads = cache.stats()["loads"]
            assert restored.fetch(("a", 1, 120)).bases.tobytes() == bases[1:120]
            assert cache.stats()["loads"] == loads
        finally:
            restored.close()


@pytest.mark.parametrize("kind", ["track", "features", "contact", "2bit"])
def test_format_blocks_are_shared_by_clones_and_cleared(tmp_path, kind):
    chroms = {"a": 1000}
    cache = BlockCache(max_bytes=1 << 20)
    region = ow.Region("a", 10, 30)
    if kind == "track":
        path = tmp_path / "x.bw"
        with ow.formats.BigWigWriter(path, chroms) as writer:
            writer.add("a", [0, 20], [20, 40], [1.5, 3.0])
        reader = ow.query.TrackReader(path, cache=cache)
        query = region

        def fetch(r):
            return r.signal(query).values.copy()
    elif kind == "features":
        path = tmp_path / "x.bb"
        with ow.formats.BigBedWriter(path, chroms, ow.formats.BigBedSchema(4)) as writer:
            writer.add("a", [0, 20], [20, 40], ["one", "two"])
        reader = ow.query.IntervalReader(path, cache=cache)
        query = region

        def fetch(r):
            return r.fetch(query).start.copy()
    elif kind == "contact":
        path = tmp_path / "x.hic"
        header = ow.formats.HicHeader("test", chroms, (10,))
        with ow.formats.HicWriter(path, header) as writer:
            writer.write_matrix("a", "a", 10, [1, 2], [1, 2], [1.5, 3.0])
        reader = ow.query.ContactReader(path, cache=cache)
        query = ow.types.ContactWindow(region, region, 10)

        def fetch(r):
            return r.matrix(query).values.copy()
    else:
        reference_path, _ = reference(tmp_path, ".fa")
        path = tmp_path / "x.2bit"
        ow.convert.run(
            ow.formats.Fasta(reference_path), ow.formats.TwoBit(path), threads=1, memory="64M"
        )
        reader = ow.query.SequenceReader(path, cache=cache)
        query = region

        def fetch(r):
            return r.fetch(query).bases.copy()

    try:
        expected = fetch(reader)
        assert cache.stats()["resident_entries"] > 0
        loads = cache.stats()["loads"]
        clone = reader.clone()
        try:
            np.testing.assert_array_equal(fetch(clone), expected)
            assert cache.stats()["loads"] == loads
            cache.clear()
            assert cache.stats()["resident_entries"] == 0
            np.testing.assert_array_equal(fetch(clone), expected)
            assert cache.stats()["loads"] > loads
        finally:
            clone.close()
    finally:
        reader.close()


def test_bigbed_clones_share_the_file_and_take_a_block_cache(tmp_path):
    path = tmp_path / "x.bb"
    schema = ow.formats.BigBedSchema(field_count=4)
    with ow.formats.BigBedWriter(path, {"a": 100}, schema) as writer:
        writer.add("a", [1, 5], [10, 12], ["one", "two"])
    reader = ow.query.IntervalReader(path)
    clone = reader.clone()
    reader.close()
    try:
        assert clone.fetch(("a", 7, 8)).annotation.to_list() == ["one", "two"]
        clone.set_cache(BlockCache(max_bytes=1 << 20))
        assert clone.fetch(("a", 7, 8)).annotation.to_list() == ["one", "two"]
    finally:
        clone.close()


def test_default_sequence_arrays_own_their_storage(tmp_path):
    path, bases = reference(tmp_path)
    with ow.query.SequenceReader(path) as reader:
        result = reader.fetch(("a", 0, 12)).bases
        assert result.flags.owndata and result.flags.writeable
        result[0] = ord("T")
        assert reader.fetch(("a", 0, 12)).bases.tobytes() == bases[:12]
        result.resize(30)
        assert result.shape == (30,)


def test_consensus_clone_shares_variants_and_independent_reference_context(tmp_path):
    path, _ = reference(tmp_path)
    variants = tmp_path / "v.vcf"
    variants.write_text(
        "##fileformat=VCFv4.3\n##contig=<ID=a,length=200000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "a\t1\t.\tA\tC\t.\t.\t.\n"
    )
    ow.convert.run(variants, tmp_path / "v.bcf", index=True, threads=1)
    cache = BlockCache(max_bytes=1 << 20)
    reference_reader = ow.query.SequenceReader(path, cache=cache)
    original = ow.query.ConsensusReader(tmp_path / "v.bcf", reference_reader)
    clone = original.clone()
    original.close()
    reference_reader.close()
    try:
        assert clone.fetch(("a", 0, 5)).bases.tobytes() == b"CCGTN"
        assert cache.stats()["loads"] > 0
    finally:
        clone.close()


def _spawn_fetch(reader, result):
    try:
        result.send(reader.fetch(("a", 0, 5)).bases.tobytes())
    except BaseException as error:
        result.send(repr(error))
    finally:
        result.close()


def test_spawn_reopens_reader_without_serializing_cache(tmp_path):
    path, _ = reference(tmp_path)
    cache = BlockCache(max_bytes=1 << 20)
    with ow.query.SequenceReader(path, cache=cache) as reader:
        reader.fetch(("a", 0, 5))
        before = cache.stats()
        context = multiprocessing.get_context("spawn")
        parent, child = context.Pipe(duplex=False)
        process = context.Process(target=_spawn_fetch, args=(reader, child))
        process.start()
        child.close()
        try:
            assert parent.poll(30)
            value = parent.recv()
            assert value == b"ACGTN", value
        finally:
            process.join(30)
            parent.close()
        assert process.exitcode == 0
        assert cache.stats() == before


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork required")
def test_inherited_reader_and_cache_reject_before_touching_locks(tmp_path):
    path, _ = reference(tmp_path)
    cache = BlockCache(max_bytes=1 << 20)
    with ow.query.SequenceReader(path, cache=cache) as reader:
        reader.fetch(("a", 0, 10))
        read_fd, write_fd = os.pipe()
        pid = os.fork()
        if pid == 0:
            os.close(read_fd)
            failures = 0
            for operation in [lambda: reader.fetch(("a", 0, 10)), cache.stats]:
                try:
                    operation()
                except RuntimeError:
                    failures += 1
            fresh_cache = BlockCache(max_bytes=1 << 20)
            with ow.query.SequenceReader(path, cache=fresh_cache) as fresh:
                if fresh.fetch(("a", 0, 5)).bases.tobytes() == b"ACGTN":
                    failures += 1
            os.write(write_fd, str(failures).encode())
            os._exit(0)
        os.close(write_fd)
        try:
            assert os.read(read_fd, 16) == b"3"
            _, status = os.waitpid(pid, 0)
            assert os.waitstatus_to_exitcode(status) == 0
        finally:
            os.close(read_fd)


def test_stranded_tracks_share_a_cache_and_clone(tmp_path):
    plus = track(tmp_path / "p.bw", np.random.default_rng(6))
    minus = track(tmp_path / "m.bw", np.random.default_rng(7))
    cache = BlockCache(max_bytes=1 << 24)
    with StrandedTrackReader(plus, minus, cache=cache) as rna:
        rna.signal(("chr1", 0, 10_000, "-"))
        loads = cache.stats()["loads"]
        assert loads > 0
        clone = rna.clone()
        np.testing.assert_array_equal(
            clone.signal(("chr1", 0, 10_000, "-")).values,
            rna.signal(("chr1", 0, 10_000, "-")).values,
        )
        assert cache.stats()["loads"] == loads
        clone.close()
        assert not rna.closed
        with TrackReader(plus) as t:
            assert rna.plus.summary == t.summary
