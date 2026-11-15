import multiprocessing

import numpy as np
import pytest

import omicsweave as ow


@pytest.fixture(params=["alignment", "variant"])
def indexed_file(tmp_path, request):
    if request.param == "alignment":
        source, target = tmp_path / "in.sam", tmp_path / "in.bam"
        source.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n"
            "a\t0\tchr1\t11\t60\t3M1I2M\t*\t0\t0\tACGTAC\tIIIIII\tNM:i:1\n"
            "b\t16\tchr1\t31\t20\t2S4M3N2M\t*\t0\t0\tGGACGTAC\t*\n"
            "c\t0\tchr2\t5\t30\t4M\t*\t0\t0\tTTAA\tHHHH\n"
        )
        ow.convert.run(source, target, index="bai", threads=1, memory="64M")
        reader_type = ow.query.AlignmentReader
    else:
        source, target = tmp_path / "in.vcf", tmp_path / "in.bcf"
        source.write_text(
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n"
            "##contig=<ID=chr2,length=500>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfirst\tsecond\tthird\n"
            "chr1\t11\tv1\tA\tC,G\t50\tPASS\t.\tGT\t1|0\t2\t./.\n"
            "chr1\t31\t.\tT\tTAA\t50\tPASS\t.\tGT\t./1\t0|0\t1/1\n"
            "chr1\t70\t.\tA\tG\t50\tPASS\t.\tDP\t4\t5\t6\n"
            "chr2\t5\tv2\tG\tA\t50\tPASS\t.\tGT\t1\t0\t.\n"
        )
        ow.convert.run(source, target, index=True, threads=1, memory="64M")
        reader_type = ow.query.VariantReader
    return request.param, target, reader_type


def payload(reader, kind):
    region = ("chr1", 0, 100)
    if kind == "alignment":
        result = reader.fetch(region, "all")
        return (
            result.start.tobytes(),
            result.end.tobytes(),
            result.flag.tobytes(),
            result.mapq.tobytes(),
            result.names.offsets.tobytes(),
            result.names.data.tobytes(),
            result.seqs.data.tobytes(),
            result.quals.data.tobytes(),
            result.cigars.offsets.tobytes(),
            result.cigars.ops.tobytes(),
            result.aux.data.tobytes(),
        )
    result = reader.genotypes(region)
    return (
        result.start.tobytes(),
        result.end.tobytes(),
        result.gt.tobytes(),
        result.phased.tobytes(),
        result.ploidy,
        result.samples,
    )


def test_indexed_prefetch_clone_and_eviction(indexed_file):
    kind, path, reader_type = indexed_file
    cache = ow.cache.BlockCache(max_bytes=1 << 20)
    reader = reader_type(path, cache=cache)
    region = ("chr1", 0, 100)
    try:
        dependencies = reader.describe_cache(region)
        assert dependencies and all(k.source == reader.cache_namespace for k in dependencies)
        reader.prefetch_cache(region)
        assert set(dependencies).issubset(cache.keys())
        reader.reset_stream_stats()
        expected = payload(reader, kind)
        stats = reader.stream_stats()
        assert stats["inflate_calls"] == stats["read_calls"] == 0
        clone = reader.clone()
        assert clone.cache_namespace == reader.cache_namespace
    finally:
        reader.close()
    try:
        assert payload(clone, kind) == expected
        assert clone.stream_stats()["inflate_calls"] == clone.stream_stats()["read_calls"] == 0
        cache.clear(source=clone.cache_namespace)
        clone.reset_stream_stats()
        assert payload(clone, kind) == expected
        assert clone.stream_stats()["inflate_calls"] > 0
        clone.set_cache(None)
        clone.reset_stream_stats()
        assert payload(clone, kind) == expected
        work = clone.stream_stats()["inflate_calls"]
        assert payload(clone, kind) == expected
        assert clone.stream_stats()["inflate_calls"] > work
    finally:
        clone.close()


def test_zero_budget_indexed_queries_do_not_change_output(indexed_file):
    kind, path, reader_type = indexed_file
    cache = ow.cache.BlockCache(max_bytes=0, max_entries=0)
    with reader_type(path) as baseline, reader_type(path, cache=cache) as cached:
        assert payload(cached, kind) == payload(baseline, kind)
        assert cache.stats()["resident_entries"] == 0
        assert not cache.keys()
        assert cached.stream_stats()["inflate_calls"] > 0


def test_alignment_fields_and_genotype_selection_survive_cache_reuse(indexed_file):
    kind, path, reader_type = indexed_file
    cache = ow.cache.BlockCache(max_bytes=1 << 20)
    with reader_type(path, cache=cache) as reader:
        region = ("chr1", 0, 100)
        if kind == "alignment":
            full = reader.fetch(region, "all")
            assert full.names.to_list() == ["a", "b"]
            assert full.seqs.to_list() == ["ACGTAC", "GGACGTAC"]
            assert full.cigars.string(0) == "3M1I2M"
            bare = reader.fetch(region)
            assert bare.names is bare.seqs is bare.cigars is bare.aux is None
            partial = reader.fetch(region, ("seq", "qual"))
            assert partial.names is partial.cigars is partial.aux is None
            np.testing.assert_array_equal(partial.start, full.start)
            assert partial.seqs.data.tobytes() == full.seqs.data.tobytes()
        else:
            all_samples = reader.genotypes(region)
            explicit = reader.genotypes(region, reader.samples)
            np.testing.assert_array_equal(explicit.gt, all_samples.gt)
            np.testing.assert_array_equal(explicit.phased, all_samples.phased)
            selection = ["second", "first", "second"]
            reordered = reader.genotypes(region, selection)
            assert reordered.samples == tuple(selection)
            np.testing.assert_array_equal(reordered.gt, all_samples.gt[:, [1, 0, 1], :])
            np.testing.assert_array_equal(reordered.phased, all_samples.phased[:, [1, 0, 1]])
            haploid = reader.genotypes(("chr2", 0, 20), "second")
            assert haploid.ploidy == 1 and haploid.gt.tolist() == [[[0]]]
            empty_samples = reader.genotypes(region, [])
            assert empty_samples.gt.shape[:2] == (3, 0)
            assert len(reader.genotypes(("chr1", 500, 600))) == 0
            with pytest.raises(KeyError, match="no sample"):
                reader.genotypes(region, "missing")


def _query_in_process(reader, kind, pipe):
    try:
        cache_absent = False
        try:
            reader.prefetch_cache(("chr1", 0, 100))
        except ow.errors.OmicsWeaveError:
            cache_absent = True
        pipe.send(("ok", payload(reader, kind), cache_absent))
    except BaseException as error:
        pipe.send(("error", type(error).__name__))
    finally:
        pipe.close()


@pytest.mark.parametrize("method", ["spawn", "fork"])
def test_indexed_readers_spawn_reopen_and_fork_reject(indexed_file, method):
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"{method} start method unavailable")
    kind, path, reader_type = indexed_file
    cache = ow.cache.BlockCache(max_bytes=1 << 20)
    with reader_type(path, cache=cache) as reader:
        expected = payload(reader, kind)
        context = multiprocessing.get_context(method)
        parent, child = context.Pipe(duplex=False)
        process = context.Process(target=_query_in_process, args=(reader, kind, child))
        process.start()
        child.close()
        try:
            assert parent.poll(30), "child reader did not finish"
            received = parent.recv()
            if method == "spawn":
                assert received == ("ok", expected, True), received
            else:
                assert received == ("error", "RuntimeError"), received
        finally:
            process.join(5)
            if process.is_alive():
                process.terminate()
                process.join(5)
            parent.close()
        assert process.exitcode == 0
        assert payload(reader, kind) == expected
