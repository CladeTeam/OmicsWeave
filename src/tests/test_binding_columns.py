import numpy as np
import pytest

from omicsweave import _native, errors
from omicsweave.formats import (
    BbiOptions,
    BigBedSchema,
    BigBedWriter,
    BigWigWriter,
    FastaWriter,
    HicHeader,
    HicWriter,
)
from omicsweave.query import ConsensusReader, ContactReader, IntervalReader, TrackReader
from omicsweave.types import AllAlt, ContactWindow, Haplotype


def test_interval_and_track_columns_keep_order_bits_and_ownership(tmp_path):
    n = 4097
    starts = np.arange(n, dtype=np.uint32) * 3
    ends = starts + 2
    values = np.resize(np.array([0.0, -0.0, 1.5, -7.25], np.float32), n)
    notes = [f"note-{i}" if i % 3 else "" for i in range(n)]
    chroms = {"c": int(ends[-1]) + 1}
    options = BbiOptions(zoom_levels=0)
    track_path = tmp_path / "columns.bw"
    bed_path = tmp_path / "columns.bb"
    text_path = tmp_path / "columns.bed"
    with BigWigWriter(track_path, chroms, options) as writer:
        writer.add("c", starts, ends, values)
    with BigBedWriter(bed_path, chroms, BigBedSchema(field_count=4), options) as writer:
        writer.add("c", starts, ends, notes)
    text_path.write_text(
        "".join(f"c\t{s}\t{e}\t{note}\n" for s, e, note in zip(starts, ends, notes, strict=True))
    )
    with TrackReader(track_path) as track:
        full = track.intervals(("c", 0, chroms["c"]))
        np.testing.assert_array_equal(full.start, starts)
        np.testing.assert_array_equal(full.end, ends)
        np.testing.assert_array_equal(full.value.view(np.uint32), values.view(np.uint32))
        assert len(track.intervals(("c", 0, 0)).start) == 0
        assert len(track.intervals(("c", 0, 1)).start) == 1
        again = track.intervals(("c", 0, chroms["c"]))
        assert not np.shares_memory(full.start, again.start)
        assert not np.shares_memory(full.value, again.value)
        full.value[:] = 99
        np.testing.assert_array_equal(again.value.view(np.uint32), values.view(np.uint32))
    for path in (bed_path, text_path):
        with IntervalReader(path) as reader:
            full = reader.fetch(("c", 0, chroms["c"]))
            np.testing.assert_array_equal(full.start, starts)
            np.testing.assert_array_equal(full.end, ends)
            assert full.annotation.to_list() == notes
            assert len(reader.fetch(("c", 0, 0))) == 0
            assert len(reader.fetch(("c", 0, 1))) == 1
            again = reader.fetch(("c", 0, chroms["c"]))
            assert not np.shares_memory(full.start, again.start)
            assert not np.shares_memory(full.annotation.data, again.annotation.data)
            full.start[:] = 99
            np.testing.assert_array_equal(again.start, starts)


def test_contact_columns_and_matrix_shape_validation_leave_outputs_intact(tmp_path):
    n = 4097
    x = np.arange(n, dtype=np.uint32)
    y = np.zeros(n, dtype=np.uint32)
    values = (np.arange(n, dtype=np.uint32) % 17 + 1).astype(np.float32)
    path = tmp_path / "columns.hic"
    header = HicHeader("columns", {"a": n * 10, "b": 30}, (10,))
    with HicWriter(path, header) as writer:
        writer.write_matrix("a", "b", 10, x, y, values)
    with ContactReader(path) as reader:
        window = ContactWindow(("a", 0, n * 10), ("b", 0, 30), 10)
        cells = reader.pixels(window)
        np.testing.assert_array_equal(cells.row, x * 10)
        np.testing.assert_array_equal(cells.col, y)
        np.testing.assert_array_equal(cells.value, values)
        assert len(reader.pixels(ContactWindow(("a", 0, 0), ("b", 0, 30), 10))) == 0
        assert len(reader.pixels(ContactWindow(("a", 0, 10), ("b", 0, 30), 10))) == 1
        again = reader.pixels(window)
        assert not np.shares_memory(cells.row, again.row)
        assert not np.shares_memory(cells.value, again.value)
        cells.value[:] = 99
        np.testing.assert_array_equal(again.value, values)

        small = ContactWindow(("a", 0, 100), ("b", 0, 30), 10)
        wrong = np.full((3, 10), -123, dtype=np.float32)
        with pytest.raises(ValueError, match="out shape"):
            reader.matrix(small, out=wrong)
        np.testing.assert_array_equal(wrong, -123)
        strided = np.full((10, 6), -123, dtype=np.float32)
        with pytest.raises(ValueError, match="contiguous"):
            reader.matrix(small, out=strided[:, ::2])
        np.testing.assert_array_equal(strided, -123)
        column_major = np.full((10, 3), -123, dtype=np.float32, order="F")
        with pytest.raises(ValueError, match="C-contiguous"):
            reader.matrix(small, out=column_major)
        np.testing.assert_array_equal(column_major, -123)
        out = np.full((10, 3), -123, dtype=np.float32)
        assert reader.matrix(small, out=out).values is out
        np.testing.assert_array_equal(out[:, 0], values[:10])
        np.testing.assert_array_equal(out[:, 1:], 0)
        one_column = ContactWindow(("a", 0, 100), ("b", 0, 10), 10)
        both_orders = np.full((10, 1), -123, dtype=np.float32, order="F")
        assert both_orders.flags.c_contiguous and both_orders.flags.f_contiguous
        assert reader.matrix(one_column, out=both_orders).values is both_orders
        np.testing.assert_array_equal(both_orders[:, 0], values[:10])

    native = _native.ContactReader(str(path))
    try:
        tile = ("a", 0, 100, "b", 0, 30, 10, "BP")
        wrong = np.full((3, 10), -123, dtype=np.float32)
        with pytest.raises(ValueError, match="out shape"):
            native.matrix(tile, wrong)
        np.testing.assert_array_equal(wrong, -123)
        strided = np.full((10, 6), -123, dtype=np.float32)
        with pytest.raises(ValueError, match="contiguous"):
            native.matrix(tile, strided[:, ::2])
        np.testing.assert_array_equal(strided, -123)
        column_major = np.full((10, 3), -123, dtype=np.float32, order="F")
        with pytest.raises(ValueError, match="C-contiguous"):
            native.matrix(tile, column_major)
        np.testing.assert_array_equal(column_major, -123)
        out = np.full((10, 3), -123, dtype=np.float32)
        native.matrix(tile, out)
        np.testing.assert_array_equal(out[:, 0], values[:10])
        np.testing.assert_array_equal(out[:, 1:], 0)
        with pytest.raises(errors.InvalidError):
            native.matrix(("a", 0, 100, "b", 0, 30, 0, "BP"), np.empty((10, 3), np.float32))
        with pytest.raises(errors.InvalidRegionError):
            native.matrix(("a", 100, 0, "b", 0, 30, 10, "BP"), np.empty((10, 3), np.float32))
    finally:
        native.close()


def test_cached_consensus_selectors_and_region_changes_recover_after_errors(tmp_path):
    reference = tmp_path / "cache.fa"
    with FastaWriter(reference) as writer:
        writer.add("c", b"ACGTACGT")
    vcf = tmp_path / "cache.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n##contig=<ID=c,length=8>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
        "c\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0|1\n"
    )
    with ConsensusReader(vcf, reference) as builder:
        selectors = [Haplotype("s", "1"), Haplotype("s", "2"), AllAlt()]
        for key, expected in [
            (("c", 0, 4), ["ACGT", "ATGT", "ATGT"]),
            (("c", 1, 4), ["CGT", "TGT", "TGT"]),
            (("c", 0, 4), ["ACGT", "ATGT", "ATGT"]),
        ]:
            assert [builder.fetch(key, selector).text for selector in selectors] == expected
            assert [r.text for r in builder.fetch_batch(key, selectors)] == expected
            with pytest.raises(errors.InvalidError):
                builder.fetch(key, Haplotype("unknown"))
            assert builder.fetch(key, selectors[1]).text == expected[1]
        with pytest.raises(errors.RegionOutOfRangeError):
            builder.fetch(("c", 0, 9))
        assert builder.fetch(("c", 0, 4)).text == "ATGT"
