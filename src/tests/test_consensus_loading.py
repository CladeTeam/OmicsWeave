"""Only selectors in one explicit batch share a loaded consensus region."""

import sys
from pathlib import Path

import pytest

import omicsweave as ow


def make_reader(tmp_path):
    reference = tmp_path / "reference.fa"
    with ow.formats.FastaWriter(reference) as writer:
        writer.add("c", b"ACGTACGT")
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c,length=8>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    bcf = tmp_path / "empty.bcf"
    ow.convert.run(vcf, bcf, index=True, threads=1)
    ref = ow.query.SequenceReader(reference)
    return ref, ow.query.ConsensusReader(bcf, ref)


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux per-process IO counters")
def test_consensus_reloads_between_calls_but_shares_one_selector_batch(tmp_path):
    reference, reader = make_reader(tmp_path)
    region = ("c", 0, 8)

    def reads():
        counters = dict(line.split(": ") for line in Path("/proc/self/io").read_text().splitlines())
        return int(counters["syscr"])

    try:
        assert reader.fetch(region).text == "ACGTACGT"
        first = reads()
        for _ in range(3):
            assert reader.fetch(region).text == "ACGTACGT"
        singles = reads() - first
        first = reads()
        batch = reader.fetch_batch(region, [ow.types.AllAlt()] * 3)
        together = reads() - first
        assert batch.offsets.tolist() == [0, 8, 16, 24]
        assert batch.bases.tobytes() == b"ACGTACGT" * 3
        # Both measurements include the same /proc counter read. Each
        # additional public call must issue a real reference read.
        assert singles - together >= 2, (singles, together)
    finally:
        reader.close()
        reference.close()


def test_consensus_does_not_return_old_window_after_reference_close(tmp_path):
    reference, reader = make_reader(tmp_path)
    try:
        assert reader.fetch(("c", 0, 8)).text == "ACGTACGT"
        reference.close()
        with pytest.raises(ValueError, match="closed"):
            reader.fetch(("c", 0, 8))
    finally:
        reader.close()
        reference.close()
