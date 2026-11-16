import gzip

import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.formats import FastaWriter
from omicsweave.query import ConsensusOptions, ConsensusReader
from omicsweave.types import AllAlt, Haplotype, Record, Variant

REF = "ACGTAACGTCACGTACGTAGAACGTACGTACGTACGTACG"

# Sites-only records, as a variant-effect benchmark lists them: each is
# scored alone, and two sit close enough to share any window.
VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=40>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t5\t.\tA\tG,T\t.\t.\t.
chr1\t10\t.\tC\tCTT\t.\t.\t.
chr1\t20\t.\tGA\tG\t.\t.\t.
"""


def edit(at, ref, alt):
    return REF[:at] + alt + REF[at + len(ref) :]


@pytest.fixture
def reader(tmp_path):
    with FastaWriter(tmp_path / "ref.fa") as w:
        w.add("chr1", REF)
    vcf = tmp_path / "v.vcf.gz"
    vcf.write_bytes(gzip.compress(VCF.encode()))
    bcf = tmp_path / "v.bcf"
    ow.convert.run(vcf, bcf, index=True, threads=1)
    with ConsensusReader(bcf, tmp_path / "ref.fa") as c:
        yield c


def test_records_list_what_record_indexes(reader):
    records = reader.records(("chr1", 0, 40))
    assert records.start.tolist() == [4, 9, 19]
    assert records.end.tolist() == [5, 10, 21]
    assert records.ref.to_list() == ["A", "C", "GA"]
    assert records.alt.to_list() == ["G,T", "CTT", "G"]
    assert records.masked.tolist() == [False, False, False]
    assert len(reader.records(("chr1", 6, 12))) == 1


@pytest.mark.parametrize(
    ("record", "variant", "expected"),
    [
        (Record(0), Variant(4, "A", "G"), edit(4, "A", "G")),
        (Record(0, alt=2), Variant(4, "a", "t"), edit(4, "A", "T")),
        (Record(1), Variant(9, "C", "CTT"), edit(9, "C", "CTT")),
        (Record(2), Variant(19, "GA", "G"), edit(19, "GA", "G")),
    ],
)
def test_one_variant_leaves_the_others_out(reader, record, variant, expected):
    region = ("chr1", 0, 40)
    assert reader.fetch(region, record).text == expected
    assert reader.fetch(region, variant).text == expected
    minus = reader.fetch(("chr1", 0, 40, "-"), record).text
    assert minus == expected.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_batches_mix_genome_and_variant_selections(reader):
    batch = reader.fetch_batch(("chr1", 0, 40), [AllAlt(), Record(0), Record(2)])
    everything = edit(4, "A", "G")
    everything = everything[:9] + "CTT" + everything[10:19] + "G" + everything[21:]
    assert [s.text for s in batch] == [everything, edit(4, "A", "G"), edit(19, "GA", "G")]
    assert batch.alleles[1] == Record(0)


def test_missing_records_are_errors(reader):
    region = ("chr1", 0, 12)
    for selection in [Record(2), Record(0, alt=3), Variant(4, "A", "C"), Variant(19, "GA", "G")]:
        with pytest.raises(errors.InvalidError):
            reader.fetch(region, selection)
    with pytest.raises(ValueError):
        Record(0, alt=0)
    with pytest.raises(ValueError):
        Variant(-1, "A", "G")
    with pytest.raises(TypeError):
        reader.fetch(region, ("record", 0))


def test_a_masked_record_stays_out(tmp_path):
    with FastaWriter(tmp_path / "ref.fa") as w:
        w.add("chr1", REF)
    vcf = tmp_path / "v.vcf.gz"
    vcf.write_bytes(gzip.compress(VCF.encode()))
    bcf = tmp_path / "v.bcf"
    ow.convert.run(vcf, bcf, index=True, threads=1)
    mask = tmp_path / "mask.bed"
    mask.write_text("chr1\t4\t5\n")
    with ConsensusReader(bcf, tmp_path / "ref.fa", options=ConsensusOptions(mask=mask)) as c:
        assert c.records(("chr1", 0, 40)).masked.tolist() == [True, False, False]
        assert c.fetch(("chr1", 0, 40), Record(0)).text == REF[:4] + "N" + REF[5:]


def test_genome_selections_still_apply_every_record(reader):
    assert reader.fetch(("chr1", 0, 6), AllAlt()).text == "ACGTG" + REF[5]
    with pytest.raises(errors.InvalidError):
        reader.fetch(("chr1", 0, 6), Haplotype("nobody"))
