import gzip
import pickle

import numpy as np
import pytest

from omicsweave import _native, errors
from omicsweave.formats import Bed, Fasta, FastaWriter, GzipVcf
from omicsweave.query import ConsensusOptions, ConsensusReader, SequenceReader
from omicsweave.types import AllAlt, Haplotype, Region, Samples

REF = "ACGTAACGTCACGTACGTAGAACGTACGTACGTACGTACG"
assert REF[4] == "A" and REF[9] == "C" and REF[19:21] == "GA"

VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=40>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
chr1\t5\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t1|1
chr1\t10\t.\tC\tCTT\t50\tPASS\tDP=3\tGT\t1|0\t0|0
chr1\t20\t.\tGA\tG\t50\tPASS\tDP=20\tGT\t1|1\t0|1
"""


def edit(snv=False, ins=False, dele=False):
    s = list(REF)
    if dele:
        del s[20]
    if ins:
        s[9] = "CTT"
    if snv:
        s[4] = "G"
    return "".join(s)


@pytest.fixture
def files(tmp_path):
    with FastaWriter(tmp_path / "ref.fa") as w:
        w.add("chr1", REF)
    vcf = tmp_path / "v.vcf.gz"
    vcf.write_bytes(gzip.compress(VCF.encode()))
    return vcf, tmp_path / "ref.fa"


def test_haplotypes_and_all_alt(files):
    vcf, ref = files
    c = ConsensusReader(vcf, ref)
    assert c.samples == ("s1", "s2")
    assert (c.records, c.skipped_lines) == (3, 0)
    whole = Region("chr1", 0, 40)
    assert c.fetch(whole).text == edit(snv=True, ins=True, dele=True)
    expect = {
        ("s1", "1"): edit(ins=True, dele=True),
        ("s1", "2"): edit(snv=True, dele=True),
        ("s2", "1"): edit(snv=True),
        ("s2", "2"): edit(snv=True, dele=True),
    }
    for (sample, hap), want in expect.items():
        got = c.fetch(whole, Haplotype(sample, hap))
        assert got.text == want
        assert got.alleles == Haplotype(sample, hap)
    part = c.fetch(("chr1", 2, 12), Haplotype("s1", "1"))
    assert part.text == REF[2:9] + "CTT" + REF[10:12]


def test_batch_equals_single_queries(files):
    vcf, ref = files
    c = ConsensusReader(vcf, ref, ConsensusOptions(chain=True))
    selectors = [Haplotype(s, h) for s in c.samples for h in "12"] + [
        AllAlt(),
        Samples(["s1", "s2"]),
    ]
    batch = c.fetch_batch(("chr1", 0, 40), selectors)
    assert len(batch) == len(selectors)
    assert batch.offsets[0] == 0 and batch.offsets[-1] == len(batch.bases)
    for got, alleles in zip(batch, selectors, strict=True):
        single = c.fetch(("chr1", 0, 40), alleles)
        assert got.text == single.text
        assert got.chain == single.chain and got.chain.startswith("chain")


def test_options(files, tmp_path):
    vcf, ref = files
    shallow = ConsensusReader(vcf, ref, ConsensusOptions(include="DP>5"))
    assert shallow.records == 2
    assert shallow.fetch(("chr1", 0, 40)).text == edit(snv=True, dele=True)
    deep = ConsensusReader(vcf, ref, ConsensusOptions(exclude="DP>5"))
    assert deep.fetch(("chr1", 0, 40)).text == edit(ins=True)
    iupac = ConsensusReader(vcf, ref, ConsensusOptions(iupac_codes=True))
    assert iupac.fetch(("chr1", 0, 10), Samples(["s1"])).text[4] == "R"
    marked = ConsensusReader(vcf, ref, ConsensusOptions(mark_snv="lower", mark_del="-"))
    text = marked.fetch(("chr1", 0, 40)).text
    assert text[4] == "g" and "-" in text
    mask = tmp_path / "mask.bed"
    mask.write_text("chr1\t0\t6\n")
    masked = ConsensusReader(vcf, ref, ConsensusOptions(mask=mask))
    assert masked.fetch(("chr1", 0, 8)).text == "NNNNNN" + REF[6:8]
    with pytest.raises(ValueError):
        ConsensusOptions(include="DP>1", exclude="DP>2")


def test_errors_and_pickle(files):
    vcf, ref = files
    c = ConsensusReader(vcf, SequenceReader(ref))
    with pytest.raises(errors.InvalidError):
        c.fetch(("chr1", 0, 10), Haplotype("nobody"))
    with pytest.raises(errors.InvalidError):
        c.fetch(("chr1", 0, 10), Haplotype("s1", "9x"))
    with pytest.raises(errors.RegionOutOfRangeError):
        c.fetch(("chr1", 0, 41))
    with pytest.raises(errors.OmicsWeaveError):
        ConsensusReader(vcf, ref, ConsensusOptions(include="NOSUCH>1"))
    d = pickle.loads(pickle.dumps(c))
    assert d.reference is not c.reference and not d.reference.closed
    assert d.fetch(("chr1", 0, 40)).text == c.fetch(("chr1", 0, 40)).text
    assert isinstance(d.fetch(("chr1", 0, 5)).bases, np.ndarray)


def test_file_objects_and_ownership(files):
    vcf, ref = files
    shared = SequenceReader(ref)
    with ConsensusReader(GzipVcf(vcf), shared) as c:
        assert c.source == GzipVcf(vcf) and c.reference is shared
    assert c.closed and not shared.closed
    with ConsensusReader(GzipVcf(vcf), Fasta(ref)) as owned:
        assert owned.reference.source == Fasta(ref)
    assert owned.reference.closed
    with pytest.raises(TypeError, match="ConsensusReader reads"):
        ConsensusReader(Bed(vcf), shared)


@pytest.mark.parametrize("owns_reference", [False, True])
def test_failed_constructor_closes_only_its_own_reference(files, monkeypatch, owns_reference):
    vcf, ref = files
    references = []

    def fail(*args):
        references.append(args[1])
        raise RuntimeError("builder failed")

    monkeypatch.setattr(_native, "ConsensusBuilder", fail)
    shared = None if owns_reference else SequenceReader(ref)
    try:
        with pytest.raises(RuntimeError, match="builder failed"):
            ConsensusReader(vcf, ref if owns_reference else shared)
        if owns_reference:
            with pytest.raises(ValueError, match="closed"):
                references[0].sequence("chr1", 0, 4, np.empty(4, dtype=np.uint8))
        else:
            assert not shared.closed
            assert shared.fetch(("chr1", 0, 4)).text == REF[:4]
    finally:
        if shared is not None:
            shared.close()
