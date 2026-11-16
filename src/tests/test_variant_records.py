import numpy as np
import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.query import VariantReader
from omicsweave.types import FLOAT_END, INT_END, INT_MISSING, Variants

VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FILTER=<ID=q10,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="Known">
##INFO=<ID=NOTE,Number=1,Type=String,Description="Note">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb
chr1\t5\trs1\tA\tG\t50\tPASS\tDP=10;AF=0.5;DB;NOTE=x\tGT:AD:FT\t0/1:3,7:ok\t1/1:0,9:ok
chr1\t10\t.\tC\tCT,G\t.\tq10\tDP=3;AF=0.1,0.2\tGT:AD\t0/2:1,0,2\t./.:.
chr1\t20\t.\tG\tT\t20\t.\t.\tGT\t0/0\t0/1
"""


@pytest.fixture
def bcf(tmp_path):
    (tmp_path / "v.vcf").write_text(VCF)
    return ow.convert.run(tmp_path / "v.vcf", tmp_path / "v.bcf", index=True, threads=1).target


def bits(a):
    return np.asarray(a, dtype=np.float32).view(np.uint32)


def test_records_carry_qual_filter_info_and_format(bcf):
    with VariantReader(bcf) as reader:
        v = reader.fetch(("chr1", 0, 1000), info=["DP", "AF", "DB", "NOTE"], format=["AD", "FT"])
        g = reader.genotypes(("chr1", 0, 1000))
    assert isinstance(v, Variants) and len(v) == 3
    np.testing.assert_array_equal(v.genotypes.gt, g.gt)
    np.testing.assert_array_equal(v.genotypes.phased, g.phased)
    assert v.genotypes.start.tolist() == [4, 9, 19]
    assert v.qual[0] == 50 and np.isnan(v.qual[1]) and v.qual[2] == 20
    assert v.filters.to_list() == ["PASS", "q10", ""]

    assert v.info["DP"].tolist() == [[10], [3], [INT_END]]
    af = v.info["AF"]
    assert af.shape == (3, 2)
    np.testing.assert_allclose(af[1], [0.1, 0.2], rtol=1e-6)
    assert af[0, 0] == 0.5 and bits(af[0, 1]) == FLOAT_END
    assert (bits(af[2]) == FLOAT_END).all()
    assert v.info["DB"].tolist() == [True, False, False]
    assert v.info["NOTE"].to_list() == ["x", "", ""]

    e, m = INT_END, INT_MISSING
    assert v.format["AD"].shape == (3, 2, 3)
    assert v.format["AD"].tolist() == [
        [[3, 7, e], [0, 9, e]],
        [[1, 0, 2], [m, e, e]],
        [[e, e, e], [e, e, e]],
    ]
    assert v.format["FT"].to_list() == ["ok", "ok", "", "", "", ""]


def test_samples_and_filter_expressions(bcf):
    with VariantReader(bcf) as reader:
        v = reader.fetch(("chr1", 0, 1000), ["b"], format=["AD"], include="QUAL>30")
        assert v.genotypes.samples == ("b",)
        assert v.genotypes.start.tolist() == [4]
        assert v.format["AD"].tolist() == [[[0, 9]]]
        kept = reader.fetch(("chr1", 0, 1000), exclude="FILTER='q10'")
        assert kept.genotypes.start.tolist() == [4, 19]
        assert reader.fetch(("chr1", 0, 8)).filters.to_list() == ["PASS"]


def test_field_errors(bcf):
    with VariantReader(bcf) as reader:
        region = ("chr1", 0, 1000)
        with pytest.raises(errors.InvalidError, match="not declared"):
            reader.fetch(region, info=["XX"])
        with pytest.raises(errors.InvalidError, match="genotype matrix"):
            reader.fetch(region, format=["GT"])
        with pytest.raises(errors.InvalidError, match="include or exclude"):
            reader.fetch(region, include="QUAL>1", exclude="QUAL>2")
        with pytest.raises(ValueError, match="twice"):
            reader.fetch(region, info=["DP", "DP"])
        with pytest.raises(KeyError, match="nobody"):
            reader.fetch(region, ["nobody"])


def test_records_and_genotypes_agree_on_mixed_ploidy(tmp_path):
    from test_variant import VCF as MIXED

    (tmp_path / "m.vcf").write_text(MIXED)
    bcf = ow.convert.run(tmp_path / "m.vcf", tmp_path / "m.bcf", index=True, threads=1).target
    with VariantReader(bcf) as reader:
        for region in [("chr1", 0, 1000), ("chr2", 0, 500), ("chr1", 15, 25)]:
            g = reader.genotypes(region, ["c", "a"])
            v = reader.fetch(region, ["c", "a"]).genotypes
            for name in ("start", "end", "gt", "phased"):
                np.testing.assert_array_equal(getattr(v, name), getattr(g, name))
            assert v.alt.to_list() == g.alt.to_list()
