import pickle

import numpy as np
import pytest

import omicsweave as ow
from omicsweave import errors
from omicsweave.formats import Bcf, Csi
from omicsweave.query import VariantReader
from omicsweave.types import GT_ABSENT, GT_MISSING, Genotypes

from ..support.variants import VCF


@pytest.fixture
def bcf(tmp_path):
    (tmp_path / "v.vcf").write_text(VCF)
    report = ow.convert.run(tmp_path / "v.vcf", tmp_path / "v.bcf", index=True, threads=1)
    assert report.companions == (Csi(tmp_path / "v.bcf.csi"),)
    return report.target


def test_genotypes_of_a_region(bcf):
    with VariantReader(bcf) as reader:
        assert reader.samples == ("a", "b", "c")
        assert dict(reader.chroms) == {"chr1": 1000, "chr2": 500}
        g = reader.genotypes(("chr1", 0, 1000))
    assert isinstance(g, Genotypes) and len(g) == 4 and g.ploidy == 2
    assert g.start.tolist() == [4, 9, 19, 299]
    assert g.end.tolist() == [5, 10, 120, 300]
    assert g.ids.to_list() == ["rs1", ".", ".", "."]
    assert g.ref.to_list() == ["A", "C", "G", "T"]
    assert g.alt.to_list() == ["G", "CTT,G", "<DEL>", ""]
    assert g.gt.dtype == np.int8 and g.gt.shape == (4, 3, 2)
    m, x = GT_MISSING, GT_ABSENT
    np.testing.assert_array_equal(
        g.gt,
        [
            [[0, 1], [1, 1], [m, m]],
            [[2, 0], [0, 0], [1, x]],
            [[0, 1], [m, x], [1, 1]],
            [[x, x], [x, x], [x, x]],
        ],
    )
    np.testing.assert_array_equal(
        g.phased,
        # A lone allele takes the phase htslib gives it: before VCF 4.4, as phased
        # as the padding after it in a diploid record.
        [[True, False, False], [True, True, True], [False, True, True], [False] * 3],
    )


def test_regions_and_sample_order(bcf):
    reader = VariantReader(bcf)
    g = reader.genotypes(("chr1", 50, 60), ["c", "a"])
    # Only the deletion spans [50, 60).
    assert g.start.tolist() == [19] and g.samples == ("c", "a")
    np.testing.assert_array_equal(g.gt, [[[1, 1], [0, 1]]])
    haploid = reader.genotypes(("chr2", 0, 500), "a")
    assert haploid.ploidy == 1 and haploid.gt.tolist() == [[[1]]]
    assert haploid.phased.tolist() == [[True]]
    empty = reader.genotypes(("chr1", 500, 600))
    assert len(empty) == 0 and empty.gt.shape == (0, 3, 0)
    with pytest.raises(KeyError, match="no sample 'z'"):
        reader.genotypes(("chr1", 0, 10), ["z"])
    with pytest.raises(errors.UnknownChromError):
        reader.genotypes(("chrX", 0, 10))
    again = pickle.loads(pickle.dumps(reader))
    assert again.genotypes(("chr1", 0, 6)).start.tolist() == [4]
    assert isinstance(ow.open(bcf), VariantReader)


def test_index_is_required(tmp_path):
    (tmp_path / "v.vcf").write_text(VCF)
    report = ow.convert.run(tmp_path / "v.vcf", tmp_path / "v.bcf", threads=1)
    assert report.target == Bcf(tmp_path / "v.bcf") and report.companions == ()
    with pytest.raises(errors.IoError, match="no index"):
        VariantReader(report.target)
    with pytest.raises(errors.InvalidError):
        VariantReader(tmp_path / "v.vcf")
