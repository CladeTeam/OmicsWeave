import pickle

import pytest

from omicsweave import quant
from omicsweave.cli.app import main
from omicsweave.errors import OmicsWeaveError

from ._support import CHROMS, write_bam
from ._support import data as data


def test_refusals(tmp_path, data):
    _, recs, _, _ = data
    swapped = list(recs)
    swapped[10], swapped[500] = swapped[500], swapped[10]
    bam = write_bam(tmp_path, swapped, "swapped")
    task = quant.Task([bam], tmp_path / "o", [quant.Coverage()])
    with pytest.raises(OmicsWeaveError, match="not sorted by coordinate"):
        quant.run(task, threads=1)
    unsorted = write_bam(tmp_path, recs, "unsorted", "@HD\tVN:1.6\tSO:unsorted\n")
    with pytest.raises(OmicsWeaveError, match="--sort"):
        quant.plan(quant.Task([unsorted], tmp_path / "o", [quant.Coverage()]))
    with pytest.raises(OmicsWeaveError, match="needs at least"):
        quant.plan(task, threads=4, memory="4M")
    with pytest.raises(OmicsWeaveError, match="annotation"):
        quant.plan(quant.Task([bam], tmp_path / "o", [quant.Counts.rna_single_cell()]))
    assert main(["quant", str(bam), "-o", str(tmp_path / "x")]) == 2


def test_presets_come_from_the_native_layer():
    single, bulk = quant.Counts.rna_single_cell(), quant.Counts.rna_bulk(unit="pair")
    assert (single.rows, single.dedup, single.strand) == ("CB", "umi", "forward")
    assert single.filter.unique and isinstance(single.features, quant.GeneExons)
    assert (bulk.rows, bulk.unit, bulk.table, bulk.strand) == ("@file", "pair", True, "unstranded")
    chip = quant.Counts.chip(quant.Peaks("p.bb"))
    assert chip.features == quant.Peaks("p.bb") and chip.multi_overlap
    assert quant.Atac().gene_activity == quant.GeneActivity(2000, 0)
    assert quant.Alleles.bulk("s.bcf").unit == "fragment"
    assert quant.JunctionMatrix.bulk().rows == "@file"
    assert quant.Contacts().resolutions[-1] == 5_000


def test_products_are_checked_by_the_native_layer(data):
    tmp, _, bam, _ = data

    def plan(*products):
        return quant.plan(quant.Task([bam], tmp / "x", products), memory="256M")

    with pytest.raises(ValueError, match="split strands"):
        plan(quant.Coverage(strand="forward"))
    with pytest.raises(ValueError, match="fragment unit"):
        plan(quant.Coverage(max_fragment=500))
    with pytest.raises(ValueError, match="strand must be"):
        plan(quant.Coverage(strands="split", strand="sense"))
    with pytest.raises(ValueError, match="two characters"):
        plan(quant.Counts(quant.Bins(1000), umi_tag="UMI"))
    with pytest.raises(ValueError, match="none or umi"):
        plan(quant.Junctions(matrix=quant.JunctionMatrix(dedup="position")))
    with pytest.raises(AttributeError):
        quant.plan(quant.Task([bam], tmp / "x", [object()]))


def test_a_plan_runs_without_planning_again(data):
    tmp, recs, bam, _ = data
    p = quant.plan(quant.Task([bam], tmp / "once", [quant.Coverage()]), threads=1, memory="256M")
    ((path, refs, groups, size),) = [
        (i.path, i.references, i.read_groups, i.size) for i in p.inputs
    ]
    assert path == bam and refs == tuple(CHROMS) and groups == () and size > 0
    report = p.run(progress=lambda *a: None)
    assert report.records == len(recs) and 0 <= report.scan_seconds <= report.seconds
    assert [f.kind for f in report.files if f.path.name == "coverage.bw"] == ["bigwig"]
    with pytest.raises(TypeError, match="process-local"):
        pickle.dumps(p)
