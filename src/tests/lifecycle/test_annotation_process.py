"""Lazy native annotation indexes belong to the process that created them."""

import os
import subprocess
import sys

import pytest

ANNOTATION_FORK = r"""
import os
from pathlib import Path
import pickle
import sys
import omicsweave as ow
from omicsweave.annotation import Annotation

directory = Path(sys.argv[1])
path = directory / "genes.gtf"
path.write_text('c\tx\texon\t1\t4\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
ann = Annotation.read(path, threads=1)
intervals = ann.exonic_regions()
region = ("c", 0, 4)
assert ann.overlaps(region).tolist() == [0]

pid = os.fork()
if pid == 0:
    try:
        operations = [
            lambda: ann.overlaps(region),
            lambda: ann.gene_rows("g"),
            lambda: ann.transcripts,
            lambda: intervals.overlaps(region),
            lambda: ow.quant.plan(ow.quant.Task([directory / "missing.bam"],
                                               directory / "output", annotation=ann)),
        ]
        for operation in operations:
            try:
                operation()
            except RuntimeError as error:
                assert "fork" in str(error)
            else:
                raise AssertionError("inherited lazy annotation state was queried")
        # Exporting ordinary columns does not enter lazy index initialization.
        assert intervals.start.tolist() == [0]
        restored = pickle.loads(pickle.dumps(ann))
        assert restored.overlaps(region).tolist() == [0]
        assert restored.exonic_regions().overlaps(region).tolist() == [0]
        os._exit(0)
    except BaseException:
        import traceback
        traceback.print_exc()
        os._exit(1)
_, status = os.waitpid(pid, 0)
assert os.waitstatus_to_exitcode(status) == 0
assert ann.overlaps(region).tolist() == [0]
assert intervals.overlaps(region).tolist() == [0]
"""


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork required")
def test_inherited_annotation_queries_reopen_before_lazy_initialization(tmp_path):
    process = subprocess.run(
        [sys.executable, "-c", ANNOTATION_FORK, str(tmp_path)],
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert process.returncode == 0, (process.stdout.decode(), process.stderr.decode())
