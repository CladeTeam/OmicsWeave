"""A running plan remains inspectable and owns its outputs exclusively."""

import os
import subprocess
import sys

import pytest

PLAN = r"""
import os
from pathlib import Path
import sys

import omicsweave as ow

directory, kind, action = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
if kind == "convert":
    source = directory / "in.fa"
    source.write_bytes(b">c\nACGTACGT\n")
    plan = ow.convert.plan(source, directory / "out.fa.gz", threads=1, memory="64M")
else:
    source = directory / "in.sam"
    source.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100\n"
                      "r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\t*\n")
    bam = directory / "in.bam"
    ow.convert.run(source, bam, threads=1, memory="64M")
    source = bam
    task = ow.quant.Task([bam], directory / "quant", [ow.quant.Coverage()])
    plan = ow.quant.plan(task, threads=1, memory="64M")

if action == "callbacks":
    seen = []
    def progress(*args):
        assert plan.budget.threads == 1
        assert 0 < plan.memory.total <= plan.budget.memory
        if kind == "quant":
            assert plan.inputs[0].path == source
            assert plan.stages
            assert plan.needs
        else:
            assert plan._native.formats == ("fasta", "bgzip-fasta")
            assert plan._native.paths[0] == str(source)
        try:
            plan.run()
        except RuntimeError as error:
            assert "already running" in str(error)
        else:
            raise AssertionError("reentrant run acquired the same output files")
        seen.append(args)
    plan.run(progress=progress)
    assert seen, "no callback completed"
    plan.run()
elif action == "error":
    saved = source.read_bytes()
    source.unlink()
    try:
        plan.run()
    except ow.errors.OmicsWeaveError:
        pass
    else:
        raise AssertionError("missing source was accepted")
    source.write_bytes(saved)
    plan.run()
else:
    child = os.fork()
    if child == 0:
        try:
            plan.run()
        except RuntimeError as error:
            os._exit(0 if "after fork" in str(error) else 2)
        except BaseException:
            os._exit(3)
        os._exit(1)
    _, status = os.waitpid(child, 0)
    assert os.waitstatus_to_exitcode(status) == 0, status
    plan.run()
"""


@pytest.mark.parametrize("kind", ["convert", "quant"])
@pytest.mark.parametrize("action", ["callbacks", "error", "fork"])
def test_plan_run_lifecycle(tmp_path, kind, action):
    if action == "fork" and not hasattr(os, "fork"):
        pytest.skip("POSIX fork required")
    process = subprocess.run(
        [sys.executable, "-c", PLAN, str(tmp_path), kind, action],
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert process.returncode == 0, (process.stdout.decode(), process.stderr.decode())
