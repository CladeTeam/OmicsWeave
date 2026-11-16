"""Inherited buffered sinks must not flush into files owned by the parent."""

import os
import signal
import subprocess
import sys
from pathlib import Path

import pytest

WRITER_FORK = r"""
import gc
import os
from pathlib import Path
import sys
import traceback

directory, kind, action = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
threads = int(sys.argv[4])
w = None
early_writer = None
early_cache = None

def collect_before_native_hook():
    global w, early_writer, early_cache
    w = None
    gc.collect()
    early_writer = ow.formats.FastaWriter(directory / "early.fa")
    early_cache = ow.cache.BlockCache(max_bytes=4096)

if action == "callback_gc":
    # A user callback can run before the extension's process marker hook.
    os.register_at_fork(after_in_child=collect_before_native_hook)

import numpy as np
import omicsweave as ow

bases = b"ACGT" * 250000
chroms = {"c": 100}
if kind in {"fa_plain", "fa_bgzf"}:
    path = directory / ("parent.fa.gz" if kind == "fa_bgzf" else "parent.fa")
    w = ow.formats.FastaWriter(path)
    w.add("c", bases)
elif kind == "bigwig":
    path = directory / "parent.bw"
    w = ow.formats.BigWigWriter(path, chroms, threads=threads)
    w.add("c", [0, 4], [4, 8], [1.5, 3.0])
elif kind == "bigbed":
    path = directory / "parent.bb"
    w = ow.formats.BigBedWriter(path, chroms, ow.formats.BigBedSchema(4), threads=threads)
    w.add("c", [0, 4], [4, 8], ["x", "y"])
else:
    path = directory / "parent.hic"
    header = ow.formats.HicHeader("test", chroms, (10,))
    if kind == "hic_builder":
        w = ow.formats.HicBuilder(path, header, threads=threads)
        w.add("c", 10, "c", 20, 7.0)
    else:
        w = ow.formats.HicWriter(path, header, threads=threads)
        w.write_matrix("c", "c", 10, [1], [2], [7.0])

before = path.read_bytes() if path.exists() else None
pid = os.fork()
if pid == 0:
    try:
        if action == "close":
            w.close()
        w = None
        gc.collect()
        if early_writer is not None:
            early_writer.add("early", b"TGCA")
            early_writer.close()
            assert early_cache.stats()["resident_entries"] == 0
        with ow.formats.FastaWriter(directory / "child.fa") as fresh:
            fresh.add("child", b"ACGT")
        os._exit(0)
    except BaseException:
        traceback.print_exc()
        os._exit(1)

_, status = os.waitpid(pid, 0)
assert os.waitstatus_to_exitcode(status) == 0
after = path.read_bytes() if path.exists() else None
assert after == before, ("child changed the parent's file", kind, action, before, after)
w.close()
with ow.query.SequenceReader(directory / "child.fa") as child_reference:
    assert child_reference.fetch(("child", 0, 4)).bases.tobytes() == b"ACGT"
if action == "callback_gc":
    with ow.query.SequenceReader(directory / "early.fa") as early_reference:
        assert early_reference.fetch(("early", 0, 4)).bases.tobytes() == b"TGCA"
if kind in {"fa_plain", "fa_bgzf"}:
    with ow.query.SequenceReader(path) as reader:
        assert reader.fetch(("c", 0, len(bases))).bases.tobytes() == bases
elif kind == "bigwig":
    with ow.query.TrackReader(path) as reader:
        np.testing.assert_array_equal(reader.signal(("c", 0, 8)).values, [1.5]*4 + [3.0]*4)
elif kind == "bigbed":
    with ow.query.IntervalReader(path) as reader:
        assert reader.fetch(("c", 0, 8)).annotation.to_list() == ["x", "y"]
else:
    q = ow.types.ContactWindow(("c", 0, 40), ("c", 0, 40), 10)
    expected = np.zeros((4, 4), dtype=np.float32)
    expected[1, 2] = expected[2, 1] = 7.0
    with ow.query.ContactReader(path) as reader:
        np.testing.assert_array_equal(reader.matrix(q).values, expected)
"""


def exercise_fork_writer(directory: Path, kind: str, action: str, threads: int = 1) -> None:
    process = subprocess.Popen(
        [sys.executable, "-c", WRITER_FORK, str(directory), kind, action, str(threads)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=45)
    except subprocess.TimeoutExpired:
        # The session contains only this test and its own fork child.
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        pytest.fail(f"{kind}/{action}: forked writer did not finish within 45s")
    assert process.returncode == 0, (
        stdout.decode(errors="replace"),
        stderr.decode(errors="replace"),
    )


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork required")
@pytest.mark.parametrize("kind", ["fa_plain", "fa_bgzf", "bigwig", "bigbed", "hic"])
@pytest.mark.parametrize("action", ["gc", "close"])
def test_forked_writer_close_and_gc_leave_parent_file_unchanged(tmp_path, kind, action):
    exercise_fork_writer(tmp_path, kind, action)


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork required")
def test_buffered_writer_gc_before_native_fork_hook_does_not_flush(tmp_path):
    exercise_fork_writer(tmp_path, "fa_bgzf", "callback_gc")


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork required")
@pytest.mark.parametrize("kind", ["bigwig", "bigbed", "hic", "hic_builder"])
@pytest.mark.parametrize("action", ["gc", "close", "callback_gc"])
def test_forked_writer_pool_does_not_join_inherited_threads(tmp_path, kind, action):
    exercise_fork_writer(tmp_path, kind, action, threads=2)
