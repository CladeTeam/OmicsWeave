"""Process checks run before user callbacks can touch inherited state."""

import os
import subprocess
import sys

import pytest

EARLY_CALLBACK = r"""
import json
import os
from pathlib import Path
import sys

def before_library_callback():
    outcomes = {}
    for name, operation in [("reader", lambda: reader.fetch(("c", 0, 4))),
                            ("cache", cache.stats)]:
        try:
            operation()
            outcomes[name] = "allowed"
        except BaseException as error:
            outcomes[name] = type(error).__name__
    fresh = ow.cache.BlockCache(max_bytes=1024)
    assert fresh.stats()["resident_entries"] == 0
    os.write(child_pipe, json.dumps(outcomes).encode())

# Python runs child callbacks in registration order. This callback therefore
# runs before the one installed while importing the extension.
os.register_at_fork(after_in_child=before_library_callback)
import omicsweave as ow

path = Path(sys.argv[1]) / "reference.fa"
with ow.formats.FastaWriter(path) as writer:
    writer.add("c", b"ACGTACGT")
reader = ow.query.SequenceReader(path)
cache = ow.cache.BlockCache(max_bytes=1024)
assert reader.fetch(("c", 0, 4)).bases.tobytes() == b"ACGT"
parent_pipe, child_pipe = os.pipe()
pid = os.fork()
if pid == 0:
    os._exit(0)
os.close(child_pipe)
outcomes = json.loads(os.read(parent_pipe, 8192))
os.close(parent_pipe)
_, status = os.waitpid(pid, 0)
assert os.waitstatus_to_exitcode(status) == 0
assert outcomes == {"reader": "RuntimeError", "cache": "RuntimeError"}, outcomes
assert reader.fetch(("c", 0, 4)).bases.tobytes() == b"ACGT"
assert cache.stats()["resident_entries"] == 0
reader.close()
"""


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork required")
def test_early_child_callback_cannot_query_inherited_reader_or_cache(tmp_path):
    process = subprocess.run(
        [sys.executable, "-c", EARLY_CALLBACK, str(tmp_path)],
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert process.returncode == 0, (process.stdout.decode(), process.stderr.decode())
