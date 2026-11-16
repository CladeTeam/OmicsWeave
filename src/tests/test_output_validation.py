"""Dense output allocation follows key and source validation."""

import json
import subprocess
import sys

import numpy as np
import pytest

import omicsweave as ow

CHECK_BEFORE_ALLOCATION = r"""
import json
import os
from pathlib import Path
import resource
import sys
import omicsweave as ow

path=Path(sys.argv[1])/"reference.fa"
with ow.formats.FastaWriter(path) as writer:
    writer.add("c",b"ACGTACGT")
reader=ow.query.SequenceReader(path)
reader.fetch(("c",0,4))
large=Path(sys.argv[1])/"large.fa"
with large.open("wb") as file:
    file.write(b">long\n")
    file.truncate((1<<30)+7)
large.with_name("large.fa.fai").write_text("long\t1073741824\t6\t1073741824\t1073741825\n")
large_reader=ow.query.SequenceReader(large)
pages=int(Path("/proc/self/statm").read_text().split()[0])
limit=pages*os.sysconf("SC_PAGE_SIZE")+(64<<20)
resource.setrlimit(resource.RLIMIT_AS,(limit,limit))
results={}
for name,key in [
    ("overlap",("c",[(0,2**32-1),(0,2**32-1)],".")),
    ("reversed",("c",[(4,2)],".")),
    ("unknown",("missing",[(0,2**32-1)],".")),
    ("outside",("c",[(0,2**32-1)],".")),
]:
    try:
        reader._native.fetch(key)
        results[name]="allowed"
    except BaseException as error:
        results[name]=type(error).__name__
reader.close()
try:
    reader._native.fetch(("c",[(0,4)],"."))
except BaseException as error:
    results["closed"]=type(error).__name__
try:
    large_reader._native.fetch(("long",[(0,1<<30)],"."))
except MemoryError:
    results["valid_oom"]="MemoryError"
except BaseException as error:
    results["valid_oom"]=type(error).__name__
print(json.dumps(results),flush=True)
"""


@pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="Linux address-space resource guard"
)
def test_native_sequence_rejects_invalid_keys_before_large_allocation(tmp_path):
    result = subprocess.run(
        [sys.executable, "-c", CHECK_BEFORE_ALLOCATION, str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert json.loads(result.stdout) == {
        "overlap": "InvalidError",
        "reversed": "InvalidRegionError",
        "unknown": "UnknownChromError",
        "outside": "RegionOutOfRangeError",
        "closed": "ValueError",
        "valid_oom": "MemoryError",
    }


def test_unknown_dense_keys_fail_before_default_storage(tmp_path, monkeypatch):
    fasta, track, contacts = tmp_path / "ref.fa", tmp_path / "track.bw", tmp_path / "contacts.hic"
    with ow.formats.FastaWriter(fasta) as writer:
        writer.add("c", b"ACGTACGT")
    with ow.formats.BigWigWriter(track, {"c": 16}) as writer:
        writer.add("c", [0], [16], [1.0])
    with ow.formats.HicWriter(contacts, ow.formats.HicHeader("test", {"c": 16}, (1,))) as writer:
        writer.write_matrix("c", "c", 1, [0], [0], [2.0])
    with (
        ow.query.SequenceReader(fasta) as sequence,
        ow.query.TrackReader(track) as signal,
        ow.query.ContactReader(contacts) as contact,
    ):

        def allocate(*args, **kwargs):
            pytest.fail("invalid key allocated default output")

        monkeypatch.setattr(np, "empty", allocate)
        with pytest.raises(ow.errors.RegionOutOfRangeError):
            sequence.fetch(("c", 0, 2**32 - 1))
        with pytest.raises(ow.errors.UnknownChromError):
            signal.signal(("missing", 0, 2**32 - 1))
        with pytest.raises(ow.errors.UnknownChromError):
            signal.bins(("missing", 0, 2**32 - 1), 1)
        window = ow.types.ContactWindow(("missing", 0, 2**32 - 1), ("c", 0, 16), 1)
        with pytest.raises(ow.errors.UnknownChromError):
            contact.matrix(window)


@pytest.mark.parametrize("length", [0, 8, 4096, 8192])
def test_sequence_default_storage_owns_resizable_memory(tmp_path, length):
    path = tmp_path / "reference.fa"
    bases = b"ACGT" * 2048
    with ow.formats.FastaWriter(path) as writer:
        writer.add("c", bases)
    with ow.query.SequenceReader(path) as reader:
        output = reader.fetch(("c", 0, length)).bases
        assert output.flags.owndata and output.base is None
        output.resize((length + 1,))
        assert output[:length].tobytes() == bases[:length]
        assert output[-1] == 0
