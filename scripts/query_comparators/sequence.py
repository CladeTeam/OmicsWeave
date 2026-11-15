"""Real sequence APIs on identical on-disk FASTA/BGZF/2bit fixtures.

Fixtures: ref.fa, ref.fa.fai; ref.fa.gz, ref.fa.gz.fai, ref.fa.gz.gzi;
ref.2bit. Every file must contain identical chr1/chr2 sequences of 2**22
uppercase A/C/G/T/N bases. All indexes must already exist. No fixture or
index is built here. CLI measurements include subprocess and FASTA formatting.

A reused Omics output is borrowed until the next query on that backend;
normalize each raw result before issuing the next query when validating.
"""

from __future__ import annotations

import os
import random
import subprocess
from pathlib import Path
from typing import Any

Query = tuple[str, int, int]
_SEQUENCE_LENGTH = 1 << 22


def _checked(bases: bytes, query: Query) -> bytes:
    expected = query[2] - query[1]
    if len(bases) != expected:
        raise ValueError(f"{query}: expected {expected} bases, received {len(bases)}")
    return bases


def _native_bytes(raw: Any, query: Query) -> bytes:
    return _checked(raw.bases.tobytes(), query)


def _ascii_bytes(raw: str, query: Query) -> bytes:
    return _checked(raw.encode("ascii"), query)


def _fasta_bytes(raw: bytes, query: Query) -> bytes:
    # splitlines handles LF/CRLF, without removing any sequence characters.
    bases = b"".join(line for line in raw.splitlines() if not line.startswith(b">"))
    return _checked(bases, query)


def _queries(length: int) -> list[Query]:
    rng = random.Random(0xFA57A008 + length)
    queries = []
    for i in range(16):
        chrom = f"chr{1 + i % 2}"
        start = rng.randrange(_SEQUENCE_LENGTH - length + 1)
        queries.append((chrom, start, start + length))
    return queries


def build_cases(ow: Any, fixtures: Path, tools: dict) -> list[dict]:
    """Return factory-driven cases; query calls fully materialize native results."""
    import numpy as np
    import py2bit
    import pyfaidx
    import pysam

    fixtures = Path(fixtures)
    paths = {
        "fasta": fixtures / "ref.fa",
        "bgzf": fixtures / "ref.fa.gz",
        "2bit": fixtures / "ref.2bit",
    }
    required = [
        paths["fasta"],
        Path(f"{paths['fasta']}.fai"),
        paths["bgzf"],
        Path(f"{paths['bgzf']}.fai"),
        Path(f"{paths['bgzf']}.gzi"),
        paths["2bit"],
    ]
    missing = [os.fspath(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"sequence fixtures/indexes missing: {missing}")

    samtools = os.fspath(tools["samtools"])
    twobit_to_fa = os.fspath(tools["twoBitToFa"])

    def omics_factory(path: Path, length: int, reuse_out: bool):
        def factory():
            reader = ow.query.SequenceReader(path, soft_mask=False)
            out = np.empty(length, dtype=np.uint8) if reuse_out else None

            def query(q: Query):
                return reader.fetch(q, out=out)

            return {
                "query": query,
                "normalize": _native_bytes,
                "close": reader.close,
                "reader": reader,
            }

        return factory

    def pysam_factory(path: Path):
        def factory():
            reader = pysam.FastaFile(os.fspath(path))

            def query(q: Query):
                return reader.fetch(q[0], q[1], q[2])

            return {"query": query, "normalize": _ascii_bytes, "close": reader.close}

        return factory

    def pyfaidx_factory(path: Path):
        def factory():
            reader = pyfaidx.Fasta(
                os.fspath(path),
                as_raw=True,
                sequence_always_upper=False,
                build_index=False,
                rebuild=False,
            )

            def query(q: Query):
                # as_raw returns the complete sequence string at the API boundary.
                return reader[q[0]][q[1] : q[2]]

            return {"query": query, "normalize": _ascii_bytes, "close": reader.close}

        return factory

    def py2bit_factory(path: Path):
        def factory():
            reader = py2bit.open(os.fspath(path))

            def query(q: Query):
                return reader.sequence(q[0], q[1], q[2])

            return {"query": query, "normalize": _ascii_bytes, "close": reader.close}

        return factory

    def cli_factory(path: Path, is_twobit: bool):
        def factory():
            def query(q: Query):
                chrom, start, end = q
                if is_twobit:
                    # UCSC flags use zero-based start and exclusive end.
                    argv = [
                        twobit_to_fa,
                        os.fspath(path),
                        "stdout",
                        f"-seq={chrom}",
                        f"-start={start}",
                        f"-end={end}",
                        "-noMask",
                    ]
                else:
                    # samtools region strings are one-based, end-inclusive.
                    argv = [
                        samtools,
                        "faidx",
                        "-n",
                        "0",
                        os.fspath(path),
                        f"{chrom}:{start + 1}-{end}",
                    ]
                return subprocess.run(argv, check=True, capture_output=True, timeout=30).stdout

            return {"query": query, "normalize": _fasta_bytes, "close": lambda: None}

        return factory

    cases = []
    for format_name, path in paths.items():
        for length in (128, 4096, 131072):
            backends = [
                {
                    "name": "omics.allocate",
                    "kind": "omics",
                    "factory": omics_factory(path, length, False),
                },
                {
                    "name": "omics.reuse_out",
                    "kind": "omics",
                    "factory": omics_factory(path, length, True),
                },
            ]
            if format_name == "2bit":
                backends.extend(
                    [
                        {
                            "name": "py2bit.sequence",
                            "kind": "python",
                            "factory": py2bit_factory(path),
                        },
                        {"name": "twoBitToFa", "kind": "cli", "factory": cli_factory(path, True)},
                    ]
                )
            else:
                backends.extend(
                    [
                        {"name": "pysam.fetch", "kind": "python", "factory": pysam_factory(path)},
                        {
                            "name": "pyfaidx.slice",
                            "kind": "python",
                            "factory": pyfaidx_factory(path),
                        },
                        {
                            "name": "samtools.faidx",
                            "kind": "cli",
                            "factory": cli_factory(path, False),
                        },
                    ]
                )
            cases.append(
                {
                    "id": f"sequence.{format_name}.random{length}",
                    "interface": "SequenceReader.fetch",
                    "format": format_name,
                    "workload": f"16 fixed-seed alternating-chromosome queries, {length} bases each",
                    "queries": _queries(length),
                    "backends": backends,
                }
            )
    return cases
