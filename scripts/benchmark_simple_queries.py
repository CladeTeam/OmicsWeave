"""One loaded OmicsWeave query vs loaded bindings and CLI -> file -> read."""

from __future__ import annotations

import argparse
import bisect
import ctypes
import gc
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
import platform
import resource
import shutil
import statistics
import struct
import subprocess
import sys
import time
from pathlib import Path


class FileSubprocess:
    """Adapter-local subprocess interface; the process-wide module is unchanged."""

    def __init__(self, directory: Path, name: str):
        self.path = directory / (name + ".out")
        self.calls = 0

    def __getattr__(self, name):
        return getattr(subprocess, name)

    def run(self, argv, **kwargs):
        argv = [str(self.path) if a in {"stdout", "/dev/stdout"} else a for a in argv]
        kwargs.pop("capture_output", None)
        kwargs.pop("stdout", None)
        kwargs.setdefault("stderr", subprocess.PIPE)
        with self.path.open("wb") as output:
            result = subprocess.run(argv, stdout=output, **kwargs)
        data = self.path.read_bytes()
        self.calls += 1
        return subprocess.CompletedProcess(result.args, result.returncode, data, result.stderr)


def atomic_json(path, value):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False))
    temp.replace(path)


def cpu_ns():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time_ns() + round((child.ru_utime + child.ru_stime) * 1e9)


def pyfaidx_without_block_reuse(path, normalize):
    import pyfaidx

    reader = pyfaidx.Fasta(
        os.fspath(path),
        as_raw=True,
        sequence_always_upper=False,
        build_index=False,
        rebuild=False,
        read_ahead=0,
    )
    if reader.faidx._bgzf:
        # Bio.bgzf requires at least one block. Alternating chromosomes make
        # consecutive queries disjoint; the diagnostic verifies actual inflate
        # calls for every query, including across repetitions of the workload.
        reader.faidx.file.max_cache = 1
        reader.faidx.file._buffers.clear()

    def query(q):
        return reader[q[0]][q[1] : q[2]]

    return {"query": query, "normalize": normalize, "close": reader.close}


def bgzf_query_blocks(fixtures, queries):
    fai = {}
    for line in (fixtures / "ref.fa.gz.fai").read_text().splitlines():
        chrom, _, offset, bases, width, *_ = line.split("\t")
        fai[chrom] = (int(offset), int(bases), int(width))
    data = (fixtures / "ref.fa.gz.gzi").read_bytes()
    count = struct.unpack_from("<Q", data)[0]
    assert len(data) == 8 + count * 16
    starts = [(0, 0), *struct.iter_unpack("<QQ", data[8:])]
    offsets = [point[1] for point in starts]
    spans = []
    for chrom, start, end in queries:
        offset, bases, width = fai[chrom]
        lo = offset + start // bases * width + start % bases
        hi = offset + (end - 1) // bases * width + (end - 1) % bases
        first = bisect.bisect_right(offsets, lo) - 1
        last = bisect.bisect_right(offsets, hi) - 1
        spans.append([starts[i][0] for i in range(first, last + 1)])
    for previous, current in zip(spans, [*spans[1:], spans[0]], strict=True):
        assert not set(previous).intersection(current), "adjacent BGZF queries share a block"
    return spans


def audit_queries(obj, queries, reference, probe, pyfaidx_bgzf, expected_blocks):
    counts, loaded_blocks = [], []
    original = None
    if pyfaidx_bgzf:
        from Bio import bgzf

        original = bgzf._load_bgzf_block

        def counted(handle, *args, **kwargs):
            loaded_blocks.append(handle.tell())
            return original(handle, *args, **kwargs)

        bgzf._load_bgzf_block = counted
    try:
        for _ in range(2):
            for i, (q, expected) in enumerate(zip(queries, reference, strict=True)):
                before = len(loaded_blocks) if original else probe.ow_probe_counts()
                raw = obj["query"](q)
                after = len(loaded_blocks) if original else probe.ow_probe_counts()
                counts.append(after - before)
                if original:
                    assert loaded_blocks[before:after] == expected_blocks[i], (
                        "BGZF query did not freshly decode every required block",
                        q,
                        loaded_blocks[before:after],
                        expected_blocks[i],
                    )
                assert obj["normalize"](raw, q) == expected
    finally:
        if original:
            bgzf._load_bgzf_block = original
    return {
        "method": "Bio.bgzf._load_bgzf_block" if original else "zlib inflate/uncompress interpose",
        "queries_after_warmup": len(counts),
        "decompression_calls_per_query": counts,
        "minimum": min(counts),
        "all_required_blocks_freshly_loaded": True if original else None,
    }


def ncls_direct(path, normalize):
    def factory():
        import numpy as np
        from ncls import NCLS

        groups = {}
        with path.open("rb") as source:
            for line in source:
                c, s, e, note = line.rstrip(b"\n").split(b"\t", 3)
                groups.setdefault(c.decode(), []).append((int(s), int(e), note))
        indexes = {
            c: NCLS(
                np.array([r[0] for r in rr], dtype=np.int64),
                np.array([r[1] for r in rr], dtype=np.int64),
                np.arange(len(rr), dtype=np.int64),
            )
            for c, rr in groups.items()
        }

        def query(q):
            return list(indexes[q[0]].find_overlap(q[1], q[2]))

        def check(raw, q):
            ids = sorted(int(row[2]) for row in raw)
            return normalize([groups[q[0]][i] for i in ids], q)

        def close():
            indexes.clear()
            groups.clear()

        return {"query": query, "normalize": check, "close": close}

    return factory


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--tools", type=Path, default=Path("/tmp/omicsweave-final-r8-tools"))
    p.add_argument("--cli-output-dir", type=Path, required=True)
    p.add_argument("--filter", default="")
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--sample-ms", type=int, default=100)
    p.add_argument("--cpu", type=int)
    p.add_argument("--validate-only", action="store_true")
    p.add_argument("--cache-audit-library", type=Path)
    p.add_argument("--ours-cache-bytes", type=int)
    p.add_argument("--ours-cache-entries", type=int, default=8192)
    args = p.parse_args()
    probe = None
    if args.cache_audit_library:
        if not args.validate_only:
            p.error("cache audit must run separately from timing with --validate-only")
        if args.ours_cache_bytes is not None:
            p.error("decoded-cache audit is for the uncached protocol only")
        probe = ctypes.CDLL(str(args.cache_audit_library.resolve()))
        probe.ow_probe_counts.restype = ctypes.c_uint64
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    import omicsweave as ow

    args.cli_output_dir.mkdir(parents=True, exist_ok=True)
    tools = {
        n: str(args.tools / n) for n in ["twoBitToFa", "bigWigToBedGraph", "bigBedToBed", "straw"]
    }
    tools.update({n: shutil.which(n) for n in ["samtools", "bcftools", "bedtools"]})
    cases, proxies = [], []
    for name in ["sequence", "track_interval", "consensus", "contact", "alignment_variant"]:
        spec = importlib.util.spec_from_file_location(
            "simple_" + name, Path(__file__).parent / "query_comparators" / (name + ".py")
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        proxy = FileSubprocess(args.cli_output_dir, name)
        module.subprocess = proxy
        proxies.append(proxy)
        cases += module.build_cases(ow, args.fixtures, tools)
    cases = [c for c in cases if args.filter in c["id"]]
    cache_enabled = args.ours_cache_bytes is not None
    if cache_enabled:
        if args.ours_cache_bytes < 0 or args.ours_cache_entries < 0:
            p.error("cache limits must be nonnegative")
        cases = [
            c
            for c in cases
            if c["format"] in {"bgzf", "2bit", "bigwig", "bigbed", "bam", "bcf"}
            or (c["format"].startswith("hic-") and not c["id"].endswith(".sites"))
        ]
    result = {
        "protocol": {
            "ours": (
                "one default public query; persistent loaded reader; "
                "preconstructed region/window; default allocation"
            ),
            "python_bindings": (
                "already imported; persistent reader/index where API supports one; "
                "direct native output"
            ),
            "cli": "each call subprocess -> ordinary output file -> wait -> read_bytes; no fsync",
            "cli_output_directory": str(args.cli_output_dir.resolve()),
            "excluded": [
                "imports",
                "reader/index/store construction",
                "query-argument construction",
                "extra output normalization",
            ],
            "validation": (
                "full logical payload normalization and byte equality before timing "
                "and after every sample, outside timed interval"
            ),
            "repeats": args.repeats,
            "sample_ms": args.sample_ms,
            "cache": (
                "explicit external OmicsWeave block cache; full sequence warmup; "
                "binding-native default caches retained; these cache policies differ"
                if cache_enabled
                else "warm OS cache; loaded metadata/indexes; no cross-query decoded-block hits; "
                "pyfaidx read_ahead=0 and Bio.bgzf max_cache=1 "
                "with disjoint consecutive query blocks"
            ),
            "ours_cache_limits": (
                {"bytes": args.ours_cache_bytes, "entries": args.ours_cache_entries}
                if cache_enabled
                else None
            ),
            "cache_audit_library": str(args.cache_audit_library) if probe else None,
            "affinity": sorted(os.sched_getaffinity(0)),
            "exclusive_cpu": False,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "benchmark_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "omics_native": ow._native.__file__,
            "omics_native_sha256": hashlib.sha256(
                Path(ow._native.__file__).read_bytes()
            ).hexdigest(),
            "versions": {
                n: metadata.version(n)
                for n in [
                    "numpy",
                    "pysam",
                    "pyfaidx",
                    "py2bit",
                    "pyBigWig",
                    "hic-straw",
                    "ncls",
                    "biopython",
                ]
            },
        },
        "cases": [],
        "rows": [],
        "errors": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    for case in cases:
        print("case", case["id"], flush=True)
        queries = case["queries"]
        expected_blocks = (
            bgzf_query_blocks(args.fixtures, queries) if case["format"] == "bgzf" else None
        )
        selected = []
        for b in case["backends"]:
            if b["kind"] != "omics" or not any(x["kind"] == "omics" for x in selected):
                selected.append(b)
        assert sum(b["kind"] == "omics" for b in selected) == 1
        contexts, reference, samples = {}, None, {}
        cache = (
            ow.cache.BlockCache(
                max_bytes=args.ours_cache_bytes, max_entries=args.ours_cache_entries
            )
            if cache_enabled
            else None
        )
        info = {k: v for k, v in case.items() if k not in {"queries", "backends"}}
        info.update(queries=queries, backends=[])
        if expected_blocks is not None:
            info["bgzf_required_compressed_block_offsets"] = expected_blocks
            info["cyclically_adjacent_queries_share_no_bgzf_blocks"] = True
        for b in selected:
            obj = None
            entry = {k: v for k, v in b.items() if k != "factory"}
            try:
                obj = b["factory"]()
                if b["kind"] == "omics" and cache is not None:
                    obj["reader"].set_cache(cache)
                if not cache_enabled and b["kind"] == "python" and b["name"].startswith("pyfaidx"):
                    normalize = obj["normalize"]
                    obj["close"]()
                    path = args.fixtures / (
                        "ref.fa.gz" if case["format"].lower() == "bgzf" else "ref.fa"
                    )
                    obj = pyfaidx_without_block_reuse(path, normalize)
                    entry["cache_policy"] = (
                        "read_ahead=0; BGZF max_cache=1; "
                        "consecutive queries touch disjoint chromosome blocks"
                    )
                if b["kind"] == "omics" and len(queries[0]) == 3:
                    regions = {q: ow.Region(*q) for q in queries}
                    fn = obj["query"]
                    obj["query"] = lambda q, fn=fn, regions=regions: fn(regions[q])
                if b["kind"] == "python" and b["name"].startswith("NCLS."):
                    norm = obj["normalize"]
                    obj["close"]()
                    obj = ncls_direct(args.fixtures / "features.bed", norm)()
                    entry["note"] = (
                        "direct find_overlap iterator fully consumed to list; "
                        "annotation lookup and stable order only in untimed validation"
                    )
                # pyBigWig native results are already materialized; no adapter list copy.
                if b["kind"] == "python" and b["name"].startswith("pyBigWig."):
                    import pyBigWig

                    norm = obj["normalize"]
                    obj["close"]()
                    is_bigbed = "entries" in b["name"]
                    reader = pyBigWig.open(
                        str(args.fixtures / ("features.bb" if is_bigbed else "signal.bw"))
                    )
                    if is_bigbed:

                        def fn(q, reader=reader):
                            return reader.entries(*q)
                    elif "intervals" in b["name"]:

                        def fn(q, reader=reader):
                            return reader.intervals(*q)
                    else:

                        def fn(q, reader=reader):
                            return reader.values(*q, numpy=True)

                    obj = {"query": fn, "normalize": norm, "close": reader.close}
                payloads = [obj["normalize"](obj["query"](q), q) for q in queries]
                if reference is None:
                    reference = payloads
                assert payloads == reference, ("logical payload mismatch", case["id"], b["name"])
                if probe and b["kind"] == "python":
                    entry["cache_audit"] = audit_queries(
                        obj,
                        queries,
                        reference,
                        probe,
                        b["name"].startswith("pyfaidx") and expected_blocks is not None,
                        expected_blocks,
                    )
                    format_name = case["format"].lower()
                    compressed = format_name in {
                        "bgzf",
                        "bigwig",
                        "bigbed",
                        "bed.gz",
                        "bam",
                        "bcf",
                    } or format_name.startswith("hic-")
                    if compressed:
                        assert entry["cache_audit"]["minimum"] > 0, (
                            "query bypassed decompression",
                            case["id"],
                            b["name"],
                        )
                entry["status"] = "ok"
                entry["checksums"] = [hashlib.sha256(v).hexdigest() for v in payloads]
                if b["kind"] == "omics" and cache is not None:
                    entry["cache_after_warmup"] = cache.stats()
                contexts[b["name"]] = obj
                samples[b["name"]] = []
            except Exception as e:
                entry.update(status="error", error=repr(e))
                result["errors"].append(
                    {"case": case["id"], "backend": b["name"], "error": repr(e)}
                )
                if obj:
                    obj["close"]()
                print("error", entry, flush=True)
            info["backends"].append(entry)
        result["cases"].append(info)
        try:
            if not args.validate_only:
                for repeat in range(args.repeats):
                    order = selected[repeat % len(selected) :] + selected[: repeat % len(selected)]
                    if repeat % 2:
                        order = list(reversed(order))
                    for b in order:
                        name = b["name"]
                        if name not in contexts:
                            continue
                        a = contexts[name]
                        fn = a["query"]
                        gc.collect()
                        enabled = gc.isenabled()
                        gc.disable()
                        count = 0
                        last = None
                        start = time.perf_counter_ns()
                        cpu_start = cpu_ns()
                        try:
                            while (
                                count == 0
                                or time.perf_counter_ns() - start < args.sample_ms * 1_000_000
                            ):
                                for q in queries:
                                    last = fn(q)
                                    count += 1
                            elapsed = time.perf_counter_ns() - start
                            cpu = cpu_ns() - cpu_start
                        finally:
                            if enabled:
                                gc.enable()
                        assert a["normalize"](last, queries[-1]) == reference[-1]
                        samples[name].append(
                            {"wall_ns": elapsed / count, "cpu_ns": cpu / count, "calls": count}
                        )
                for b in selected:
                    v = samples.get(b["name"], [])
                    if not v:
                        continue
                    result["rows"].append(
                        {
                            "case": case["id"],
                            "interface": case["interface"],
                            "format": case["format"],
                            "workload": case["workload"],
                            "backend": b["name"],
                            "kind": b["kind"],
                            "wall_ns": statistics.median(x["wall_ns"] for x in v),
                            "cpu_ns": statistics.median(x["cpu_ns"] for x in v),
                            "samples": v,
                            "sample_spread": max(x["wall_ns"] for x in v)
                            / min(x["wall_ns"] for x in v),
                            "cache_stats": cache.stats()
                            if b["kind"] == "omics" and cache is not None
                            else None,
                        }
                    )
            atomic_json(args.output, result)
        finally:
            for a in contexts.values():
                a["close"]()
    result["completed"] = True
    result["cli_file_calls"] = {p.path.name: p.calls for p in proxies}
    atomic_json(args.output, result)
    print("complete", len(result["cases"]), len(result["rows"]), len(result["errors"]), flush=True)
    return bool(result["errors"])


if __name__ == "__main__":
    raise SystemExit(main())
