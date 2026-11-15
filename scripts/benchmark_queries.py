"""Compare public OmicsWeave query APIs with real bindings and CLI tools.

Run in the isolated comparator uv environment, with candidate src on PYTHONPATH.
All fixture/index generation is separate. Every result is checked before timing.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import platform
import resource
import statistics
import sys
import time
import traceback
from pathlib import Path


def cpu_ns():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time_ns() + round((r.ru_utime + r.ru_stime) * 1e9)


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False))
    temp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--tools", type=Path, default=Path("/tmp/omicsweave-final-r8-tools"))
    p.add_argument("--adapters", type=Path, default=Path(__file__).parent / "query_comparators")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--validate-only", action="store_true")
    p.add_argument("--filter", default="")
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--sample-ms", type=int, default=100)
    p.add_argument("--cpu", type=int)
    args = p.parse_args()
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    import importlib.metadata as md

    import omicsweave as ow

    tools = {
        n: str(args.tools / n) for n in ["twoBitToFa", "bigWigToBedGraph", "bigBedToBed", "straw"]
    }
    import shutil

    for n in ["samtools", "bcftools", "bedtools"]:
        tools[n] = shutil.which(n)
    cases = []
    for name in ["sequence", "track_interval", "consensus", "contact"]:
        path = args.adapters / (name + ".py")
        spec = importlib.util.spec_from_file_location("query_compare_" + name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        cases.extend(module.build_cases(ow, args.fixtures, tools))
    cases = [c for c in cases if args.filter in c["id"]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "protocol": {
            "repeats": args.repeats,
            "sample_ms": args.sample_ms,
            "cache": "warm OS page cache; no global cache flushing",
            "raw_reuse": "opened reader; raw native fully materialized output",
            "canonical_reuse": "opened reader; equivalent normalized complete payload",
            "canonical_open_query_close": "open + query + normalize + close per region",
            "cli": "subprocess startup/capture included; no persistent CLI reader",
            "cpu_time": "self plus waited child process CPU time",
            "cpu_affinity": sorted(os.sched_getaffinity(0)),
            "cpu_exclusive": False,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "omics_native": ow._native.__file__,
            "omics_native_sha256": hashlib.sha256(
                Path(ow._native.__file__).read_bytes()
            ).hexdigest(),
            "versions": {
                n: md.version(n)
                for n in [
                    "numpy",
                    "pysam",
                    "pyfaidx",
                    "py2bit",
                    "pyBigWig",
                    "hic-straw",
                    "ncls",
                    "pyranges",
                    "biopython",
                ]
            },
        },
        "fixtures": str(args.fixtures),
        "cases": [],
        "rows": [],
        "errors": [],
    }
    atomic_json(args.output, result)
    for case in cases:
        print("validate", case["id"], flush=True)
        qs = case["queries"]
        reference = None
        valid = []
        case_info = {k: v for k, v in case.items() if k not in {"backends", "queries"}}
        case_info["queries"] = qs
        case_info["backends"] = []
        for backend in case["backends"]:
            info = {k: v for k, v in backend.items() if k != "factory"}
            obj = None
            try:
                t = time.perf_counter_ns()
                obj = backend["factory"]()
                info["open_once_ns"] = time.perf_counter_ns() - t
                hashes = []
                sizes = []
                for i, q in enumerate(qs):
                    raw = obj["query"](q)
                    data = obj["normalize"](raw, q)
                    if not isinstance(data, bytes):
                        raise TypeError("normalizer must materialize bytes")
                    hashes.append(hashlib.sha256(data).hexdigest())
                    sizes.append(len(data))
                    info["raw_type"] = type(raw).__module__ + "." + type(raw).__qualname__
                    if reference is not None and data != reference[i]:
                        a, b = data, reference[i]
                        j = next(
                            (j for j, (x, y) in enumerate(zip(a, b, strict=False)) if x != y),
                            min(len(a), len(b)),
                        )
                        raise AssertionError(
                            f"payload mismatch query={q} bytes={len(a)}/{len(b)} first_difference={j}"
                        )
                    if reference is None:
                        info.setdefault("_payloads", []).append(data)
                if reference is None:
                    reference = info.pop("_payloads")
                info["checksums"] = hashes
                info["payload_bytes"] = sizes
                info["status"] = "ok"
                valid.append(backend)
            except Exception as e:
                info.pop("_payloads", None)
                info.update(status="error", error=repr(e))
                result["errors"].append(
                    {
                        "case": case["id"],
                        "backend": backend["name"],
                        "stage": "validation",
                        "error": repr(e),
                        "traceback": traceback.format_exc(),
                    }
                )
                print("validation error", case["id"], backend["name"], repr(e), flush=True)
            finally:
                if obj is not None:
                    obj["close"]()
            case_info["backends"].append(info)
        result["cases"].append(case_info)
        atomic_json(args.output, result)
        if args.validate_only or case.get("correctness_only"):
            continue
        for mode in ["raw_reuse", "canonical_reuse", "canonical_open_query_close"]:
            contexts = {}
            samples = {}
            opens = {}
            for b in valid:
                try:
                    if mode != "canonical_open_query_close":
                        t = time.perf_counter_ns()
                        a = b["factory"]()
                        opens[b["name"]] = time.perf_counter_ns() - t
                        contexts[b["name"]] = a
                        a["normalize"](a["query"](qs[0]), qs[0])
                    samples[b["name"]] = []
                except Exception as e:
                    result["errors"].append(
                        {"case": case["id"], "backend": b["name"], "stage": mode, "error": repr(e)}
                    )
            try:
                for repeat in range(args.repeats):
                    ordered = (
                        valid[repeat % len(valid) :] + valid[: repeat % len(valid)] if valid else []
                    )
                    if repeat % 2:
                        ordered = list(reversed(ordered))
                    for b in ordered:
                        name = b["name"]
                        if name not in samples:
                            continue
                        factory = b["factory"]
                        if mode == "canonical_open_query_close":

                            def op(q, factory=factory):
                                a = factory()
                                try:
                                    return a["normalize"](a["query"](q), q)
                                finally:
                                    a["close"]()
                        elif mode == "raw_reuse":
                            op = contexts[name]["query"]
                        else:
                            a = contexts[name]

                            def op(q, a=a):
                                return a["normalize"](a["query"](q), q)

                        gc.collect()
                        was_enabled = gc.isenabled()
                        gc.disable()
                        count = 0
                        last = None
                        begin = time.perf_counter_ns()
                        cpu_begin = cpu_ns()
                        try:
                            while (
                                count == 0
                                or time.perf_counter_ns() - begin < args.sample_ms * 1_000_000
                            ):
                                for q in qs:
                                    last = op(q)
                                    count += 1
                            cpu = cpu_ns() - cpu_begin
                            wall = time.perf_counter_ns() - begin
                        finally:
                            if was_enabled:
                                gc.enable()
                        samples[name].append(
                            {
                                "wall_ns_per_query": wall / count,
                                "cpu_ns_per_query": cpu / count,
                                "queries": count,
                            }
                        )
                        # Keep the returned object alive through completion; timing already ended.
                        if last is None:
                            raise ValueError("query returned no payload")
                        checked = (
                            contexts[name]["normalize"](last, qs[-1])
                            if mode == "raw_reuse"
                            else last
                        )
                        if checked != reference[-1]:
                            raise AssertionError(
                                f"repeated query payload changed: {case['id']} {name} {mode}"
                            )
                for b in valid:
                    values = samples.get(b["name"], [])
                    if not values:
                        continue
                    row = {
                        "case": case["id"],
                        "interface": case["interface"],
                        "format": case["format"],
                        "workload": case["workload"],
                        "backend": b["name"],
                        "kind": b["kind"],
                        "mode": mode,
                        "samples": values,
                        "wall_ns": statistics.median(v["wall_ns_per_query"] for v in values),
                        "cpu_ns": statistics.median(v["cpu_ns_per_query"] for v in values),
                        "open_once_ns": opens.get(b["name"]),
                        "note": b.get("note", case.get("note", "")),
                    }
                    row["wall_spread"] = max(v["wall_ns_per_query"] for v in values) / min(
                        v["wall_ns_per_query"] for v in values
                    )
                    result["rows"].append(row)
                print("measured", case["id"], mode, flush=True)
                atomic_json(args.output, result)
            finally:
                for a in contexts.values():
                    a["close"]()
    result["completed"] = True
    atomic_json(args.output, result)
    print(
        "complete",
        len(result["cases"]),
        "cases",
        len(result["rows"]),
        "rows",
        len(result["errors"]),
        "errors",
        flush=True,
    )
    return 1 if result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
