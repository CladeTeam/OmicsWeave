"""Loaded BAM/BCF query APIs against pysam and actual file-producing CLIs."""

from __future__ import annotations

import json
import random
import re
import struct
import subprocess


def canonical(value):
    return json.dumps(value, separators=(",", ":"), ensure_ascii=True).encode()


def queries(width):
    rng = random.Random(20261115 + width)
    return [
        (f"chr{1 + i % 2}", start, start + width)
        for i in range(16)
        for start in [rng.randrange((1 << 22) - width + 1)]
    ]


def aux_bytes(data):
    out, at = [], 0
    scalars = {"c": "b", "C": "B", "s": "h", "S": "H", "i": "i", "I": "I", "f": "f", "d": "d"}
    while at < len(data):
        tag, kind = data[at : at + 2].decode(), chr(data[at + 2])
        at += 3
        if kind in scalars:
            fmt = "<" + scalars[kind]
            value = struct.unpack_from(fmt, data, at)[0]
            at += struct.calcsize(fmt)
            category = "float" if kind in {"f", "d"} else "int"
        elif kind in {"Z", "H"}:
            end = data.index(0, at)
            value = data[at:end].decode()
            at = end + 1
            category = kind
        elif kind == "A":
            value = chr(data[at])
            at += 1
            category = kind
        elif kind == "B":
            subtype = chr(data[at])
            count = struct.unpack_from("<I", data, at + 1)[0]
            at += 5
            fmt = "<" + scalars[subtype] * count
            value = list(struct.unpack_from(fmt, data, at))
            at += struct.calcsize(fmt)
            category = "array"
        else:
            raise ValueError(f"unknown BAM aux type {kind}")
        out.append([tag, category, value])
    return out


def bam_native(raw, all_fields):
    columns = [raw.start, raw.end, raw.flag, raw.mapq, raw.mate_tid, raw.mate_start, raw.tlen]
    rows = []
    for i, row in enumerate(zip(*columns, strict=True)):
        item = [int(v) for v in row]
        if all_fields:
            item += [
                raw.names.bytes(i).decode(),
                raw.seqs.bytes(i).decode(),
                raw.quals.bytes(i).hex(),
                [int(v) for v in raw.cigars[i]],
                aux_bytes(raw.aux.bytes(i)),
            ]
        rows.append(item)
    return canonical(rows)


def bam_pysam(raw, all_fields):
    rows = []
    for value in raw:
        if all_fields:
            r, name, sequence, qualities, cigar, tags = value
        else:
            r = value
            cigar = r.cigartuples or []
        span = sum(n for op, n in cigar if op in {0, 2, 3, 7, 8}) if not r.is_unmapped else 0
        item = [
            r.reference_start,
            r.reference_start + max(span, 1),
            r.flag,
            r.mapping_quality,
            r.next_reference_id,
            r.next_reference_start,
            r.template_length,
        ]
        if all_fields:
            aux = []
            for tag, value, kind in tags:
                category = "int" if kind in "cCsSiI" else "float" if kind in "fd" else kind
                if kind == "B":
                    category, value = "array", list(value)
                aux.append([tag, category, value])
            qualities = bytes(qualities) if qualities is not None else b"\xff" * len(sequence)
            item += [
                name,
                sequence,
                qualities.hex(),
                [(n << 4) | op for op, n in cigar],
                aux,
            ]
        rows.append(item)
    return canonical(rows)


def bam_sam(raw, all_fields):
    refs, rows = {}, []
    for line in raw.decode().splitlines():
        if line.startswith("@SQ\t"):
            name = next(v[3:] for v in line.split("\t") if v.startswith("SN:"))
            refs[name] = len(refs)
        if line.startswith("@") or not line:
            continue
        fields = line.split("\t")
        name, flag, chrom, pos, mapq, cigar, mate, mate_pos, tlen, sequence, qual = fields[:11]
        flag, start = int(flag), int(pos) - 1
        ops = (
            []
            if cigar == "*"
            else [
                ("MIDNSHP=X".index(op), int(n))
                for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
            ]
        )
        span = sum(n for op, n in ops if op in {0, 2, 3, 7, 8}) if not flag & 4 else 0
        item = [
            start,
            start + max(span, 1),
            flag,
            int(mapq),
            refs.get(chrom if mate == "=" else mate, -1),
            int(mate_pos) - 1,
            int(tlen),
        ]
        if all_fields:
            aux = []
            for field in fields[11:]:
                tag, kind, value = field.split(":", 2)
                if kind == "i":
                    category, value = "int", int(value)
                elif kind == "f":
                    category, value = "float", float(value)
                else:
                    category = kind
                aux.append([tag, category, value])
            sequence = "" if sequence == "*" else sequence
            quals = b"\xff" * len(sequence) if qual == "*" else bytes(ord(c) - 33 for c in qual)
            item += [name, sequence, quals.hex(), [(n << 4) | op for op, n in ops], aux]
        rows.append(item)
    return canonical(rows)


def genotype_rows(records, selected):
    out = []
    for record, calls in records:
        out.append(
            [
                record.start,
                record.stop,
                record.id or ".",
                record.ref,
                ",".join(record.alts or []),
                [[-1 if v is None else v for v in (gt or ())] for gt, _ in calls],
                [bool(phased) for _, phased in calls],
            ]
        )
    return canonical([selected, out])


def genotype_native(raw):
    out = []
    for i in range(len(raw)):
        out.append(
            [
                int(raw.start[i]),
                int(raw.end[i]),
                raw.ids.bytes(i).decode() or ".",
                raw.ref.bytes(i).decode(),
                raw.alt.bytes(i).decode(),
                raw.gt[i].tolist(),
                raw.phased[i].tolist(),
            ]
        )
    return canonical([list(raw.samples), out])


def genotype_vcf(raw, selected):
    order, out = [], []
    for line in raw.decode().splitlines():
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            order = [samples.index(s) for s in selected]
        if line.startswith("#") or not line:
            continue
        fields = line.split("\t")
        gt_at = fields[8].split(":").index("GT")
        gts, phases = [], []
        for index in order:
            text = fields[9 + index].split(":")[gt_at]
            gts.append([-1 if v == "." else int(v) for v in re.split(r"[/|]", text)])
            phases.append("|" in text)
        out.append(
            [
                int(fields[1]) - 1,
                int(fields[1]) - 1 + len(fields[3]),
                fields[2] or ".",
                fields[3],
                "" if fields[4] == "." else fields[4],
                gts,
                phases,
            ]
        )
    return canonical([selected, out])


def build_cases(ow, fixtures, tools):
    import pysam

    # Historical snapshots before these query interfaces remain runnable.
    if not hasattr(ow.query, "AlignmentReader") or not hasattr(ow.query, "VariantReader"):
        return []
    required = ["reads.bam", "reads.bam.bai", "reads.bam.csi", "cohort.bcf", "cohort.bcf.csi"]
    for name in required:
        if not (fixtures / name).is_file():
            raise FileNotFoundError(fixtures / name)
    cases = []

    def bam_ours(index, fields):
        def factory():
            reader = ow.query.AlignmentReader(fixtures / "reads.bam", index=index)
            return {
                "query": (lambda q: reader.fetch(q, fields="all")) if fields else reader.fetch,
                "normalize": lambda raw, q: bam_native(raw, fields),
                "close": reader.close,
                "reader": reader,
            }

        return factory

    def bam_binding(index, fields):
        def factory():
            reader = pysam.AlignmentFile(
                str(fixtures / "reads.bam"), "rb", index_filename=str(index)
            )

            def query(q):
                if not fields:
                    return list(reader.fetch(*q))
                return [
                    (
                        r,
                        r.query_name or "",
                        r.query_sequence or "",
                        r.query_qualities,
                        r.cigartuples or [],
                        r.get_tags(with_value_type=True),
                    )
                    for r in reader.fetch(*q)
                ]

            return {
                "query": query,
                "normalize": lambda raw, q: bam_pysam(raw, fields),
                "close": reader.close,
            }

        return factory

    def bam_cli(index, fields):
        def factory():
            def query(q):
                return subprocess.run(
                    [
                        tools["samtools"],
                        "view",
                        "-h",
                        "-X",
                        str(fixtures / "reads.bam"),
                        str(index),
                        f"{q[0]}:{q[1] + 1}-{q[2]}",
                    ],
                    capture_output=True,
                    check=True,
                    timeout=30,
                ).stdout

            return {
                "query": query,
                "normalize": lambda raw, q: bam_sam(raw, fields),
                "close": lambda: None,
            }

        return factory

    samples = [f"S{i:02d}" for i in range(64)]

    def bcf_ours(selected):
        def factory():
            reader = ow.query.VariantReader(fixtures / "cohort.bcf")
            selection = None if selected == samples else selected
            return {
                "query": lambda q: reader.genotypes(q, samples=selection),
                "normalize": lambda raw, q: genotype_native(raw),
                "close": reader.close,
                "reader": reader,
            }

        return factory

    def bcf_binding(selected):
        def factory():
            reader = pysam.VariantFile(str(fixtures / "cohort.bcf"))

            def query(q):
                result = []
                for record in reader.fetch(*q):
                    calls = [record.samples[s] for s in selected]
                    result.append((record, [(c.get("GT"), c.phased) for c in calls]))
                return result

            return {
                "query": query,
                "normalize": lambda raw, q: genotype_rows(raw, selected),
                "close": reader.close,
            }

        return factory

    def bcf_cli(selected):
        def factory():
            def query(q):
                return subprocess.run(
                    [
                        tools["bcftools"],
                        "view",
                        "-Ov",
                        "-r",
                        f"{q[0]}:{q[1] + 1}-{q[2]}",
                        "-s",
                        ",".join(selected),
                        str(fixtures / "cohort.bcf"),
                    ],
                    capture_output=True,
                    check=True,
                    timeout=30,
                ).stdout

            return {
                "query": query,
                "normalize": lambda raw, q: genotype_vcf(raw, selected),
                "close": lambda: None,
            }

        return factory

    for width in [128, 4096, 131072]:
        for index_kind in ["bai", "csi"]:
            index = fixtures / ("reads.bam." + index_kind)
            for fields in [False, True]:
                mode = "all" if fields else "fixed"
                cases.append(
                    {
                        "id": f"alignment.{index_kind}.{mode}.w{width}",
                        "interface": "AlignmentReader.fetch",
                        "format": "bam",
                        "workload": (
                            f"{index_kind}; {mode} fields; width={width}; "
                            "16 alternating-chromosome queries"
                        ),
                        "queries": queries(width),
                        "backends": [
                            {
                                "name": "omics.fetch",
                                "kind": "omics",
                                "factory": bam_ours(index, fields),
                            },
                            {
                                "name": "pysam.AlignmentFile.fetch.list",
                                "kind": "python",
                                "factory": bam_binding(index, fields),
                            },
                            {
                                "name": "samtools.view.SAM",
                                "kind": "cli",
                                "factory": bam_cli(index, fields),
                            },
                        ],
                        "note": (
                            "Binding iterator fully consumed; all mode materializes requested "
                            "sequence/quality/CIGAR/tag properties in timing; canonical conversion "
                            "excluded. CLI uses the same explicit BAI/CSI, full SAM with header."
                        ),
                    }
                )
        for label, selected in [
            ("all64", samples),
            ("one", [samples[0]]),
            ("two_reordered", [samples[-1], samples[0]]),
        ]:
            cases.append(
                {
                    "id": f"genotypes.{label}.w{width}",
                    "interface": "VariantReader.genotypes",
                    "format": "bcf",
                    "workload": (
                        f"64 source samples; selected={len(selected)}; width={width}; "
                        "16 alternating queries"
                    ),
                    "queries": queries(width),
                    "backends": [
                        {"name": "omics.genotypes", "kind": "omics", "factory": bcf_ours(selected)},
                        {
                            "name": "pysam.fetch.selected_GT",
                            "kind": "python",
                            "factory": bcf_binding(selected),
                        },
                        {
                            "name": "bcftools.view.selected_GT",
                            "kind": "cli",
                            "factory": bcf_cli(selected),
                        },
                    ],
                    "note": (
                        "Binding materializes requested GT tuples and phasing; "
                        "NumPy/canonical conversion excluded. CLI writes VCF with "
                        "selected samples and header."
                    ),
                }
            )
    return cases
