import re

CHROMS = {"chr1": 5000, "chr2": 3000}


RECORDS = [
    ("a", 99, "chr1", 101, 60, "10M", "=", 251, 160, "ACGTACGTAC", "IIIIIIIIII", "NM:i:0"),
    ("b", 0, "chr1", 105, 30, "3S5M2D4M", "*", 0, 0, "GGGACGTAAAAA", "*", ""),
    ("c", 16, "chr1", 150, 10, "4M100N4M", "*", 0, 0, "ACGTACGT", "!!!!####", "XS:A:-"),
    ("d", 1024, "chr1", 200, 60, "8M", "*", 0, 0, "ACGTACGT", "*", ""),
    ("u", 69, "chr1", 220, 0, "*", "=", 101, 0, "ACGT", "*", ""),
    ("a", 147, "chr1", 251, 60, "10M", "=", 101, -160, "ACGTACGTAC", "IIIIIIIIII", ""),
    ("e", 0, "chr1", 4995, 60, "6M", "*", 0, 0, "ACGTAC", "*", ""),
    ("f", 0, "chr2", 10, 60, "5M", "*", 0, 0, "ACGTA", "*", ""),
    ("x", 4, "*", 0, 0, "*", "*", 0, 0, "ACGT", "*", ""),
]


def sam_text() -> str:
    head = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
        f"@SQ\tSN:{n}\tLN:{ln}\n" for n, ln in CHROMS.items()
    )
    body = ""
    for r in RECORDS:
        fields = [str(x) for x in r[:11]]
        if r[11]:
            fields.append(r[11])
        body += "\t".join(fields) + "\n"
    return head + body


def ref_len(cigar: str) -> int:
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X")


def span(r) -> tuple[int, int]:
    start = r[3] - 1
    rlen = 0 if r[1] & 4 or r[5] == "*" else ref_len(r[5])
    return start, start + (rlen or 1)


def expected(chrom: str, start: int, end: int):
    return [r for r in RECORDS if r[2] == chrom and span(r)[1] > start and span(r)[0] < end]
