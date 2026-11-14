//! htslib as the reference SAM → BAM encoder, this crate's encoder run the
//! same way, and random SAM fixtures.

#![allow(dead_code)]

use std::ffi::{CString, c_void};
use std::path::Path;

use omicsweave_sam::{BamEncoder, BamScratch, BamSpan, SamHeader, bam_header_bytes};

/// splitmix64, for deterministic fixtures from a proptest-chosen seed.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }

    pub fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.range(0, items.len() as u64 - 1) as usize]
    }
}

/// One record as an encoder produced it.
pub type Record = Option<(Vec<u8>, BamSpan)>;

/// A whole SAM text encoded: the BAM header (`None` if the header is
/// rejected) and, per data line, the record (`None` if the line is).
#[derive(Debug)]
pub struct Encoded {
    pub header: Option<Vec<u8>>,
    pub records: Vec<Record>,
}

fn le_i32(b: &[u8], at: usize) -> i32 {
    i32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

/// Reads all of a BGZF file.
fn inflate(path: &Path) -> Vec<u8> {
    let c = |s: &str| CString::new(s).unwrap();
    let mut out = Vec::new();
    // SAFETY: the handle is closed below; reads go into `buf`.
    unsafe {
        let fp = htslib_sys::bgzf_open(c(path.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let mut buf = vec![0u8; 1 << 16];
        loop {
            let n = htslib_sys::bgzf_read(fp, buf.as_mut_ptr() as *mut c_void, buf.len());
            assert!(n >= 0);
            if n == 0 {
                break;
            }
            out.extend_from_slice(&buf[..n as usize]);
        }
        htslib_sys::bgzf_close(fp);
    }
    out
}

/// htslib's encoding: `sam_read1` each record of the SAM text and
/// `sam_write1` it to an uncompressed BAM (`"wbu"`), then cut the inflated
/// stream at the sizes `sam_write1` reports. With `exact`, the stream must
/// hold nothing else (a failed write can leave partial bytes behind).
pub fn htslib(sam: &[u8], exact: bool) -> Encoded {
    let c = |s: &str| CString::new(s).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.sam");
    let output = dir.path().join("out.bam");
    std::fs::write(&input, sam).unwrap();
    let lines = sam.iter().filter(|&&b| b == b'\n').count() + 1;
    let mut written: Vec<Option<(usize, i64)>> = Vec::new();
    // SAFETY: handles come from htslib and are released below.
    unsafe {
        htslib_sys::hts_set_log_level(htslib_sys::HTS_LOG_OFF);
        let fp = htslib_sys::hts_open(c(input.to_str().unwrap()).as_ptr(), c("r").as_ptr());
        assert!(!fp.is_null());
        let hdr = htslib_sys::sam_hdr_read(fp);
        if hdr.is_null() {
            htslib_sys::hts_close(fp);
            return Encoded {
                header: None,
                records: Vec::new(),
            };
        }
        let out = htslib_sys::hts_open(c(output.to_str().unwrap()).as_ptr(), c("wbu").as_ptr());
        assert!(!out.is_null());
        assert_eq!(htslib_sys::sam_hdr_write(out, hdr), 0);
        let b = htslib_sys::bam_init1();
        loop {
            let r = htslib_sys::sam_read1(fp, hdr, b);
            if r == -1 {
                break;
            }
            assert!(written.len() <= lines, "htslib keeps failing at the end");
            if r < -1 {
                written.push(None);
                continue;
            }
            let w = htslib_sys::sam_write1(out, hdr, b);
            written.push((w >= 0).then(|| (w as usize, htslib_sys::bam_endpos(b))));
        }
        htslib_sys::bam_destroy1(b);
        assert_eq!(htslib_sys::hts_close(out), 0);
        htslib_sys::sam_hdr_destroy(hdr);
        htslib_sys::hts_close(fp);
    }
    let stream = inflate(&output);
    // BAM header: magic, l_text, text, n_ref, then l_name, name, l_ref each.
    let mut at = 8 + le_i32(&stream, 4) as usize;
    let n_ref = le_i32(&stream, at);
    at += 4;
    for _ in 0..n_ref {
        at += 8 + le_i32(&stream, at) as usize;
    }
    let header = stream[..at].to_vec();
    let records = written
        .into_iter()
        .map(|w| {
            let (size, end) = w?;
            let rec = stream[at..at + size].to_vec();
            at += size;
            let flag = u16::from_le_bytes([rec[18], rec[19]]);
            let span = BamSpan {
                tid: le_i32(&rec, 4),
                beg: i64::from(le_i32(&rec, 8)),
                end,
                mapped: flag & 4 == 0,
            };
            Some((rec, span))
        })
        .collect();
    if exact {
        assert_eq!(at, stream.len(), "bytes after the last record");
    }
    Encoded {
        header: Some(header),
        records,
    }
}

/// The data lines after the header, as a line splitter hands them over:
/// without their `\n`, and no empty line after a final `\n`.
pub fn data_lines(rest: &[u8]) -> Vec<&[u8]> {
    if rest.is_empty() {
        return Vec::new();
    }
    let body = rest.strip_suffix(b"\n").unwrap_or(rest);
    body.split(|&c| c == b'\n').collect()
}

/// This crate's encoding, appending every record to one buffer as a
/// worker does.
pub fn ours(sam: &[u8]) -> Encoded {
    let mut input = sam;
    let Ok(header) = SamHeader::read(&mut input) else {
        return Encoded {
            header: None,
            records: Vec::new(),
        };
    };
    assert_eq!(header.raw(), &sam[..sam.len() - input.len()]);
    let encoder = BamEncoder::new(&header);
    let mut scratch = BamScratch::new();
    let mut all = b"prefix".to_vec();
    let mut records = Vec::new();
    for (i, line) in data_lines(input).into_iter().enumerate() {
        let before = all.len();
        let line_no = header.lines() + 1 + i as u64;
        match encoder.encode_line(&mut scratch, line, line_no, &mut all) {
            Ok(span) => {
                let rec = all[before..].to_vec();
                assert_eq!(
                    u32::from_le_bytes(rec[..4].try_into().unwrap()) as usize,
                    rec.len() - 4
                );
                records.push(Some((rec, span)));
            }
            Err(e) => {
                assert_eq!(all.len(), before, "a failed line left bytes behind");
                assert!(
                    matches!(e, omicsweave_core::Error::Format { line, .. } if line == line_no)
                );
                records.push(None);
            }
        }
    }
    assert!(all.starts_with(b"prefix"));
    Encoded {
        header: Some(bam_header_bytes(&header)),
        records,
    }
}

/// A random SAM text and, per data line, whether it names a reference the
/// header lacks (which this crate rejects and htslib does not).
pub struct Fixture {
    pub text: Vec<u8>,
    pub unknown: Vec<bool>,
}

const NAMES: &[&str] = &[
    "chr1",
    "chr2",
    "chrX",
    "chrM",
    "scaffold_17",
    "HLA-A*01:01:01:01",
    "1",
    "chrUn_KI270302v1",
];
const ALT_NAMES: &[&str] = &["alt1", "alt2", "AltX", "NC_000001.11"];
const LENGTHS: &[&str] = &[
    "1",
    "1000",
    "248956422",
    "536870912",
    "2147483647",
    "4294967295",
    "4294967296",
    "99999999999999999999",
    " 500",
    "+500",
    "500abc",
];
const POSITIONS: &[&str] = &[
    "0",
    "1",
    "2",
    "16384",
    "16385",
    "131072",
    "536870912",
    "536870913",
    "2147483647",
];
const INTS: &[&str] = &[
    "-2147483649",
    "-2147483648",
    "-32769",
    "-32768",
    "-129",
    "-128",
    "-1",
    "-0",
    "0",
    "1",
    "127",
    "128",
    "255",
    "256",
    "32767",
    "32768",
    "65535",
    "65536",
    "2147483647",
    "2147483648",
    "4294967295",
    "4294967296",
    "+12",
    "007",
];
const FLOATS: &[&str] = &[
    "0",
    "1.5",
    "-0",
    "-2.75",
    "0.1",
    "1e10",
    "3.4028235e38",
    "3.4028236e38",
    "1e39",
    "-1e39",
    "1e-46",
    "1.4e-45",
    "nan",
    "-nan",
    "inf",
    "-Infinity",
    "0x1p-3",
    "0x1.8p1",
    "1.",
    ".5",
    " 2.5",
    "5e",
    "12abc",
    "123456789.123456789",
];
const BAD_AUX: &[&str] = &[
    "XX:i:",
    "XQ:Q:1",
    "X",
    "XX",
    "XX:Z",
    "XX:A:",
    "XX:f:",
    "XX:f:abc",
    "XX:A:ab",
    "XX:i:12abc",
    "XYZi:5",
    "XB:B:c,1x",
    "XB:B:q,1",
    "XB:B:d,1",
    "XB:B:Z,1",
    "XH:H:ABC",
    "XB:B:A",
    "XB:B:A,1",
    "XB:B",
    "XB:B:C,-",
    "XB:B:c,-",
    "XB:B:C,-1,5",
    "XB:B:I,-5",
    "XB:B:S,70000",
    "XB:B:c,9223372036854775808",
    "XB:B:i,5,,7,",
    " X:i:1",
    "\u{e9}X:i:1",
    "XA:A:\u{e9}",
    "XI:i:5\u{e9}x",
    "XI:i:5\u{e9}xyz:i:7",
    "XB:B:C,-\u{e9},1",
];

/// B arrays whose values overflow their subtype, which htslib re-types
/// from the values' range (note `max < 255` for `C`).
const RETYPED: &[&str] = &[
    "XB:B:c,255",
    "XB:B:c,254,0",
    "XB:B:c,128,1",
    "XB:B:C,-1,127",
    "XB:B:C,-128,127,0",
    "XB:B:C,-129",
    "XB:B:s,65535",
    "XB:B:s,32768,-1",
    "XB:B:S,-32768,32767",
    "XB:B:S,65536",
    "XB:B:i,4294967295",
    "XB:B:i,2147483648,-1",
    "XB:B:I,-2147483648",
    "XB:B:I,4294967296",
    "XB:B:c,-9223372036854775808",
];

fn word(rng: &mut Rng, alphabet: &[u8], lo: u64, hi: u64) -> Vec<u8> {
    let n = rng.range(lo, hi);
    (0..n).map(|_| *rng.pick(alphabet)).collect()
}

fn header(rng: &mut Rng, names: &mut Vec<String>) -> Vec<Vec<u8>> {
    let mut lines: Vec<Vec<u8>> = Vec::new();
    let n_ref = rng.range(0, 4) as usize;
    let mut alts = ALT_NAMES.to_vec();
    for _ in 0..n_ref {
        let name = loop {
            let n = rng.pick(NAMES).to_string();
            if !names.contains(&n) {
                break n;
            }
        };
        let mut tags = vec![format!("SN:{name}"), format!("LN:{}", rng.pick(LENGTHS))];
        if rng.chance(20) {
            tags.swap(0, 1);
        }
        if rng.chance(20) {
            tags.push("M5:6aef897c3d6ff0c78aff06ac189178dd".into());
        }
        if rng.chance(20) && !alts.is_empty() {
            let k = rng.range(0, alts.len() as u64 - 1) as usize;
            let alt = alts.remove(k);
            tags.push(format!("AN:{alt},,{name}"));
            names.push(alt.to_string());
        }
        names.push(name);
        lines.push(format!("@SQ\t{}", tags.join("\t")).into_bytes());
    }
    let others: [(u64, &[u8]); 6] = [
        (40, b"@HD\tVN:1.6\tSO:coordinate"),
        (10, b"@HD\tVN:1.4"),
        (30, b"@RG\tID:rg1\tSM:sample one\tPL:ILLUMINA"),
        (
            30,
            b"@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem ref.fa r1.fq",
        ),
        (20, b"@PG\tID:samtools\tPN:samtools\tPP:bwa"),
        (30, b"@CO\tfree text\twith a tab"),
    ];
    for (percent, line) in others {
        if rng.chance(percent) {
            let at = rng.range(0, lines.len() as u64) as usize;
            lines.insert(at, line.to_vec());
        }
    }
    if !lines.is_empty() && rng.chance(10) {
        let at = rng.range(1, lines.len() as u64) as usize;
        lines.insert(at, b"@CO".to_vec());
    }
    if !lines.is_empty() && rng.chance(5) {
        let bad: &[&[u8]] = &[
            b"@SQ\tSN:nolength",
            b"@SQ\tSN:zero\tLN:0",
            b"@SQ\tSN:neg\tLN:-5",
            b"@SQ\tLN:10",
            b"@SQ\tSN:two\tLN:10\tLN:11",
            b"@SQ\tSN:same\tLN:10\tLN:10",
            b"@SQ\tSN:chr1\tLN:5\t",
            b"@SQ\tSNchr1\tLN:5",
            b"@RG\tSM:x",
            b"@PG\tPN:x",
            b"@XY\tAB:c",
            b"@HD",
            b"@Co\tx",
            b"@COx",
            b"@",
            b"@1Q\tSN:x",
        ];
        let at = rng.range(1, lines.len() as u64) as usize;
        lines.insert(at, rng.pick(bad).to_vec());
        if rng.chance(30) && !names.is_empty() {
            let dup = format!("@SQ\tSN:{}\tLN:9", names[0]);
            lines.push(dup.into_bytes());
        }
    }
    lines
}

fn number_field(rng: &mut Rng, plain: bool) -> String {
    match rng.range(0, 9) {
        0..=3 => rng.range(0, 3_000_000).to_string(),
        4..=7 => rng.pick(POSITIONS).to_string(),
        _ if plain => "100".into(),
        _ => rng
            .pick(&["", "+7", "2147483648", "4611686018427387904", "-1", "7x"])
            .to_string(),
    }
}

const QUERY_OPS: &[u8] = b"MIS=X";
const REF_OPS: &[u8] = b"DNPH";

/// A CIGAR consuming exactly `len` query bases, as (length, op) pairs.
fn cigar_ops(rng: &mut Rng, len: u64, max_ops: u64) -> Vec<(u64, u8)> {
    let mut ops = Vec::new();
    let mut left = len;
    let n = rng.range(1, max_ops);
    for i in 0..n {
        if left > 0 && (i == n - 1 || rng.chance(70)) {
            let l = if i == n - 1 { left } else { rng.range(1, left) };
            ops.push((l, *rng.pick(QUERY_OPS)));
            left -= l;
        } else {
            ops.push((rng.range(0, 60), *rng.pick(REF_OPS)));
        }
    }
    if left > 0 {
        ops.push((left, b'M'));
    }
    ops
}

fn cigar_text(ops: &[(u64, u8)]) -> Vec<u8> {
    let mut s = Vec::new();
    for &(l, op) in ops {
        s.extend_from_slice(l.to_string().as_bytes());
        s.push(op);
    }
    s
}

fn cigar_code(op: u8) -> u32 {
    b"MIDNSHP=XB".iter().position(|&c| c == op).unwrap() as u32
}

fn aux(rng: &mut Rng, plain: bool) -> Vec<u8> {
    let id = rng.pick(&["XA", "XB", "NM", "MD", "AS", "RG", "XS", "ZZ", "a1", "tp"]);
    let mut f = format!("{id}:").into_bytes();
    match rng.range(0, if plain { 9 } else { 11 }) {
        0 => {
            f.extend_from_slice(b"A:");
            f.push(*rng.pick(b"AcZ=!~0"));
        }
        1 | 2 => {
            let v = if rng.chance(70) {
                let v = rng.pick(INTS);
                if plain && matches!(*v, "-2147483649" | "4294967296") {
                    "-7".to_string()
                } else {
                    v.to_string()
                }
            } else {
                (rng.range(0, 10_000_000_000) as i64 - 5_000_000_000).to_string()
            };
            f.extend_from_slice(format!("{}:{v}", rng.pick(&["i", "I"])).as_bytes());
        }
        3 => f.extend_from_slice(format!("f:{}", rng.pick(FLOATS)).as_bytes()),
        4 | 5 => {
            f.extend_from_slice(b"Z:");
            f.extend(word(rng, b"ACGT acgt:;=,.-_0123456789", 0, 20));
            if !plain && rng.chance(10) {
                f.extend_from_slice("café".as_bytes());
            }
        }
        6 => {
            f.extend_from_slice(b"H:");
            let n = 2 * rng.range(0, 6);
            f.extend(word(rng, b"0123456789ABCDEFabcdef", n, n));
        }
        7 | 8 if rng.chance(15) => f = rng.pick(RETYPED).as_bytes().to_vec(),
        7 | 8 => {
            let sub = *rng.pick(b"cCsSiIf");
            f.extend_from_slice(format!("B:{}", sub as char).as_bytes());
            let n = rng.range(0, 8);
            for _ in 0..n {
                let v = if sub == b'f' {
                    let v = rng.pick(FLOATS).trim();
                    if plain && matches!(v, "5e" | "12abc") {
                        "2".to_string()
                    } else {
                        v.to_string()
                    }
                } else if plain || rng.chance(70) {
                    let (lo, hi): (i64, i64) = match sub {
                        b'c' => (-128, 127),
                        b'C' => (0, 255),
                        b's' => (-32768, 32767),
                        b'S' => (0, 65535),
                        b'i' => (-2147483648, 2147483647),
                        _ => (0, 4294967295),
                    };
                    let edges = [lo, hi, 0, lo + 1, hi - 1];
                    if rng.chance(50) {
                        rng.pick(&edges).to_string()
                    } else {
                        (lo + (rng.next() % (hi - lo + 1) as u64) as i64).to_string()
                    }
                } else {
                    rng.pick(INTS).to_string()
                };
                f.push(b',');
                f.extend_from_slice(v.as_bytes());
            }
            if !plain && rng.chance(5) {
                f.push(b',');
            }
        }
        9 => f.extend_from_slice(format!("d:{}", rng.pick(FLOATS)).as_bytes()),
        _ => f = rng.pick(BAD_AUX).as_bytes().to_vec(),
    }
    f
}

/// One data line. `plain` lines stay clear of deliberately malformed
/// fields, and so are what htslib's format detection recognises as SAM, as
/// the first line of a headerless file must be.
fn record(rng: &mut Rng, names: &[String], plain: bool) -> (Vec<u8>, bool) {
    let plain = plain || rng.chance(60);
    let mut unknown = false;
    let mut f: Vec<Vec<u8>> = Vec::new();
    // QNAME
    f.push(if plain || rng.chance(95) {
        let mut q = word(
            rng,
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
            1,
            1,
        );
        q.extend(word(rng, b"ABCabc0123456789_:.-/#", 0, 24));
        q
    } else {
        rng.pick(&[String::new(), "q".repeat(254), "q".repeat(255)])
            .as_bytes()
            .to_vec()
    });
    // FLAG
    let flag = *rng.pick(&[
        0u32, 4, 16, 83, 99, 147, 163, 256, 2048, 0x900, 0xffff, 1, 8, 0x40,
    ]);
    let flag = if rng.chance(20) {
        rng.range(0, 65535) as u32
    } else {
        flag
    };
    f.push(
        match rng.range(0, if plain { 0 } else { 19 }) {
            0..=13 => flag.to_string(),
            14 => format!("0x{flag:x}"),
            15 => format!("0{flag:o}"),
            _ => rng
                .pick(&["", "65536", "0x10000", "1a", "-1", "0x", "08"])
                .to_string(),
        }
        .into_bytes(),
    );
    // RNAME
    let rname = if !plain && rng.chance(3) {
        unknown = true;
        "chrUnknown".to_string()
    } else if !names.is_empty() && rng.chance(85) {
        rng.pick(names).clone()
    } else {
        "*".to_string()
    };
    let placed = rname != "*";
    f.push(rname.into_bytes());
    // POS
    let pos = number_field(rng, plain);
    f.push(pos.clone().into_bytes());
    // MAPQ
    f.push(
        match rng.range(0, if plain { 3 } else { 5 }) {
            0..=3 => rng.range(0, 255).to_string(),
            _ => rng.pick(&["", "256", "60x"]).to_string(),
        }
        .into_bytes(),
    );
    // CIGAR, SEQ, QUAL
    let len = if rng.chance(10) { 0 } else { rng.range(1, 40) };
    let mut cg: Option<Vec<u8>> = None;
    let cigar = match rng.range(0, 99) {
        0..=14 => b"*".to_vec(),
        15..=79 => cigar_text(&cigar_ops(rng, len, 6)),
        80..=84 => {
            let wrong = len + rng.range(1, 3);
            cigar_text(&cigar_ops(rng, wrong, 4))
        }
        85..=91 if len > 0 && placed && pos != "0" => {
            // A CIGAR stored in a CG tag behind a soft-clip placeholder.
            let real = cigar_ops(rng, len, 8);
            let mut tag = format!("CG:B:{},", rng.pick(&["I", "I", "i"])).into_bytes();
            let codes: Vec<String> = real
                .iter()
                .map(|&(l, op)| ((l as u32) << 4 | cigar_code(op)).to_string())
                .collect();
            tag.extend_from_slice(codes.join(",").as_bytes());
            cg = Some(tag);
            if rng.chance(70) {
                format!("{len}S").into_bytes()
            } else {
                format!("{len}S{}N", rng.range(1, 100)).into_bytes()
            }
        }
        _ if plain => cigar_text(&cigar_ops(rng, len, 6)),
        _ => rng
            .pick(&[
                "10M5",
                "M",
                "0M",
                "268435455N",
                "268435456M",
                "10m",
                "4B",
                "+4M",
            ])
            .as_bytes()
            .to_vec(),
    };
    f.push(cigar);
    // RNEXT
    f.push(if !plain && rng.chance(2) {
        unknown = true;
        b"chrMate".to_vec()
    } else {
        match rng.range(0, 2) {
            0 => b"=".to_vec(),
            1 => b"*".to_vec(),
            _ if names.is_empty() => b"*".to_vec(),
            _ => rng.pick(names).clone().into_bytes(),
        }
    });
    // PNEXT, TLEN
    f.push(number_field(rng, plain).into_bytes());
    let tlens: &[&str] = &["0", "-1", "150", "-150", "+5", "2147483647", "-2147483648"];
    let wild: &[&str] = &[
        "2147483648",
        "-2147483649",
        "4611686018427387903",
        "-4611686018427387904",
        "4611686018427387904",
    ];
    let tlens = if plain || rng.chance(50) { tlens } else { wild };
    f.push(rng.pick(tlens).as_bytes().to_vec());
    let seq = if len == 0 && rng.chance(80) {
        b"*".to_vec()
    } else {
        word(rng, b"ACGTNacgtn=MRWSYKVHDB.*", len, len)
    };
    let l_seq = if seq == b"*" { 0 } else { seq.len() as u64 };
    f.push(seq);
    f.push(match rng.range(0, if plain { 89 } else { 99 }) {
        0..=79 => word(rng, &(33u8..=126).collect::<Vec<_>>(), l_seq, l_seq),
        80..=89 => b"*".to_vec(),
        90..=92 => word(rng, b"I", l_seq + 1, l_seq + 1),
        93 => vec![b' '; l_seq as usize],
        94 => vec![0xa0; l_seq as usize],
        95 => vec![0xa1; l_seq as usize],
        _ => Vec::new(),
    });
    for _ in 0..rng.range(0, 6) {
        f.push(aux(rng, plain));
    }
    if let Some(tag) = cg {
        let at = rng.range(11, f.len() as u64) as usize;
        f.insert(at, tag);
    }
    (f.join(&b'\t'), unknown)
}

/// A random SAM text: possibly no header, references with alternative
/// names, all record fields at their boundaries, and malformed fields.
pub fn random_sam(rng: &mut Rng) -> Fixture {
    let mut names = Vec::new();
    let header = if rng.chance(15) {
        Vec::new()
    } else {
        header(rng, &mut names)
    };
    let eol: &[u8] = if rng.chance(20) { b"\r\n" } else { b"\n" };
    let mut text = Vec::new();
    for l in &header {
        text.extend_from_slice(l);
        text.extend_from_slice(eol);
    }
    let mut unknown = Vec::new();
    let n = rng.range(if header.is_empty() { 1 } else { 0 }, 12);
    for i in 0..n {
        let (line, u) = record(rng, &names, header.is_empty() && i == 0);
        text.extend_from_slice(&line);
        text.extend_from_slice(eol);
        unknown.push(u);
    }
    Fixture { text, unknown }
}

/// SAM with records whose CIGARs exceed BAM's 65535 operations, given as
/// CIGAR text or as a CG tag behind a soft-clip placeholder.
pub fn long_cigar_sam(rng: &mut Rng) -> Fixture {
    let mut text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n".to_vec();
    let mut unknown = Vec::new();
    for _ in 0..rng.range(1, 3) {
        let n_ops = rng.range(65_530, 66_600);
        let mut ops: Vec<(u64, u8)> = (0..n_ops)
            .map(|i| {
                let op = if i % 2 == 0 {
                    b'M'
                } else {
                    *rng.pick(b"IDN=X")
                };
                (rng.range(1, 3), op)
            })
            .collect();
        if rng.chance(30) {
            ops.insert(0, (rng.range(1, 9), b'S'));
        }
        let qlen: u64 = ops
            .iter()
            .filter(|(_, op)| b"MIS=X".contains(op))
            .map(|(l, _)| l)
            .sum();
        let flag = *rng.pick(&[0, 16, 4, 99]);
        let seq = word(rng, b"ACGT", qlen, qlen);
        let qual = if rng.chance(50) {
            b"*".to_vec()
        } else {
            word(rng, b"#5?I", qlen, qlen)
        };
        let pos = rng.pick(&["1", "1000", "100000000"]).to_string();
        let mut line = format!("long{}\t{flag}\tchr1\t{pos}\t60\t", rng.range(0, 999)).into_bytes();
        let as_tag = rng.chance(40);
        if as_tag {
            line.extend_from_slice(format!("{qlen}S").as_bytes());
        } else {
            line.extend_from_slice(&cigar_text(&ops));
        }
        line.extend_from_slice(b"\t=\t1\t0\t");
        line.extend_from_slice(&seq);
        line.push(b'\t');
        line.extend_from_slice(&qual);
        line.extend_from_slice(b"\tNM:i:3");
        if as_tag {
            line.extend_from_slice(b"\tCG:B:I");
            for &(l, op) in &ops {
                line.extend_from_slice(format!(",{}", (l as u32) << 4 | cigar_code(op)).as_bytes());
            }
        }
        line.extend_from_slice(b"\tXZ:Z:after\n");
        text.extend_from_slice(&line);
        unknown.push(false);
    }
    Fixture { text, unknown }
}
