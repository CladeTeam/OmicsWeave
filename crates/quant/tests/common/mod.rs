//! Shared helpers of the quantification tests: random sorted reads and
//! gene models, BAMs written by htslib, the reference tools run in
//! process (samtools depth, featureCounts), and readers of our outputs.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::ffi::{CString, c_char, c_int};
use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use htslib_sys as h;
use omicsweave_core::{FileSource, TrackRead};

/// The reference tools keep process-global state (getopt, statics).
pub static TOOLS: Mutex<()> = Mutex::new(());

/// A small deterministic generator (splitmix64).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `lo..hi`.
    pub fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo).max(1)
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }

    pub fn pick<'a, T>(&mut self, v: &'a [T]) -> &'a T {
        &v[self.range(0, v.len() as u64) as usize]
    }
}

/// One alignment record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    pub flag: u16,
    pub tid: usize,
    /// 0-based.
    pub pos: u32,
    pub mapq: u8,
    /// `(len, op)`; empty for `*`.
    pub cigar: Vec<(u32, char)>,
    pub mtid: i32,
    pub mpos: u32,
    pub tlen: i64,
    pub tags: Vec<String>,
    /// Bases and qualities (Phred + 33) as written: random bases and
    /// qualities without `seq`, no qualities (`*`) with `seq` but not
    /// `qual`.
    pub seq: Option<Vec<u8>>,
    pub qual: Option<Vec<u8>>,
}

impl Read {
    pub fn ref_len(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|(_, o)| matches!(o, 'M' | 'D' | 'N' | '=' | 'X'))
            .map(|(l, _)| l)
            .sum()
    }

    pub fn query_len(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|(_, o)| matches!(o, 'M' | 'I' | 'S' | '=' | 'X'))
            .map(|(l, _)| l)
            .sum()
    }

    /// `bam_endpos`.
    pub fn end(&self) -> u32 {
        let r = if self.flag & 4 != 0 {
            0
        } else {
            self.ref_len()
        };
        self.pos + r.max(1)
    }

    /// Reference intervals of `M`/`=`/`X` (and `D` with `deletions`).
    pub fn blocks(&self, deletions: bool) -> Vec<(u32, u32)> {
        let mut out: Vec<(u32, u32)> = Vec::new();
        let mut at = self.pos;
        for &(l, o) in &self.cigar {
            match o {
                'M' | '=' | 'X' => {
                    out.push((at, at + l));
                    at += l;
                }
                'D' if deletions => {
                    out.push((at, at + l));
                    at += l;
                }
                'D' | 'N' => at += l,
                _ => {}
            }
        }
        out
    }

    pub fn tag(&self, key: &str) -> Option<&str> {
        self.tags.iter().find_map(|t| t.strip_prefix(key)?.get(3..))
    }

    pub fn nh(&self) -> u32 {
        self.tag("NH").map_or(1, |v| v.parse().unwrap())
    }

    pub fn reverse(&self) -> bool {
        self.flag & 0x10 != 0
    }

    /// The leftmost mate of a pair mapped with its mate on one reference.
    pub fn leader(&self) -> bool {
        self.flag & 1 != 0
            && self.flag & 0xc == 0
            && self.mtid == self.tid as i32
            && self.tlen != 0
            && (self.pos < self.mpos || (self.pos == self.mpos && self.flag & 0x40 != 0))
    }
}

/// What the generator makes.
#[derive(Debug, Clone)]
pub struct Spec {
    pub refs: Vec<(String, u32)>,
    pub fragments: usize,
    pub paired: bool,
    pub spliced: bool,
    pub barcodes: usize,
    pub umis: usize,
    /// Duplicates, secondary / supplementary / QC-failed records, low MAPQ,
    /// multi-mappers, unmapped mates.
    pub junk: bool,
    /// Percent of fragments repeated exactly (PCR duplicates).
    pub repeats: u64,
}

impl Default for Spec {
    fn default() -> Spec {
        Spec {
            refs: vec![
                ("chr1".into(), 60_000),
                ("chr2".into(), 40_000),
                ("chrM".into(), 5_000),
            ],
            fragments: 3000,
            paired: true,
            spliced: true,
            barcodes: 40,
            umis: 30,
            junk: true,
            repeats: 10,
        }
    }
}

/// A CIGAR spanning about `span` reference bases with query length `l`.
fn cigar(rng: &mut Rng, l: u32, spliced: bool) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    let mut left = l;
    if rng.chance(15) && left > 20 {
        let s = rng.range(1, 6) as u32;
        out.push((s, 'S'));
        left -= s;
    }
    let mut tail = 0;
    if rng.chance(15) && left > 20 {
        tail = rng.range(1, 6) as u32;
        left -= tail;
    }
    while left > 0 {
        let m = (rng.range(5, 60) as u32).min(left);
        let op = if rng.chance(10) {
            '='
        } else if rng.chance(5) {
            'X'
        } else {
            'M'
        };
        out.push((m, op));
        left -= m;
        if left == 0 {
            break;
        }
        match rng.range(0, 12) {
            0 | 1 => out.push((rng.range(1, 12) as u32, 'D')),
            2 if spliced => out.push((rng.range(30, 3000) as u32, 'N')),
            3 if left > 3 => {
                let i = rng.range(1, 3) as u32;
                out.push((i, 'I'));
                left -= i;
            }
            _ => {}
        }
    }
    if tail > 0 {
        out.push((tail, 'S'));
    }
    if rng.chance(3) {
        out.insert(0, (rng.range(1, 20) as u32, 'H'));
    }
    out
}

fn ref_span(c: &[(u32, char)]) -> u32 {
    c.iter()
        .filter(|(_, o)| matches!(o, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|(l, _)| l)
        .sum()
}

/// Random reads, sorted by coordinate.
pub fn synth(rng: &mut Rng, spec: &Spec) -> Vec<Read> {
    let mut reads = Vec::new();
    let barcodes: Vec<String> = (0..spec.barcodes)
        .map(|i| {
            let s: String = (0..8)
                .map(|k| b"ACGT"[(i >> (2 * k)) & 3] as char)
                .collect();
            format!("{s}-1")
        })
        .collect();
    let umis: Vec<String> = (0..spec.umis)
        .map(|i| {
            (0..10)
                .map(|k| b"ACGT"[(i * 7 + k * 3) % 4] as char)
                .collect()
        })
        .collect();
    let mut serial = 0;
    let mut previous: Option<Vec<Read>> = None;
    for _ in 0..spec.fragments {
        if let Some(prev) = &previous
            && rng.chance(spec.repeats)
        {
            serial += 1;
            for mut r in prev.clone() {
                r.name = format!("{}d{serial}", r.name);
                reads.push(r);
            }
            continue;
        }
        serial += 1;
        let tid = rng.range(0, spec.refs.len() as u64) as usize;
        let len = spec.refs[tid].1;
        let mut tags = Vec::new();
        if rng.chance(92) {
            tags.push(format!("CB:Z:{}", rng.pick(&barcodes)));
        }
        if rng.chance(92) {
            tags.push(format!("UB:Z:{}", rng.pick(&umis)));
        }
        let multi = spec.junk && rng.chance(5);
        tags.push(format!("NH:i:{}", if multi { 2 } else { 1 }));
        tags.push(format!("RG:Z:rg{}", rng.range(0, 3)));
        let mapq = if spec.junk && rng.chance(15) {
            rng.range(0, 30) as u8
        } else {
            rng.range(30, 61) as u8
        };
        let mut junk_flag = 0u16;
        if spec.junk {
            for (bit, p) in [(0x400u16, 5), (0x200, 2), (0x800, 2)] {
                if rng.chance(p) {
                    junk_flag |= bit;
                }
            }
            if multi && rng.chance(50) {
                junk_flag |= 0x100;
            }
        }
        let (l1, l2) = (rng.range(30, 120) as u32, rng.range(30, 120) as u32);
        let (mut c1, mut c2) = (cigar(rng, l1, spec.spliced), cigar(rng, l2, spec.spliced));
        // Long introns may not fit a short reference.
        if ref_span(&c1) + ref_span(&c2) + 1000 > len {
            c1 = cigar(rng, l1, false);
            c2 = cigar(rng, l2, false);
        }
        let (r1, r2) = (ref_span(&c1), ref_span(&c2));
        let s = rng.range(0, (len - r1 - r2 - 700) as u64) as u32;
        let name = format!("q{serial}");
        if !spec.paired {
            let flag = junk_flag | if rng.chance(50) { 0x10 } else { 0 };
            reads.push(Read {
                name,
                flag,
                tid,
                pos: s,
                mapq,
                cigar: c1,
                mtid: -1,
                mpos: 0,
                tlen: 0,
                tags,
                seq: None,
                qual: None,
            });
            previous = Some(vec![reads.last().unwrap().clone()]);
            continue;
        }
        let frag = rng.range(r1.max(r2) as u64, r1.max(r2) as u64 + 600) as u32;
        let p2 = s + frag - r2;
        let left_is_r1 = rng.chance(50);
        let (fl, fr) = if left_is_r1 {
            (0x40, 0x80)
        } else {
            (0x80, 0x40)
        };
        let end = (s + r1).max(p2 + r2);
        let tl = (end - s) as i64;
        let mut a = Read {
            name: name.clone(),
            flag: 0x1 | 0x2 | 0x20 | fl | junk_flag,
            tid,
            pos: s,
            mapq,
            cigar: c1,
            mtid: tid as i32,
            mpos: p2,
            tlen: tl,
            tags: tags.clone(),
            seq: None,
            qual: None,
        };
        let mut b = Read {
            name,
            flag: 0x1 | 0x2 | 0x10 | fr | junk_flag,
            tid,
            pos: p2,
            mapq: if spec.junk && rng.chance(5) { 10 } else { mapq },
            cigar: c2,
            mtid: tid as i32,
            mpos: s,
            tlen: -tl,
            tags,
            seq: None,
            qual: None,
        };
        if spec.junk && rng.chance(3) {
            // One mate unmapped, placed at the other's position.
            let (m, u) = if rng.chance(50) {
                (&mut a, &mut b)
            } else {
                (&mut b, &mut a)
            };
            m.flag = (m.flag | 0x8) & !0x2;
            u.flag = (u.flag | 0x4) & !0x2 & !0x10;
            u.flag &= !0x20;
            u.pos = m.pos;
            u.mpos = m.pos;
            u.cigar.clear();
            u.mapq = 0;
            u.tlen = 0;
            m.tlen = 0;
            m.mpos = m.pos;
        }
        reads.push(a.clone());
        reads.push(b.clone());
        previous = Some(vec![a, b]);
    }
    reads.sort_by_key(|r| (r.tid, r.pos));
    reads
}

fn cigar_text(c: &[(u32, char)]) -> String {
    if c.is_empty() {
        return "*".into();
    }
    c.iter().map(|(l, o)| format!("{l}{o}")).collect()
}

/// SAM text of `reads` with a coordinate-sorted header.
pub fn sam(refs: &[(String, u32)], reads: &[Read], rng: &mut Rng) -> String {
    let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, l) in refs {
        s.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    for g in 0..3 {
        s.push_str(&format!("@RG\tID:rg{g}\tSM:s{g}\n"));
    }
    for r in reads {
        let ql = if r.cigar.is_empty() {
            40
        } else {
            r.query_len()
        } as usize;
        let mut seq: String = (0..ql).map(|_| *rng.pick(&['A', 'C', 'G', 'T'])).collect();
        let mut qual = if rng.chance(10) {
            "*".to_string()
        } else {
            (0..ql)
                .map(|_| (b'#' + rng.range(0, 38) as u8) as char)
                .collect()
        };
        if let Some(s) = &r.seq {
            seq = String::from_utf8(s.clone()).unwrap();
            qual = match &r.qual {
                Some(q) => String::from_utf8(q.clone()).unwrap(),
                None => "*".into(),
            };
        }
        let rnext = if r.mtid < 0 {
            "*"
        } else if r.mtid == r.tid as i32 {
            "="
        } else {
            &refs[r.mtid as usize].0
        };
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.name,
            r.flag,
            refs[r.tid].0,
            r.pos + 1,
            r.mapq,
            cigar_text(&r.cigar),
            rnext,
            if r.mtid < 0 { 0 } else { r.mpos + 1 },
            r.tlen,
            seq,
            qual
        ));
        for t in &r.tags {
            s.push('\t');
            s.push_str(t);
        }
        s.push('\n');
    }
    s
}

pub fn cstr(p: &Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

/// The BAM htslib writes for the SAM at `src`.
pub fn bam_of(src: &Path, out: &Path) {
    unsafe {
        h::hts_set_log_level(h::HTS_LOG_OFF);
        let r = h::hts_open(cstr(src).as_ptr(), c"r".as_ptr());
        assert!(!r.is_null());
        let hdr = h::sam_hdr_read(r);
        assert!(!hdr.is_null());
        let w = h::hts_open(cstr(out).as_ptr(), c"wb".as_ptr());
        assert_eq!(h::sam_hdr_write(w, hdr), 0);
        let rec = h::bam_init1();
        while h::sam_read1(r, hdr, rec) >= 0 {
            assert!(h::sam_write1(w, hdr, rec) >= 0);
        }
        h::bam_destroy1(rec);
        h::sam_hdr_destroy(hdr);
        h::hts_close(r);
        assert_eq!(h::hts_close(w), 0);
    }
}

/// Writes `reads` as `dir/name.bam` (through SAM); returns its path.
pub fn write_bam(
    dir: &Path,
    name: &str,
    refs: &[(String, u32)],
    reads: &[Read],
    seed: u64,
) -> PathBuf {
    let samp = dir.join(format!("{name}.sam"));
    std::fs::write(&samp, sam(refs, reads, &mut Rng::new(seed))).unwrap();
    let bam = dir.join(format!("{name}.bam"));
    bam_of(&samp, &bam);
    bam
}

fn args(v: &[String]) -> (Vec<CString>, Vec<*mut c_char>) {
    let owned: Vec<CString> = v
        .iter()
        .map(|a| CString::new(a.as_str()).unwrap())
        .collect();
    let mut ptrs: Vec<*mut c_char> = owned.iter().map(|a| a.as_ptr() as *mut c_char).collect();
    ptrs.push(std::ptr::null_mut());
    (owned, ptrs)
}

unsafe extern "C" {
    static mut optind: c_int;
}

/// `samtools <args>`.
pub fn samtools(a: &[String]) -> i32 {
    let _g = TOOLS.lock().unwrap_or_else(|e| e.into_inner());
    let (_o, mut p) = args(a);
    unsafe { samtools_sys::ow_samtools(a.len() as c_int, p.as_mut_ptr()) }
}

/// `featureCounts <args>`.
pub fn feature_counts(a: &[String]) -> i32 {
    let _g = TOOLS.lock().unwrap_or_else(|e| e.into_inner());
    let mut full = vec!["featureCounts".to_string()];
    full.extend_from_slice(a);
    let (_o, mut p) = args(&full);
    unsafe {
        optind = 0;
        subread_sys::ow_featurecounts_main(full.len() as c_int, p.as_mut_ptr())
    }
}

/// `samtools depth -aa <extra> bam`: per reference, the depth of every base.
pub fn samtools_depth(bam: &Path, extra: &[&str], refs: &[(String, u32)]) -> Vec<Vec<u32>> {
    let out = bam.with_extension("depth.txt");
    let mut a: Vec<String> = vec!["depth".into(), "-aa".into(), "-o".into()];
    a.push(out.to_str().unwrap().into());
    a.extend(extra.iter().map(|s| s.to_string()));
    a.push(bam.to_str().unwrap().into());
    assert_eq!(samtools(&a), 0, "samtools depth failed");
    let mut depth: Vec<Vec<u32>> = refs.iter().map(|r| vec![0; r.1 as usize]).collect();
    for line in std::fs::read_to_string(&out).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let t = refs.iter().position(|r| r.0 == f[0]).unwrap();
        let p: usize = f[1].parse().unwrap();
        depth[t][p - 1] = f[2].parse().unwrap();
    }
    depth
}

/// Every base of a bigWig: NaN (no data) as 0.
pub fn bigwig_values(path: &Path, refs: &[(String, u32)]) -> Vec<Vec<f32>> {
    let mut r = omicsweave_bbi::BbiReader::open(FileSource::open(path).unwrap()).unwrap();
    refs.iter()
        .map(|(n, l)| {
            let mut v = vec![0f32; *l as usize];
            r.values(
                &omicsweave_core::Region {
                    chrom: n,
                    start: 0,
                    end: *l,
                }
                .into(),
                &mut v,
            )
            .unwrap();
            v.iter()
                .map(|x| if x.is_nan() { 0.0 } else { *x })
                .collect()
        })
        .collect()
}

/// The text of a gzip / BGZF file.
pub fn gunzip(path: &Path) -> String {
    let raw = std::fs::read(path).unwrap();
    let mut out = String::new();
    flate2::read::MultiGzDecoder::new(&raw[..])
        .read_to_string(&mut out)
        .unwrap();
    out
}

/// A 10x directory's matrix as `(row key, feature id) → count`.
pub fn read_10x(dir: &Path) -> BTreeMap<(String, String), u32> {
    read_10x_layer(dir, "matrix")
}

/// One layer of a 10x directory (`matrix`, `spliced`, `labeled`, ...).
pub fn read_10x_layer(dir: &Path, layer: &str) -> BTreeMap<(String, String), u32> {
    let rows: Vec<String> = gunzip(&dir.join("barcodes.tsv.gz"))
        .lines()
        .map(str::to_owned)
        .collect();
    let feats: Vec<String> = gunzip(&dir.join("features.tsv.gz"))
        .lines()
        .map(|l| l.split('\t').next().unwrap().to_owned())
        .collect();
    let mut out = BTreeMap::new();
    let file = |ext: &str| dir.join(format!("{layer}.{ext}"));
    if file("cpk").exists() {
        let f = cellpack::CpkFile::open(file("cpk")).unwrap();
        assert_eq!(f.n_rows() as usize, rows.len());
        assert_eq!(f.n_cols() as usize, feats.len());
        let csr = f.read_csr(0..f.n_rows()).unwrap();
        for (r, row) in rows.iter().enumerate() {
            for k in csr.indptr[r] as usize..csr.indptr[r + 1] as usize {
                let v = u32::from_le_bytes(csr.data[4 * k..4 * k + 4].try_into().unwrap());
                out.insert((row.clone(), feats[csr.indices[k] as usize].clone()), v);
            }
        }
    } else if file("mtx.gz").exists() {
        let text = gunzip(&file("mtx.gz"));
        let mut lines = text.lines().filter(|l| !l.starts_with('%'));
        let dims: Vec<usize> = lines
            .next()
            .unwrap()
            .split(' ')
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!((dims[0], dims[1]), (feats.len(), rows.len()));
        let mut n = 0;
        for l in lines {
            let v: Vec<usize> = l.split(' ').map(|x| x.parse().unwrap()).collect();
            out.insert(
                (rows[v[1] - 1].clone(), feats[v[0] - 1].clone()),
                v[2] as u32,
            );
            n += 1;
        }
        assert_eq!(n, dims[2]);
    } else {
        let text = gunzip(&file("csv.gz"));
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split(',').skip(1).collect();
        assert_eq!(head.len(), rows.len());
        for l in lines {
            let f: Vec<&str> = l.split(',').collect();
            for (i, v) in f[1..].iter().enumerate() {
                let v: u32 = v.parse().unwrap();
                if v > 0 {
                    out.insert((rows[i].clone(), f[0].to_owned()), v);
                }
            }
        }
    }
    out
}

/// A gene of the generated annotation.
#[derive(Debug, Clone)]
pub struct Gene {
    pub id: String,
    pub tid: usize,
    pub start: u32,
    pub end: u32,
    pub strand: char,
    /// Exons of each transcript, 0-based half-open.
    pub transcripts: Vec<Vec<(u32, u32)>>,
}

impl Gene {
    /// The union of the exons.
    pub fn exonic(&self) -> Vec<(u32, u32)> {
        let mut all: Vec<(u32, u32)> = self.transcripts.iter().flatten().copied().collect();
        all.sort();
        let mut out: Vec<(u32, u32)> = Vec::new();
        for (s, e) in all {
            match out.last_mut() {
                Some(l) if s <= l.1 => l.1 = l.1.max(e),
                _ => out.push((s, e)),
            }
        }
        out
    }
}

/// Random genes over `refs`, some overlapping, on both strands.
pub fn genes(rng: &mut Rng, refs: &[(String, u32)], per_ref: usize) -> Vec<Gene> {
    let mut out = Vec::new();
    let mut n = 0;
    for (tid, (_, len)) in refs.iter().enumerate() {
        for _ in 0..per_ref {
            n += 1;
            let span = rng.range(800, 8000.min(*len as u64 / 2)) as u32;
            let start = rng.range(0, (*len - span) as u64) as u32;
            let strand = if rng.chance(50) { '+' } else { '-' };
            let mut transcripts = Vec::new();
            for _ in 0..rng.range(1, 4) {
                let k = rng.range(1, 5) as u32;
                let piece = span / k;
                let mut ex = Vec::new();
                for i in 0..k {
                    let s = start + i * piece + rng.range(0, (piece / 3) as u64) as u32;
                    let e =
                        (s + rng.range(50, (piece / 2).max(51) as u64) as u32).min(start + span);
                    if e > s {
                        ex.push((s, e));
                    }
                }
                transcripts.push(ex);
            }
            let s = transcripts.iter().flatten().map(|x| x.0).min().unwrap();
            let e = transcripts.iter().flatten().map(|x| x.1).max().unwrap();
            out.push(Gene {
                id: format!("g{n}"),
                tid,
                start: s,
                end: e,
                strand,
                transcripts,
            });
        }
    }
    out
}

/// The genes as GTF.
pub fn gtf(refs: &[(String, u32)], genes: &[Gene]) -> String {
    let mut s = String::new();
    let mut t = 0;
    for g in genes {
        let c = &refs[g.tid].0;
        let attrs = format!(
            "gene_id \"{0}\"; gene_name \"N{0}\"; gene_type \"protein_coding\";",
            g.id
        );
        s.push_str(&format!(
            "{c}\tsyn\tgene\t{}\t{}\t.\t{}\t.\t{attrs}\n",
            g.start + 1,
            g.end,
            g.strand
        ));
        for ex in &g.transcripts {
            t += 1;
            let ts = ex.iter().map(|x| x.0).min().unwrap();
            let te = ex.iter().map(|x| x.1).max().unwrap();
            s.push_str(&format!(
                "{c}\tsyn\ttranscript\t{}\t{te}\t.\t{}\t.\t{attrs} transcript_id \"t{t}\";\n",
                ts + 1,
                g.strand
            ));
            for &(a, b) in ex {
                s.push_str(&format!(
                    "{c}\tsyn\texon\t{}\t{b}\t.\t{}\t.\t{attrs} transcript_id \"t{t}\";\n",
                    a + 1,
                    g.strand
                ));
            }
        }
    }
    s
}

/// Whether `[s, e)` overlaps any of `ivs`.
pub fn overlaps(ivs: &[(u32, u32)], s: u32, e: u32) -> bool {
    ivs.iter().any(|&(a, b)| a < e && s < b)
}

/// The genes (indices) whose exons the blocks overlap, on `strand` when
/// given (`Some(true)`: reverse).
pub fn gene_hits(
    genes: &[Gene],
    tid: usize,
    blocks: &[(u32, u32)],
    strand: Option<bool>,
) -> Vec<usize> {
    genes
        .iter()
        .enumerate()
        .filter(|(_, g)| g.tid == tid)
        .filter(|(_, g)| strand.is_none_or(|rev| (g.strand == '-') == rev))
        .filter(|(_, g)| {
            let ex = g.exonic();
            blocks.iter().any(|&(s, e)| overlaps(&ex, s, e))
        })
        .map(|(i, _)| i)
        .collect()
}

/// The fragment strand of a read under a library type (0 / 1 / 2).
pub fn reverse_for(r: &Read, lib: u8) -> bool {
    let read2 = r.flag & 0x81 == 0x81;
    (r.reverse() ^ read2) ^ (lib == 2)
}

/// A bigBed of `peaks` (sorted), with the chromosome sizes of `refs`.
pub fn big_bed(path: &Path, refs: &[(String, u32)], peaks: &[(usize, u32, u32)]) {
    let chroms = refs
        .iter()
        .map(|(n, l)| omicsweave_bbi::Chrom {
            name: n.clone(),
            len: *l,
        })
        .collect();
    let file = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    let mut w = omicsweave_bbi::BigBedWriter::create(
        file,
        chroms,
        omicsweave_bbi::BedSchema::bed(4),
        omicsweave_bbi::WriteOptions::default(),
    )
    .unwrap();
    for (i, &(t, s, e)) in peaks.iter().enumerate() {
        w.add(t as u32, s, e, format!("p{i}").as_bytes()).unwrap();
    }
    w.finish().unwrap();
}

/// Random bases for every reference.
pub fn genome(rng: &mut Rng, refs: &[(String, u32)]) -> Vec<Vec<u8>> {
    refs.iter()
        .map(|(_, len)| (0..*len).map(|_| *rng.pick(b"ACGT")).collect())
        .collect()
}

/// `seqs` as bgzip FASTA at `dir/name.fa.gz`, with the `.gzi` bgzip
/// writes and the `.fai` of htslib's `fai_build3`.
pub fn write_fasta(dir: &Path, name: &str, refs: &[(String, u32)], seqs: &[Vec<u8>]) -> PathBuf {
    let mut text = Vec::new();
    for ((n, _), s) in refs.iter().zip(seqs) {
        text.extend_from_slice(format!(">{n}\n").as_bytes());
        for line in s.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    let path = dir.join(format!("{name}.fa.gz"));
    let p = cstr(&path);
    unsafe {
        let fp = h::bgzf_open(p.as_ptr(), c"w".as_ptr());
        assert!(!fp.is_null());
        assert_eq!(h::bgzf_index_build_init(fp), 0);
        assert_eq!(
            h::bgzf_write(fp, text.as_ptr().cast(), text.len()),
            text.len() as isize
        );
        assert_eq!(h::bgzf_index_dump(fp, p.as_ptr(), c".gzi".as_ptr()), 0);
        assert_eq!(h::bgzf_close(fp), 0);
        assert_eq!(
            h::fai_build3(p.as_ptr(), std::ptr::null(), std::ptr::null()),
            0
        );
    }
    path
}

/// `seqs` as 2bit at `dir/name.2bit`.
pub fn write_twobit(dir: &Path, name: &str, refs: &[(String, u32)], seqs: &[Vec<u8>]) -> PathBuf {
    use omicsweave_fasta::{TwoBitWriter, encode_twobit};
    let mut w = TwoBitWriter::create(
        std::io::Cursor::new(Vec::new()),
        refs.iter().map(|r| r.0.clone()).collect(),
        0,
    )
    .unwrap();
    for ((n, _), s) in refs.iter().zip(seqs) {
        let (record, packed) = encode_twobit(s).unwrap();
        w.add(n, &record, &mut packed.as_slice()).unwrap();
    }
    let path = dir.join(format!("{name}.2bit"));
    std::fs::write(&path, w.finish().unwrap().into_inner()).unwrap();
    path
}

/// Sites `(tid, pos, ref, alt)` (sorted) as a BCF with its CSI, written by
/// htslib.
pub fn write_bcf(
    dir: &Path,
    name: &str,
    refs: &[(String, u32)],
    sites: &[(usize, u32, u8, u8)],
) -> PathBuf {
    let mut vcf = String::from("##fileformat=VCFv4.2\n");
    for (n, l) in refs {
        vcf.push_str(&format!("##contig=<ID={n},length={l}>\n"));
    }
    vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for &(t, p, r, a) in sites {
        vcf.push_str(&format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t.\n",
            refs[t].0,
            p + 1,
            r as char,
            a as char
        ));
    }
    let src = dir.join(format!("{name}.vcf"));
    std::fs::write(&src, vcf).unwrap();
    let out = dir.join(format!("{name}.bcf"));
    unsafe {
        let r = h::hts_open(cstr(&src).as_ptr(), c"r".as_ptr());
        assert!(!r.is_null());
        let hdr = h::bcf_hdr_read(r);
        let w = h::hts_open(cstr(&out).as_ptr(), c"wb".as_ptr());
        assert_eq!(h::bcf_hdr_write(w, hdr), 0);
        let rec = h::bcf_init();
        while h::bcf_read(r, hdr, rec) == 0 {
            assert_eq!(h::bcf_write(w, hdr, rec), 0);
        }
        h::bcf_destroy(rec);
        h::bcf_hdr_destroy(hdr);
        h::hts_close(r);
        assert_eq!(h::hts_close(w), 0);
        assert_eq!(
            h::bcf_index_build3(cstr(&out).as_ptr(), std::ptr::null(), 14, 1),
            0
        );
    }
    out
}

/// The complement of a base.
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

/// `(ref_pos, query_pos)` of every `M` / `=` / `X` base.
pub fn aligned_pairs(r: &Read) -> Vec<(u32, usize)> {
    let mut out = Vec::new();
    let (mut at, mut q) = (r.pos, 0usize);
    for &(l, o) in &r.cigar {
        match o {
            'M' | '=' | 'X' => {
                for i in 0..l {
                    out.push((at + i, q + i as usize));
                }
                at += l;
                q += l as usize;
            }
            'I' | 'S' => q += l as usize,
            'D' | 'N' => at += l,
            _ => {}
        }
    }
    out
}

/// The `MD` tag of `r` against `seq` (the read's reference).
pub fn md_tag(r: &Read, genome: &[u8]) -> String {
    let bases = r.seq.as_ref().unwrap();
    let mut md = String::new();
    let mut run = 0;
    let mut at = r.pos;
    let mut q = 0usize;
    for &(l, o) in &r.cigar {
        match o {
            'M' | '=' | 'X' => {
                for i in 0..l as usize {
                    let g = genome[(at as usize) + i];
                    if bases[q + i] == g {
                        run += 1;
                    } else {
                        md.push_str(&format!("{run}{}", g as char));
                        run = 0;
                    }
                }
                at += l;
                q += l as usize;
            }
            'D' => {
                md.push_str(&format!("{run}^"));
                for i in 0..l as usize {
                    md.push(genome[at as usize + i] as char);
                }
                run = 0;
                at += l;
            }
            'N' => at += l,
            'I' | 'S' => q += l as usize,
            _ => {}
        }
    }
    md.push_str(&run.to_string());
    md
}

/// htslib's `tbx_conf_t`.
#[repr(C)]
pub struct TbxConf {
    preset: i32,
    sc: i32,
    bc: i32,
    ec: i32,
    meta_char: i32,
    line_skip: i32,
}

unsafe extern "C" {
    static tbx_conf_bed: TbxConf;
    fn tbx_index_build(fname: *const c_char, min_shift: c_int, conf: *const TbxConf) -> c_int;
}

/// A tabix index as (header with names, per reference: bins with their
/// chunks and the linear index, unplaced count); bins in any order.
#[allow(clippy::type_complexity)]
pub fn parse_tbi(
    b: &[u8],
) -> (
    Vec<u8>,
    Vec<(BTreeMap<u32, Vec<(u64, u64)>>, Vec<u64>)>,
    Vec<u8>,
) {
    let u32_at = |i: usize| u32::from_le_bytes(b[i..i + 4].try_into().unwrap());
    let u64_at = |i: usize| u64::from_le_bytes(b[i..i + 8].try_into().unwrap());
    assert_eq!(&b[..4], b"TBI\x01");
    let n_ref = u32_at(4) as usize;
    let l_nm = u32_at(32) as usize;
    let mut at = 36 + l_nm;
    let head = b[..at].to_vec();
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let mut bins = BTreeMap::new();
        let n_bin = u32_at(at);
        at += 4;
        for _ in 0..n_bin {
            let bin = u32_at(at);
            let n = u32_at(at + 4) as usize;
            at += 8;
            let chunks = (0..n)
                .map(|k| (u64_at(at + 16 * k), u64_at(at + 16 * k + 8)))
                .collect();
            at += 16 * n;
            bins.insert(bin, chunks);
        }
        let n_intv = u32_at(at) as usize;
        at += 4;
        let linear = (0..n_intv).map(|k| u64_at(at + 8 * k)).collect();
        at += 8 * n_intv;
        refs.push((bins, linear));
    }
    (head, refs, b[at..].to_vec())
}

/// The bytes of a gzip / BGZF file, inflated.
pub fn gunzip_bytes(p: &Path) -> Vec<u8> {
    let raw = std::fs::read(p).unwrap();
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(&raw[..])
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Checks that `<gz>.tbi` is the index htslib's `tbx_index_build` (BED
/// preset) writes for a copy of `gz` in `scratch` (bins in any order).
pub fn assert_tabix_as_htslib(gz: &Path, scratch: &Path) {
    let copy = scratch.join(format!(
        "copy-{}",
        gz.file_name().unwrap().to_str().unwrap()
    ));
    std::fs::copy(gz, &copy).unwrap();
    {
        let _g = TOOLS.lock().unwrap_or_else(|e| e.into_inner());
        assert_eq!(
            unsafe { tbx_index_build(cstr(&copy).as_ptr(), 0, &tbx_conf_bed) },
            0
        );
    }
    let mut ours = gz.as_os_str().to_owned();
    ours.push(".tbi");
    let mut theirs = copy.as_os_str().to_owned();
    theirs.push(".tbi");
    assert_eq!(
        parse_tbi(&gunzip_bytes(Path::new(&ours))),
        parse_tbi(&gunzip_bytes(Path::new(&theirs)))
    );
}

/// VCF text (header lines included) as a BCF with its CSI, written by
/// htslib.
pub fn bcf_of_vcf(dir: &Path, name: &str, vcf: &str) -> PathBuf {
    let src = dir.join(format!("{name}.vcf"));
    std::fs::write(&src, vcf).unwrap();
    let out = dir.join(format!("{name}.bcf"));
    unsafe {
        let r = h::hts_open(cstr(&src).as_ptr(), c"r".as_ptr());
        assert!(!r.is_null());
        let hdr = h::bcf_hdr_read(r);
        let w = h::hts_open(cstr(&out).as_ptr(), c"wb".as_ptr());
        assert_eq!(h::bcf_hdr_write(w, hdr), 0);
        let rec = h::bcf_init();
        while h::bcf_read(r, hdr, rec) == 0 {
            assert_eq!(h::bcf_write(w, hdr, rec), 0);
        }
        h::bcf_destroy(rec);
        h::bcf_hdr_destroy(hdr);
        h::hts_close(r);
        assert_eq!(h::hts_close(w), 0);
        assert_eq!(
            h::bcf_index_build3(cstr(&out).as_ptr(), std::ptr::null(), 14, 1),
            0
        );
    }
    out
}

/// Every file under `dir` with its bytes, by relative path.
pub fn files(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_owned()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                stack.push(p);
            } else {
                out.push((
                    p.strip_prefix(dir).unwrap().display().to_string(),
                    std::fs::read(&p).unwrap(),
                ));
            }
        }
    }
    out.sort();
    out
}
