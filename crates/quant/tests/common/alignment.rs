use super::support::{Rng, cstr};
use htslib_sys as h;
use std::path::{Path, PathBuf};

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

/// The BAM htslib writes for the SAM at `src`.
pub fn bam_of(src: &Path, out: &Path) {
    unsafe {
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
