use super::alignment::Read;
use super::support::{Rng, cstr};
use htslib_sys as h;
use std::path::{Path, PathBuf};

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

/// A CIGAR with query length `l`, including indels, clips and optional introns.
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
    bcf_of_vcf(dir, name, &vcf)
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
