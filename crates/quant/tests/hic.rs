//! Hi-C contacts against direct computation on generated read pairs (cis
//! and trans, both strands, clips and indels, `MC` / `MQ` tags sometimes
//! missing, low MAPQ on either mate, flagged records, unflagged duplicates,
//! barcodes): the 4DN pairs file, per-barcode statistics, and the `.hic`
//! file, byte for byte the one `omicsweave_hic::HicBuilder` writes for the
//! same contacts and read back through `HicReader`; the same files at any
//! thread count and with spilled buffers.

mod common;

use std::collections::BTreeMap;
use std::io::Cursor;
use std::path::Path;

use common::{Read, Rng, files, gunzip, write_bam};
use omicsweave_core::{ContactQuery, ContactRead, FileSource, Region};
use omicsweave_hic::{BuildOptions, Chrom, HicBuilder, HicHeader, HicReader};
use omicsweave_quant::{Budget, Contacts, Product, Task, plan, quantify, run};

const REFS: [(&str, u32); 3] = [("chr1", 900_000), ("chr2", 600_000), ("chr3", 250_000)];
const RESOLUTIONS: [u32; 3] = [100_000, 25_000, 5_000];

struct Fixture {
    dir: tempfile::TempDir,
    refs: Vec<(String, u32)>,
    reads: Vec<Read>,
    bam: std::path::PathBuf,
}

fn cigar(rng: &mut Rng) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    if rng.chance(20) {
        out.push((rng.range(1, 20) as u32, 'S'));
    }
    out.push((rng.range(20, 60) as u32, 'M'));
    match rng.range(0, 6) {
        0 => out.push((rng.range(1, 4) as u32, 'D')),
        1 => out.push((rng.range(1, 4) as u32, 'I')),
        _ => {}
    }
    out.push((rng.range(20, 60) as u32, 'M'));
    if rng.chance(20) {
        out.push((rng.range(1, 20) as u32, 'S'));
    }
    out
}

fn cigar_text(c: &[(u32, char)]) -> String {
    c.iter().map(|(l, o)| format!("{l}{o}")).collect()
}

/// One end of a pair before it is written: reference, position, CIGAR,
/// strand.
#[derive(Clone)]
struct End {
    tid: usize,
    pos: u32,
    cigar: Vec<(u32, char)>,
    reverse: bool,
}

fn fixture(seed: u64) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let refs: Vec<(String, u32)> = REFS.iter().map(|(n, l)| (n.to_string(), *l)).collect();
    let mut rng = Rng::new(seed);
    let end = |rng: &mut Rng, tid: usize, near: Option<u32>| -> End {
        let c = cigar(rng);
        let span = c
            .iter()
            .filter(|x| matches!(x.1, 'M' | 'D'))
            .map(|x| x.0)
            .sum::<u32>();
        let len = refs[tid].1;
        let pos = match near {
            Some(p) => {
                let d = match rng.range(0, 3) {
                    0 => rng.range(0, 1_500),
                    1 => rng.range(0, 30_000),
                    _ => rng.range(0, 200_000),
                } as i64;
                let q = if rng.chance(50) {
                    p as i64 + d
                } else {
                    p as i64 - d
                };
                q.clamp(0, (len - span - 1) as i64) as u32
            }
            None => rng.range(0, (len - span - 1) as u64) as u32,
        };
        End {
            tid,
            pos,
            cigar: c,
            reverse: rng.chance(50),
        }
    };
    let mut pairs: Vec<(End, End, Option<String>)> = Vec::new();
    for _ in 0..4000 {
        let t = rng.range(0, 3) as usize;
        let a = end(&mut rng, t, None);
        let b = if rng.chance(70) {
            end(&mut rng, t, Some(a.pos))
        } else {
            let u = (t + rng.range(1, 3) as usize) % 3;
            end(&mut rng, u, None)
        };
        let cb = rng.chance(95).then(|| format!("C{}", rng.range(0, 5)));
        pairs.push((a.clone(), b.clone(), cb.clone()));
        // PCR duplicates: the same ends again, some with another barcode.
        if rng.chance(6) {
            let cb2 = if rng.chance(70) {
                cb
            } else {
                Some(format!("C{}", rng.range(0, 5)))
            };
            pairs.push((a, b, cb2));
        }
    }
    let mut reads = Vec::new();
    for (i, (a, b, cb)) in pairs.into_iter().enumerate() {
        let name = format!("p{i}");
        let dup = rng.chance(2);
        let mate_unmapped = rng.chance(1);
        let mapq = [
            if rng.chance(8) { 12 } else { 60 },
            if rng.chance(8) { 12 } else { 60 },
        ];
        for (k, (me, mate)) in [(&a, &b), (&b, &a)].into_iter().enumerate() {
            let mut flag = 0x1 | if k == 0 { 0x40 } else { 0x80 };
            if me.reverse {
                flag |= 0x10;
            }
            if mate.reverse {
                flag |= 0x20;
            }
            if dup {
                flag |= 0x400;
            }
            if mate_unmapped {
                flag |= 0x8;
            }
            match rng.range(0, 100) {
                0 => flag |= 0x100,
                1 => flag |= 0x800,
                2 => flag |= 0x200,
                _ => {}
            }
            let mut tags = Vec::new();
            if !rng.chance(4) {
                tags.push(format!("MC:Z:{}", cigar_text(&mate.cigar)));
            }
            if !rng.chance(3) {
                tags.push(format!("MQ:i:{}", mapq[1 - k]));
            }
            if let Some(c) = &cb {
                tags.push(format!("CB:Z:{c}"));
            }
            reads.push(Read {
                name: name.clone(),
                flag,
                tid: me.tid,
                pos: me.pos,
                mapq: mapq[k],
                cigar: me.cigar.clone(),
                mtid: mate.tid as i32,
                mpos: mate.pos,
                tlen: 0,
                tags,
                seq: None,
                qual: None,
            });
        }
    }
    reads.sort_by_key(|r| (r.tid, r.pos));
    let bam = write_bam(dir.path(), "hic", &refs, &reads, seed);
    Fixture {
        dir,
        refs,
        reads,
        bam,
    }
}

/// A contact as computed here: ends ordered, MAPQs, name, barcode, and the
/// leader's number in the BAM.
#[derive(Debug, Clone)]
struct Want {
    ends: (usize, u32, bool, u8, usize, u32, bool, u8),
    name: String,
    barcode: Option<String>,
    ordinal: usize,
}

fn ref_len(c: &[(u32, char)]) -> u32 {
    c.iter()
        .filter(|x| matches!(x.1, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|x| x.0)
        .sum()
}

fn parse_cigar(t: &str) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    let mut n = 0;
    for c in t.chars() {
        if let Some(d) = c.to_digit(10) {
            n = n * 10 + d;
        } else {
            out.push((n, c));
            n = 0;
        }
    }
    out
}

/// The contacts the definitions give, in the pairs file's order, after
/// duplicates are removed, with the barcode rows in order of first use.
fn contacts(f: &Fixture, c: &Contacts) -> (Vec<Want>, Vec<String>, u64) {
    let mut out = Vec::new();
    let mut rows: Vec<String> = Vec::new();
    for (ordinal, r) in f.reads.iter().enumerate() {
        let fl = r.flag;
        if fl & c.filter.exclude != 0 || fl & c.filter.require != c.filter.require {
            continue;
        }
        if r.mapq < c.filter.min_mapq {
            continue;
        }
        let me = (r.tid, r.pos);
        let mate = (r.mtid as usize, r.mpos);
        if me > mate || (me == mate && fl & 0x40 == 0) {
            continue;
        }
        let Some(mq) = r.tag("MQ").map(|v| v.parse::<u8>().unwrap()) else {
            continue;
        };
        if mq < c.filter.min_mapq {
            continue;
        }
        let mate_rev = fl & 0x20 != 0;
        let mate_five = if mate_rev {
            match r.tag("MC") {
                Some(t) => r.mpos + ref_len(&parse_cigar(t)),
                None => continue,
            }
        } else {
            r.mpos + 1
        };
        let own_five = if r.reverse() {
            r.pos + ref_len(&r.cigar)
        } else {
            r.pos + 1
        };
        let a = (r.tid, own_five, r.reverse(), r.mapq);
        let b = (r.mtid as usize, mate_five, mate_rev, mq);
        let (a, b) = if (a.0, a.1) <= (b.0, b.1) {
            (a, b)
        } else {
            (b, a)
        };
        if a.0 == b.0 {
            let d = b.1 - a.1;
            if d < c.min_distance || (a.2 != b.2 && d < c.self_ligation) {
                continue;
            }
        }
        let barcode = match c.barcode_tag {
            Some(_) => match r.tag("CB") {
                Some(v) => {
                    if !rows.iter().any(|x| x == v) {
                        rows.push(v.to_owned());
                    }
                    Some(v.to_owned())
                }
                None => continue,
            },
            None => None,
        };
        out.push(Want {
            ends: (a.0, a.1, a.2, a.3, b.0, b.1, b.2, b.3),
            name: r.name.clone(),
            barcode,
            ordinal,
        });
    }
    let row = |w: &Want| {
        w.barcode.as_ref().map_or(u32::MAX as usize, |b| {
            rows.iter().position(|x| x == b).unwrap()
        })
    };
    let key = |w: &Want| {
        let e = w.ends;
        (e.0, e.4, e.1, e.5, e.2, e.6, row(w))
    };
    out.sort_by_key(|w| (key(w), w.ordinal));
    let mut kept: Vec<Want> = Vec::new();
    let mut dups = 0;
    for w in out {
        if c.dedup && kept.last().is_some_and(|k| key(k) == key(&w)) {
            dups += 1;
            continue;
        }
        kept.push(w);
    }
    (kept, rows, dups)
}

fn pairs_text(f: &Fixture, want: &[Want], genome: &str, barcode: bool) -> String {
    let mut s = String::from(
        "## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n",
    );
    s.push_str(&format!("#genome_assembly: {genome}\n"));
    for (n, l) in &f.refs {
        s.push_str(&format!("#chromsize: {n} {l}\n"));
    }
    s.push_str("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type mapq1 mapq2");
    s.push_str(if barcode { " barcode\n" } else { "\n" });
    let st = |r: bool| if r { '-' } else { '+' };
    for w in want {
        let e = w.ends;
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\tUU\t{}\t{}",
            w.name,
            f.refs[e.0].0,
            e.1,
            f.refs[e.4].0,
            e.5,
            st(e.2),
            st(e.6),
            e.3,
            e.7
        ));
        if let Some(b) = &w.barcode {
            s.push('\t');
            s.push_str(b);
        }
        s.push('\n');
    }
    s
}

fn hic_bytes(f: &Fixture, want: &[Want], genome: &str) -> Vec<u8> {
    let header = HicHeader::new(
        genome,
        f.refs
            .iter()
            .map(|(n, l)| Chrom {
                name: n.clone(),
                len: *l as u64,
            })
            .collect(),
        RESOLUTIONS.to_vec(),
    );
    let mut b = HicBuilder::new(header).unwrap();
    for w in want {
        let e = w.ends;
        b.add(e.0 as u32, e.1, e.4 as u32, e.5, 1.0).unwrap();
    }
    b.finish(Cursor::new(Vec::new()), &BuildOptions::default())
        .unwrap()
        .into_inner()
}

fn task(f: &Fixture, out: &Path, c: Contacts) -> Task {
    let mut t = Task::new(vec![f.bam.clone()], out);
    t.products = vec![Product::Contacts(c)];
    t
}

fn sc_contacts() -> Contacts {
    Contacts {
        resolutions: RESOLUTIONS.to_vec(),
        genome: "syn1".into(),
        barcode_tag: Some(*b"CB"),
        self_ligation: 1_000,
        ..Contacts::default()
    }
}

#[test]
fn pairs_hic_and_cell_metrics_match_the_reads() {
    let f = fixture(81);
    let spec = sc_contacts();
    let out = f.dir.path().join("sc");
    quantify(
        task(&f, &out, spec.clone()),
        Budget::new(3, 256 << 20),
        None,
    )
    .unwrap();
    let (want, rows, dups) = contacts(&f, &spec);
    assert!(want.len() > 1500 && dups > 50, "{} {dups}", want.len());
    assert_eq!(
        gunzip(&out.join("hic/contacts.pairs.gz")),
        pairs_text(&f, &want, "syn1", true)
    );
    let hic = out.join("hic/contacts.hic");
    let got = std::fs::read(&hic).unwrap();
    let expected = hic_bytes(&f, &want, "syn1");
    assert!(
        got == expected,
        "{} vs {} bytes, first difference at {:?}",
        got.len(),
        expected.len(),
        got.iter().zip(&expected).position(|(a, b)| a != b)
    );

    // Read back: every cell of each map at 25 kb, both halves of a
    // chromosome's own map.
    let mut r = HicReader::open(FileSource::open(&hic).unwrap()).unwrap();
    let res = 25_000;
    for (i, (ni, li)) in f.refs.iter().enumerate() {
        for (j, (nj, lj)) in f.refs.iter().enumerate().skip(i) {
            let mut cells: BTreeMap<(u32, u32), f32> = BTreeMap::new();
            for w in &want {
                let e = w.ends;
                if (e.0, e.4) != (i, j) {
                    continue;
                }
                let (x, y) = (e.1.min(*li) / res * res, e.5.min(*lj) / res * res);
                *cells.entry((x, y)).or_default() += 1.0;
                if i == j && x != y {
                    *cells.entry((y, x)).or_default() += 1.0;
                }
            }
            let q = ContactQuery::new(Region::new(ni, 0, *li), Region::new(nj, 0, *lj), res);
            let n = r.query_entries(q).unwrap();
            let mut entries = vec![omicsweave_core::ContactEntry::default(); n];
            r.write_entries(&mut entries).unwrap();
            let read: BTreeMap<(u32, u32), f32> =
                entries.iter().map(|e| ((e.row, e.col), e.value)).collect();
            assert_eq!(read, cells, "{ni} x {nj}");
        }
    }

    // Per-barcode statistics.
    let csv = std::fs::read_to_string(out.join("hic/cell_metrics.csv")).unwrap();
    let mut lines = csv.lines();
    assert_eq!(
        lines.next().unwrap(),
        "barcode,pairs,duplicates,contacts,cis,trans,cis_1kb+,cis_10kb+,cis_20kb+"
    );
    let got: Vec<&str> = lines.collect();
    assert_eq!(got.len(), rows.len());
    for (line, bc) in got.iter().zip(&rows) {
        let mine: Vec<&Want> = want
            .iter()
            .filter(|w| w.barcode.as_ref() == Some(bc))
            .collect();
        let cis: Vec<u32> = mine
            .iter()
            .filter(|w| w.ends.0 == w.ends.4)
            .map(|w| w.ends.5 - w.ends.1)
            .collect();
        let fields: Vec<&str> = line.split(',').collect();
        assert_eq!(fields[0], bc);
        let n = |i: usize| fields[i].parse::<usize>().unwrap();
        assert_eq!(n(1), n(2) + n(3));
        assert_eq!(n(3), mine.len());
        assert_eq!(n(4), cis.len());
        assert_eq!(n(5), mine.len() - cis.len());
        assert_eq!(n(6), cis.iter().filter(|&&d| d >= 1_000).count());
        assert_eq!(n(7), cis.iter().filter(|&&d| d >= 10_000).count());
        assert_eq!(n(8), cis.iter().filter(|&&d| d >= 20_000).count());
    }
    let summary = std::fs::read_to_string(out.join("summary.json")).unwrap();
    assert!(summary.contains(&format!("\"contacts\": {},", want.len())));
    assert!(summary.contains(&format!("\"duplicates\": {dups},")));
}

#[test]
fn filters_and_bulk_contacts() {
    let f = fixture(82);
    let spec = Contacts {
        resolutions: RESOLUTIONS.to_vec(),
        min_distance: 2_000,
        dedup: false,
        filter: omicsweave_quant::Filter {
            min_mapq: 0,
            ..Contacts::default().filter
        },
        ..Contacts::default()
    };
    let out = f.dir.path().join("bulk");
    quantify(
        task(&f, &out, spec.clone()),
        Budget::new(2, 256 << 20),
        None,
    )
    .unwrap();
    // Without a MAPQ threshold, a record without MQ is still a contact.
    let mut relaxed = spec.clone();
    relaxed.filter.min_mapq = 0;
    let (mut want, _, _) = contacts_any_mq(&f, &relaxed);
    want.retain(|w| w.ends.0 != w.ends.4 || w.ends.5 - w.ends.1 >= 2_000);
    assert_eq!(
        gunzip(&out.join("hic/contacts.pairs.gz")),
        pairs_text(&f, &want, "unknown", false)
    );
    assert_eq!(
        std::fs::read(out.join("hic/contacts.hic")).unwrap(),
        hic_bytes(&f, &want, "unknown")
    );
    assert!(!out.join("hic/cell_metrics.csv").exists());
}

/// [`contacts`] where a missing `MQ` passes (a threshold of 0), its MAPQ
/// written as 255.
fn contacts_any_mq(f: &Fixture, c: &Contacts) -> (Vec<Want>, Vec<String>, u64) {
    let mut g = Fixture {
        dir: tempfile::tempdir().unwrap(),
        refs: f.refs.clone(),
        reads: f.reads.clone(),
        bam: f.bam.clone(),
    };
    for r in &mut g.reads {
        if r.tag("MQ").is_none() {
            r.tags.push("MQ:i:255".into());
        }
    }
    contacts(&g, c)
}

#[test]
fn contacts_do_not_depend_on_threads_or_spilling() {
    let f = fixture(83);
    let a = f.dir.path().join("a");
    quantify(task(&f, &a, sc_contacts()), Budget::new(1, 1 << 30), None).unwrap();
    let b = f.dir.path().join("b");
    let mut p = plan(task(&f, &b, sc_contacts()), Budget::new(7, 256 << 20)).unwrap();
    for s in &mut p.stages {
        s.1 = 16 << 10;
    }
    run(&p, None).unwrap();
    let summary = std::fs::read_to_string(b.join("summary.json")).unwrap();
    assert!(
        summary
            .lines()
            .any(|l| l.contains("spilled_bytes") && !l.contains(": 0,")),
        "the small buffers spill"
    );
    let strip = |v: Vec<(String, Vec<u8>)>| -> Vec<(String, Vec<u8>)> {
        v.into_iter().filter(|(n, _)| n != "summary.json").collect()
    };
    assert_eq!(strip(files(&a)), strip(files(&b)));
}
