//! Compile-once/reused Filter evaluation against the independent filter
//! baseline. Small fixtures remain within the original recursive stack.

#[path = "../tests/cfilter/mod.rs"]
mod cfilter;

use std::hint::black_box;
use std::time::{Duration, Instant};

use cfilter::{Outcome, Run};
use omicsweave_vcf::{Filter, Header, Record, VariantReader};

fn fixture(id: &str) -> Vec<u8> {
    if let Some(spec) = id.strip_prefix("equality:") {
        let (width, position) = spec.split_once(':').unwrap();
        let width: usize = width.parse().unwrap();
        let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n##INFO=<ID=X,Number=.,Type=Integer>\n##INFO=<ID=Y,Number=.,Type=Integer>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_owned();
        let mut x: Vec<_> = (0..width).map(|i| (i * 2).to_string()).collect();
        let mut y: Vec<_> = (0..width).map(|i| (i * 2 + 1).to_string()).collect();
        match position {
            "miss" => {}
            "early" => x[0] = y[0].clone(),
            "row" => x[0] = y[width - 1].clone(),
            "late" => x[width - 1] = y[width - 1].clone(),
            "uniform" => {
                x.fill("20".to_owned());
                y.fill("10".to_owned());
            }
            "same" => {
                x.fill("20".to_owned());
                y.fill("20".to_owned());
            }
            _ => panic!("unknown equality fixture"),
        }
        let (x, y) = (x.join(","), y.join(","));
        for i in 0..128 {
            text.push_str(&format!("chr1\t{}\t.\tA\tC\t.\t.\tX={x};Y={y}\n", i + 1));
        }
        return text.into_bytes();
    }
    if let Some(spec) = id.strip_prefix("vectors:") {
        let (width, position) = spec.split_once(':').unwrap();
        let width: usize = width.parse().unwrap();
        let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n##INFO=<ID=X,Number=.,Type=Integer>\n##INFO=<ID=Y,Number=.,Type=Integer>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_owned();
        let mut x = vec!["20"; width];
        let y = vec!["10"; width].join(",");
        match position {
            "early" => x[0] = "0",
            "late" => x[width - 1] = "0",
            "miss" => {}
            _ => panic!("unknown vector fixture"),
        }
        let x = x.join(",");
        for i in 0..128 {
            text.push_str(&format!("chr1\t{}\t.\tA\tC\t.\t.\tX={x};Y={y}\n", i + 1));
        }
        return text.into_bytes();
    }
    let mut text = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n##INFO=<ID=DP,Number=1,Type=Integer>\n##INFO=<ID=TXT,Number=1,Type=String>\n##FORMAT=<ID=DP,Number=1,Type=Integer>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_owned();
    for i in 0..32 {
        text.push_str(&format!("\tS{i}"));
    }
    text.push('\n');
    for i in 0..128 {
        let dp = if i % 2 == 0 { 30 } else { 10 };
        text.push_str(&format!(
            "chr1\t{}\t{id}\tA\tC\t50\tPASS\tDP={dp};TXT=sample\tDP",
            i + 1
        ));
        for _ in 0..32 {
            text.push_str(&format!("\t{dp}"));
        }
        text.push('\n');
    }
    text.into_bytes()
}

fn digest(values: impl Iterator<Item = bool>) -> (usize, u64) {
    let (mut passes, mut hash) = (0, 0u64);
    for value in values {
        passes += usize::from(value);
        hash = hash
            .wrapping_mul(0x100000001b3)
            .wrapping_add(u64::from(value));
    }
    (passes, hash)
}

fn scan(filter: &mut Filter, h: &Header, records: &[Record]) -> (usize, u64) {
    digest(
        records
            .iter()
            .map(|rec| filter.test(h, black_box(rec)).unwrap()),
    )
}

fn measure(mut f: impl FnMut()) -> f64 {
    f();
    let mut groups = [0.0f64; 5];
    for sample in &mut groups {
        let start = Instant::now();
        let mut n = 0;
        while start.elapsed() < Duration::from_millis(100) {
            f();
            n += 1;
        }
        *sample = start.elapsed().as_secs_f64() * 1e9 / n as f64;
    }
    groups.sort_by(f64::total_cmp);
    groups[2]
}

fn main() {
    let selected = std::env::args().skip(1).find(|arg| !arg.starts_with('-'));
    let dir = tempfile::tempdir().unwrap();
    for (name, id, expr) in [
        (
            "eq256-miss",
            "equality:256:miss".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq1024-miss",
            "equality:1024:miss".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq256-late",
            "equality:256:late".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq256-early",
            "equality:256:early".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq256-row",
            "equality:256:row".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq256-uniform",
            "equality:256:uniform".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq64-miss",
            "equality:64:miss".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "eq4-miss",
            "equality:4:miss".to_owned(),
            "INFO/X==INFO/Y".to_owned(),
        ),
        (
            "ne256-same",
            "equality:256:same".to_owned(),
            "INFO/X!=INFO/Y".to_owned(),
        ),
        (
            "ne256-early",
            "equality:256:miss".to_owned(),
            "INFO/X!=INFO/Y".to_owned(),
        ),
        (
            "vectors256-miss",
            "vectors:256:miss".to_owned(),
            "INFO/X<INFO/Y".to_owned(),
        ),
        (
            "vectors256-late",
            "vectors:256:late".to_owned(),
            "INFO/X<INFO/Y".to_owned(),
        ),
        (
            "vectors256-early",
            "vectors:256:early".to_owned(),
            "INFO/X<INFO/Y".to_owned(),
        ),
        (
            "vectors16-miss",
            "vectors:16:miss".to_owned(),
            "INFO/X<INFO/Y".to_owned(),
        ),
        (
            "vectors4-miss",
            "vectors:4:miss".to_owned(),
            "INFO/X<INFO/Y".to_owned(),
        ),
        ("literal", "rs12345".to_owned(), "ID~\"rs\"".to_owned()),
        (
            "literal-miss",
            format!("{}b", "a".repeat(256)),
            format!("ID~\"{}c\"", "a".repeat(32)),
        ),
        (
            "atomic-repeat",
            format!("{}b", "a".repeat(32)),
            "ID~\"^a*b$\"".to_owned(),
        ),
        ("atomic-miss", "a".repeat(32), "ID~\"^a*b$\"".to_owned()),
        (
            "numeric-string",
            "rs12345".to_owned(),
            "INFO/DP>=20 && INFO/TXT~\"sample\"".to_owned(),
        ),
        (
            "format-numeric",
            "rs12345".to_owned(),
            "FMT/DP>=20".to_owned(),
        ),
        (
            "backref-hit",
            "aa".to_owned(),
            "ID~\"^\\(a\\)\\1$\"".to_owned(),
        ),
        (
            "backref-miss",
            "ab".to_owned(),
            "ID~\"^\\(a\\)\\1$\"".to_owned(),
        ),
        (
            "backref-alt-hit",
            "aa".to_owned(),
            "ID~\"^\\(a\\|b\\)\\1$\"".to_owned(),
        ),
        (
            "backref-alt-miss",
            "ab".to_owned(),
            "ID~\"^\\(a\\|b\\)\\1$\"".to_owned(),
        ),
    ] {
        if selected.as_deref().is_some_and(|s| s != name) {
            continue;
        }
        let bytes = fixture(&id);
        let path = dir.path().join(format!("{name}.vcf"));
        std::fs::write(&path, &bytes).unwrap();
        let expected = cfilter::c_filter(&path, &expr);
        assert_eq!(
            cfilter::ours(&bytes, &expr),
            expected,
            "{name}: bcftools oracle"
        );
        let Run::Records(expected) = expected else {
            panic!("fixture must compile")
        };
        let want = digest(expected.iter().map(|result| match result {
            Outcome::Pass(b) => *b,
            other => panic!("defined fixture required: {other:?}"),
        }));
        let mut reader = VariantReader::new(&bytes[..]).unwrap();
        let mut record = Record::default();
        let mut records = Vec::new();
        while reader.read(&mut record).unwrap() {
            records.push(record.clone());
        }
        let h = reader.header();
        let mut filter = Filter::new(h, &expr).unwrap();
        assert_eq!(scan(&mut filter, h, &records), want);
        let eval_ns = measure(|| {
            assert_eq!(black_box(scan(&mut filter, h, &records)), want);
        });
        let compile_ns = measure(|| {
            black_box(Filter::new(black_box(h), black_box(&expr)).unwrap());
        });
        println!(
            "case={name} records={} pass={} digest={:016x} compile_ns={compile_ns:.0} reused_scan_ns={eval_ns:.0} ns_per_record={:.2}",
            records.len(),
            want.0,
            want.1,
            eval_ns / records.len() as f64
        );
    }
}
