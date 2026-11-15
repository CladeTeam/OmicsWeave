//! Count matrices against featureCounts (bulk genes, single-end and read
//! pairs) and against reference counts on the generated reads (single-cell
//! UMIs, read groups, whitelists, genome bins, matrix formats).

mod common;

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use common::{
    Gene, Read, Rng, Spec, feature_counts, gene_hits, genes, gtf, read_10x, reverse_for, synth,
    write_bam,
};
use omicsweave_gtf::{Annotation, ReadOptions};
use omicsweave_quant::{
    Budget, CountUnit, Counts, FeatureSource, Filter, MatrixFormat, Product, RowKey, Strandedness,
    Task, quantify,
};

struct Fixture {
    dir: tempfile::TempDir,
    spec: Spec,
    reads: Vec<Read>,
    genes: Vec<Gene>,
    bam: PathBuf,
    gtf: PathBuf,
    ann: Arc<Annotation>,
}

fn fixture(seed: u64, spec: Spec) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = Rng::new(seed);
    let reads = synth(&mut rng, &spec);
    let genes = genes(&mut rng, &spec.refs, 12);
    let bam = write_bam(dir.path(), "in", &spec.refs, &reads, seed);
    let gtf_path = dir.path().join("genes.gtf");
    std::fs::write(&gtf_path, gtf(&spec.refs, &genes)).unwrap();
    let ann = Arc::new(Annotation::read(&gtf_path, &ReadOptions::default()).unwrap());
    Fixture {
        dir,
        spec,
        reads,
        genes,
        bam,
        gtf: gtf_path,
        ann,
    }
}

fn run(
    f: &Fixture,
    out: &str,
    products: Vec<Counts>,
    matrix: MatrixFormat,
    threads: usize,
) -> PathBuf {
    let out = f.dir.path().join(out);
    let mut task = Task::new(vec![f.bam.clone()], &out);
    task.annotation = Some(Arc::clone(&f.ann));
    task.matrix = matrix;
    task.products = products.into_iter().map(Product::Counts).collect();
    quantify(task, Budget::new(threads, 512 << 20), None).unwrap();
    out
}

/// featureCounts' counts by gene id, and its summary.
fn featurecounts(f: &Fixture, extra: &[&str]) -> (HashMap<String, u64>, HashMap<String, u64>) {
    let out = f
        .dir
        .path()
        .join(format!("fc{}.txt", extra.join("").replace('-', "_")));
    let mut a: Vec<String> = [
        "-a",
        f.gtf.to_str().unwrap(),
        "-F",
        "GTF",
        "-T",
        "1",
        "--maxMOp",
        "200",
        "-o",
        out.to_str().unwrap(),
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    a.extend(extra.iter().map(|s| s.to_string()));
    a.push(f.bam.to_str().unwrap().into());
    assert_eq!(feature_counts(&a), 0, "featureCounts failed");
    let mut counts = HashMap::new();
    for l in std::fs::read_to_string(&out).unwrap().lines().skip(2) {
        let v: Vec<&str> = l.split('\t').collect();
        counts.insert(v[0].to_owned(), v[6].parse().unwrap());
    }
    let mut summary = HashMap::new();
    let mut s = out.as_os_str().to_owned();
    s.push(".summary");
    for l in std::fs::read_to_string(PathBuf::from(s))
        .unwrap()
        .lines()
        .skip(1)
    {
        let v: Vec<&str> = l.split('\t').collect();
        summary.insert(v[0].to_owned(), v[1].parse().unwrap());
    }
    (counts, summary)
}

/// Our `counts.tsv`: gene id → count (one input).
fn counts_tsv(dir: &Path) -> HashMap<String, u64> {
    std::fs::read_to_string(dir.join("counts.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let v: Vec<&str> = l.split('\t').collect();
            (v[0].to_owned(), v[6].parse().unwrap())
        })
        .collect()
}

fn summary_of(out: &Path, product: usize) -> omicsweave_quant::Json {
    let text = std::fs::read_to_string(out.join("summary.json")).unwrap();
    // The fields this test reads are plain `"key": number` lines.
    let mut fields = Vec::new();
    let mut seen = 0;
    for line in text.lines() {
        let t = line.trim().trim_end_matches(',');
        if t.starts_with("\"product\"") {
            seen += 1;
        }
        if seen == product + 1
            && let Some((k, v)) = t.split_once(": ")
            && let Ok(n) = v.parse::<u64>()
        {
            fields.push((
                k.trim_matches('"').to_owned(),
                omicsweave_quant::Json::from(n),
            ));
        }
    }
    omicsweave_quant::Json::obj(fields)
}

fn num(j: &omicsweave_quant::Json, k: &str) -> u64 {
    j.get(k)
        .and_then(|v| v.as_u64())
        .unwrap_or_else(|| panic!("no {k}"))
}

/// featureCounts' single-end rules on the generated reads.
fn naive_bulk(f: &Fixture, lib: u8, multi_overlap: bool) -> HashMap<String, u64> {
    let mut out: HashMap<String, u64> = f.genes.iter().map(|g| (g.id.clone(), 0)).collect();
    for r in &f.reads {
        if r.flag & 4 != 0 || r.nh() > 1 {
            continue;
        }
        let strand = (lib > 0).then(|| reverse_for(r, lib));
        let hits = gene_hits(&f.genes, r.tid, &r.blocks(false), strand);
        if hits.len() == 1 || (multi_overlap && !hits.is_empty()) {
            for g in hits {
                *out.get_mut(&f.genes[g].id).unwrap() += 1;
            }
        }
    }
    out
}

#[test]
fn bulk_single_end_matches_featurecounts() {
    let f = fixture(
        3,
        Spec {
            paired: false,
            fragments: 5000,
            ..Spec::default()
        },
    );
    for (lib, overlap) in [(0u8, false), (1, false), (2, false), (0, true)] {
        let mut c = Counts::bulk_rna();
        c.name = format!("g{lib}{overlap}");
        c.strand = Strandedness::from_code(lib).unwrap();
        c.multi_overlap = overlap;
        let out = run(
            &f,
            &format!("se{lib}{overlap}"),
            vec![c],
            MatrixFormat::Cpk,
            4,
        );
        let s = lib.to_string();
        let mut extra = vec!["-s", s.as_str()];
        if overlap {
            extra.push("-O");
        }
        let (fc, fs) = featurecounts(&f, &extra);
        let ours = counts_tsv(&out.join(format!("g{lib}{overlap}")));
        assert_eq!(ours, fc, "counts, -s {lib} -O {overlap}");
        assert_eq!(ours, naive_bulk(&f, lib, overlap), "reference, -s {lib}");
        let j = summary_of(&out, 0);
        assert_eq!(num(&j, "assigned"), fs["Assigned"]);
        assert_eq!(
            num(&j, "unassigned_no_features"),
            fs["Unassigned_NoFeatures"]
        );
        assert_eq!(num(&j, "unassigned_ambiguity"), fs["Unassigned_Ambiguity"]);
        assert_eq!(
            num(&j, "filtered_multimapping"),
            fs["Unassigned_MultiMapping"]
        );
        assert_eq!(num(&j, "filtered_flag"), fs["Unassigned_Unmapped"]);
        // The matrix holds the same counts as the table.
        let m = read_10x(&out.join(format!("g{lib}{overlap}")));
        for ((_, gene), v) in m {
            assert_eq!(ours[&gene], v as u64);
        }
    }
}

#[test]
fn bulk_read_pairs_match_featurecounts() {
    let f = fixture(
        4,
        Spec {
            junk: false,
            fragments: 4000,
            ..Spec::default()
        },
    );
    for lib in [0u8, 2] {
        let mut c = Counts::bulk_rna();
        c.unit = CountUnit::Pair;
        c.strand = Strandedness::from_code(lib).unwrap();
        let out = run(&f, &format!("pe{lib}"), vec![c], MatrixFormat::Mtx, 3);
        let s = lib.to_string();
        let (fc, fs) = featurecounts(&f, &["-p", "--countReadPairs", "-s", &s]);
        let ours = counts_tsv(&out.join("genes"));
        assert_eq!(ours, fc, "pair counts, -s {lib}");
        let j = summary_of(&out, 0);
        assert_eq!(num(&j, "assigned"), fs["Assigned"]);
        assert_eq!(
            num(&j, "unassigned_no_features"),
            fs["Unassigned_NoFeatures"]
        );
        assert_eq!(num(&j, "unassigned_ambiguity"), fs["Unassigned_Ambiguity"]);
    }
}

/// STARsolo-style counts on the generated reads: distinct (barcode, gene,
/// UMI) of unique, forward-stranded reads assigned to one gene.
fn naive_sc(f: &Fixture) -> (Vec<String>, BTreeMap<(String, String), u32>) {
    let mut rows = Vec::new();
    let mut seen = BTreeSet::new();
    let mut umis = BTreeSet::new();
    for r in &f.reads {
        if r.flag & 0x704 != 0 || r.nh() > 1 {
            continue;
        }
        let Some(cb) = r.tag("CB") else { continue };
        if seen.insert(cb.to_owned()) {
            rows.push(cb.to_owned());
        }
        let Some(ub) = r.tag("UB") else { continue };
        let hits = gene_hits(&f.genes, r.tid, &r.blocks(false), Some(reverse_for(r, 1)));
        if hits.len() == 1 {
            umis.insert((cb.to_owned(), f.genes[hits[0]].id.clone(), ub.to_owned()));
        }
    }
    let mut m = BTreeMap::new();
    for (cb, g, _) in umis {
        *m.entry((cb, g)).or_insert(0) += 1;
    }
    (rows, m)
}

#[test]
fn single_cell_umis_match_reference_in_every_format() {
    let f = fixture(5, Spec::default());
    let (rows, want) = naive_sc(&f);
    let mut results = Vec::new();
    for (i, fmt) in [MatrixFormat::Cpk, MatrixFormat::Mtx, MatrixFormat::Csv]
        .into_iter()
        .enumerate()
    {
        let out = run(
            &f,
            &format!("sc{i}"),
            vec![Counts::rna_single_cell()],
            fmt,
            1 + 3 * i,
        );
        let got = read_10x(&out.join("genes"));
        assert_eq!(got, want, "format {fmt:?}");
        let barcodes: Vec<String> = common::gunzip(&out.join("genes/barcodes.tsv.gz"))
            .lines()
            .map(str::to_owned)
            .collect();
        assert_eq!(barcodes, rows);
        let features = common::gunzip(&out.join("genes/features.tsv.gz"));
        let first = features.lines().next().unwrap();
        assert_eq!(first.split('\t').count(), 6);
        assert!(first.contains("\tGene Expression\t"));
        results.push(got);
    }
}

#[test]
fn whitelists_read_groups_and_bins() {
    let f = fixture(6, Spec::default());
    // A whitelist keeps its order and only its barcodes.
    let wl = f.dir.path().join("wl.txt");
    let (rows, want) = naive_sc(&f);
    let listed: Vec<String> = rows.iter().rev().step_by(2).cloned().collect();
    std::fs::write(&wl, listed.join("\n") + "\nNOTSEEN-1\n").unwrap();
    let mut c = Counts::rna_single_cell();
    c.whitelist = Some(wl);
    // By read group, all mapped unique reads.
    let mut rg = Counts::bulk_rna();
    rg.name = "rg".into();
    rg.rows = RowKey::ReadGroup;
    rg.table = false;
    // Genome bins of 1 kb, fragments × file.
    let mut bins = Counts::chip(FeatureSource::Bins(1000));
    bins.name = "bins".into();
    let out = run(&f, "mix", vec![c, rg, bins], MatrixFormat::Cpk, 4);

    let got = read_10x(&out.join("genes"));
    let keep: BTreeSet<&String> = listed.iter().collect();
    let filtered: BTreeMap<(String, String), u32> = want
        .into_iter()
        .filter(|((cb, _), _)| keep.contains(cb))
        .collect();
    assert_eq!(got, filtered);
    let barcodes = common::gunzip(&out.join("genes/barcodes.tsv.gz"));
    assert_eq!(barcodes.lines().count(), listed.len() + 1);
    assert_eq!(barcodes.lines().last(), Some("NOTSEEN-1"));

    let got = read_10x(&out.join("rg"));
    let mut want = BTreeMap::new();
    for r in &f.reads {
        if r.flag & 4 != 0 || r.nh() > 1 {
            continue;
        }
        let hits = gene_hits(&f.genes, r.tid, &r.blocks(false), None);
        if hits.len() == 1 {
            *want
                .entry((r.tag("RG").unwrap().to_owned(), f.genes[hits[0]].id.clone()))
                .or_insert(0) += 1;
        }
    }
    assert_eq!(got, want);

    let got = read_10x(&out.join("bins"));
    let mut want = BTreeMap::new();
    for r in &f.reads {
        if r.flag & Filter::DEFAULT_EXCLUDE != 0 {
            continue;
        }
        let (s, e) = if r.flag & 1 == 0 {
            (r.pos, r.end())
        } else if r.leader() {
            (r.pos, r.pos + r.tlen.unsigned_abs() as u32)
        } else {
            continue;
        };
        let len = f.spec.refs[r.tid].1;
        let chrom = &f.spec.refs[r.tid].0;
        for b in s / 1000..=(e - 1) / 1000 {
            let id = format!("{chrom}:{}-{}", b * 1000, ((b + 1) * 1000).min(len));
            *want.entry((f.bam.display().to_string(), id)).or_insert(0) += 1;
        }
    }
    assert_eq!(got, want);
}

#[test]
fn spilled_sorts_give_the_same_files() {
    let f = fixture(
        8,
        Spec {
            fragments: 20_000,
            barcodes: 400,
            umis: 200,
            ..Spec::default()
        },
    );
    let mut tiles = Counts::rna_single_cell();
    tiles.name = "tiles".into();
    tiles.features = FeatureSource::Bins(200);
    tiles.dedup = omicsweave_quant::Dedup::None;
    let products = [Counts::rna_single_cell(), tiles];
    let a = f.dir.path().join("big");
    let mut task = Task::new(vec![f.bam.clone()], &a);
    task.annotation = Some(Arc::clone(&f.ann));
    task.products = products.iter().cloned().map(Product::Counts).collect();
    quantify(task.clone(), Budget::new(1, 2 << 30), None).unwrap();
    let b = f.dir.path().join("small");
    task.out_dir = b.clone();
    // Sort buffers far smaller than the data, so that both sorts spill.
    let mut p = omicsweave_quant::plan(task, Budget::new(6, 512 << 20)).unwrap();
    for s in &mut p.stages {
        s.1 = 64 << 10;
    }
    omicsweave_quant::run(&p, None).unwrap();
    let summary = std::fs::read_to_string(b.join("summary.json")).unwrap();
    assert_eq!(
        summary
            .lines()
            .filter(|l| l.contains("spilled_bytes") && !l.contains(": 0,"))
            .count(),
        2,
        "both sorts should spill"
    );
    for name in ["genes", "tiles"] {
        for file in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.cpk"] {
            assert_eq!(
                std::fs::read(a.join(name).join(file)).unwrap(),
                std::fs::read(b.join(name).join(file)).unwrap(),
                "{name}/{file}"
            );
        }
    }
}
