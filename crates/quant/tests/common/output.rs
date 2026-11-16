use omicsweave_core::{FileSource, TrackRead};
use std::collections::BTreeMap;
use std::io::Read as _;
use std::path::Path;

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
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!((dims[0], dims[1]), (feats.len(), rows.len()));
        let mut n = 0;
        for l in lines {
            let v: Vec<usize> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            assert_eq!(v.len(), 3, "{}: Matrix Market entry {l:?}", dir.display());
            assert!(
                (1..=feats.len()).contains(&v[0]) && (1..=rows.len()).contains(&v[1]),
                "{}: Matrix Market coordinate {v:?}",
                dir.display()
            );
            assert!(
                out.insert(
                    (rows[v[1] - 1].clone(), feats[v[0] - 1].clone()),
                    v[2] as u32,
                )
                .is_none(),
                "{}: duplicate matrix coordinate {v:?}",
                dir.display()
            );
            n += 1;
        }
        assert_eq!(n, dims[2]);
    } else {
        let text = gunzip(&file("csv.gz"));
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split(',').skip(1).collect();
        assert_eq!(
            head,
            rows.iter().map(String::as_str).collect::<Vec<_>>(),
            "{}: CSV row labels",
            dir.display()
        );
        let mut n_features = 0;
        for l in lines {
            let f: Vec<&str> = l.split(',').collect();
            assert_eq!(f.len(), rows.len() + 1, "{}: CSV row {l:?}", dir.display());
            assert_eq!(
                f[0],
                feats[n_features],
                "{}: CSV feature order",
                dir.display()
            );
            n_features += 1;
            for (i, v) in f[1..].iter().enumerate() {
                let v: u32 = v.parse().unwrap();
                if v > 0 {
                    out.insert((rows[i].clone(), f[0].to_owned()), v);
                }
            }
        }
        assert_eq!(n_features, feats.len(), "{}: CSV features", dir.display());
    }
    out
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
