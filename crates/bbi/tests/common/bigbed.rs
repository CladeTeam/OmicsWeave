//! A minimal bigBed writer for test inputs. libBigWig cannot create bigBed
//! files, so tests write them here and read them with both libBigWig and the
//! crate; libBigWig is an independent reader, so agreement checks the writer
//! and the crate's reader alike.
//!
//! Section order follows the UCSC tools: header, autoSql, total summary,
//! chromosome tree, data, R-tree, trailing magic.

use omicsweave_bbi::{BIGBED_MAGIC, CHROM_TREE_MAGIC, HEADER_LEN, INDEX_MAGIC};

#[derive(Debug, Clone)]
pub struct Record {
    pub chrom: u32,
    pub start: u32,
    pub end: u32,
    pub rest: Vec<u8>,
}

pub struct Spec {
    pub chroms: Vec<(String, u32)>,
    /// Sorted by (chrom, start).
    pub records: Vec<Record>,
    pub items_per_block: usize,
    /// Start a new block at every chromosome change; otherwise blocks may
    /// span chromosomes.
    pub split_chroms: bool,
    pub chrom_fanout: usize,
    pub index_fanout: usize,
    pub auto_sql: Option<String>,
    pub compress: bool,
    pub field_count: u16,
    pub defined_field_count: u16,
}

fn put_u16(out: &mut Vec<u8>, v: u16) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_u64(out: &mut Vec<u8>, v: u64) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn put_f64(out: &mut Vec<u8>, v: f64) {
    out.extend_from_slice(&v.to_le_bytes());
}

/// Groups `n` items into tree levels, leaves first: each level lists the
/// item ranges of its nodes.
fn levels(n: usize, fanout: usize) -> Vec<Vec<std::ops::Range<usize>>> {
    let mut levels = vec![
        (0..n.max(1))
            .step_by(fanout)
            .map(|s| s..(s + fanout).min(n))
            .collect::<Vec<_>>(),
    ];
    while levels.last().unwrap().len() > 1 {
        let below = levels.last().unwrap().len();
        levels.push(
            (0..below)
                .step_by(fanout)
                .map(|s| s..(s + fanout).min(below))
                .collect(),
        );
    }
    levels
}

/// Writes a tree whose root is at `base`, laid out root level first. `leaf`
/// writes leaf item `i`; `branch` writes the key of the subtree covering leaf
/// items `range` followed by nothing (the child offset is appended here).
fn write_tree(
    out: &mut Vec<u8>,
    n: usize,
    fanout: usize,
    leaf_len: usize,
    branch_len: usize,
    leaf: &dyn Fn(&mut Vec<u8>, usize),
    branch: &dyn Fn(&mut Vec<u8>, std::ops::Range<usize>),
) {
    let levels = levels(n, fanout);
    // Leaf-item range covered by each node of each level.
    let mut covered: Vec<Vec<std::ops::Range<usize>>> = vec![levels[0].clone()];
    for level in &levels[1..] {
        let below = covered.last().unwrap();
        covered.push(
            level
                .iter()
                .map(|r| below[r.start].start..below[r.end - 1].end)
                .collect(),
        );
    }
    let node_len = |depth: usize, node: &std::ops::Range<usize>| {
        4 + node.len() * if depth == 0 { leaf_len } else { branch_len }
    };
    // File offset of each level's first node, root level first.
    let mut level_start = vec![0u64; levels.len()];
    let mut at = out.len() as u64;
    for depth in (0..levels.len()).rev() {
        level_start[depth] = at;
        at += levels[depth]
            .iter()
            .map(|r| node_len(depth, r) as u64)
            .sum::<u64>();
    }
    for depth in (0..levels.len()).rev() {
        let mut child = if depth > 0 { level_start[depth - 1] } else { 0 };
        for node in &levels[depth] {
            out.push((depth == 0) as u8);
            out.push(0);
            put_u16(out, node.len() as u16);
            for i in node.clone() {
                if depth == 0 {
                    leaf(out, i);
                } else {
                    branch(out, covered[depth - 1][i].clone());
                    put_u64(out, child);
                    child += node_len(depth - 1, &levels[depth - 1][i]) as u64;
                }
            }
        }
    }
}

pub fn write(spec: &Spec) -> Vec<u8> {
    // Blocks: ranges of records.
    let mut blocks: Vec<std::ops::Range<usize>> = Vec::new();
    let mut i = 0;
    while i < spec.records.len() {
        let mut j = i + 1;
        while j < spec.records.len()
            && j - i < spec.items_per_block
            && !(spec.split_chroms && spec.records[j].chrom != spec.records[i].chrom)
        {
            j += 1;
        }
        blocks.push(i..j);
        i = j;
    }
    let raw_blocks: Vec<Vec<u8>> = blocks
        .iter()
        .map(|r| {
            let mut b = Vec::new();
            for rec in &spec.records[r.clone()] {
                put_u32(&mut b, rec.chrom);
                put_u32(&mut b, rec.start);
                put_u32(&mut b, rec.end);
                b.extend_from_slice(&rec.rest);
                b.push(0);
            }
            b
        })
        .collect();
    let max_raw = raw_blocks.iter().map(Vec::len).max().unwrap_or(0);
    let stored: Vec<Vec<u8>> = if spec.compress {
        let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::default());
        raw_blocks
            .iter()
            .map(|b| {
                let mut out = vec![0u8; c.zlib_compress_bound(b.len())];
                let n = c.zlib_compress(b, &mut out).expect("compress");
                out.truncate(n);
                out
            })
            .collect()
    } else {
        raw_blocks
    };

    let mut out = vec![0u8; HEADER_LEN];
    let auto_sql_offset = spec.auto_sql.as_ref().map(|sql| {
        let at = out.len() as u64;
        out.extend_from_slice(sql.as_bytes());
        out.push(0);
        at
    });

    let total_summary_offset = out.len() as u64;
    let covered: u64 = spec.records.iter().map(|r| (r.end - r.start) as u64).sum();
    put_u64(&mut out, covered);
    put_f64(&mut out, 1.0);
    put_f64(&mut out, 1.0);
    put_f64(&mut out, covered as f64);
    put_f64(&mut out, covered as f64);

    // Chromosome B+ tree, keys sorted by name as the UCSC tools require.
    let chrom_tree_offset = out.len() as u64;
    let key_len = spec
        .chroms
        .iter()
        .map(|(n, _)| n.len())
        .max()
        .unwrap_or(1)
        .max(1);
    let mut by_name: Vec<usize> = (0..spec.chroms.len()).collect();
    by_name.sort_by(|&a, &b| spec.chroms[a].0.cmp(&spec.chroms[b].0));
    put_u32(&mut out, CHROM_TREE_MAGIC);
    put_u32(&mut out, spec.chrom_fanout as u32);
    put_u32(&mut out, key_len as u32);
    put_u32(&mut out, 8);
    put_u64(&mut out, spec.chroms.len() as u64);
    put_u64(&mut out, 0);
    let key = |out: &mut Vec<u8>, i: usize| {
        let name = spec.chroms[by_name[i]].0.as_bytes();
        out.extend_from_slice(name);
        out.resize(out.len() + key_len - name.len(), 0);
    };
    write_tree(
        &mut out,
        spec.chroms.len(),
        spec.chrom_fanout,
        key_len + 8,
        key_len + 8,
        &|out, i| {
            key(out, i);
            put_u32(out, by_name[i] as u32);
            put_u32(out, spec.chroms[by_name[i]].1);
        },
        &|out, range| key(out, range.start),
    );

    let full_data_offset = out.len() as u64;
    put_u64(&mut out, spec.records.len() as u64);
    let mut refs = Vec::new();
    for b in &stored {
        refs.push((out.len() as u64, b.len() as u64));
        out.extend_from_slice(b);
    }

    let bounds: Vec<((u32, u32), (u32, u32))> = blocks
        .iter()
        .map(|r| {
            let recs = &spec.records[r.clone()];
            let first = (recs[0].chrom, recs[0].start);
            let last = recs.iter().map(|x| (x.chrom, x.end)).max().unwrap();
            (first, last)
        })
        .collect();
    let span = |range: std::ops::Range<usize>| {
        let first = bounds[range.start].0;
        let last = bounds[range].iter().map(|b| b.1).max().unwrap();
        (first, last)
    };

    let full_index_offset = out.len() as u64;
    put_u32(&mut out, INDEX_MAGIC);
    put_u32(&mut out, spec.index_fanout as u32);
    put_u64(&mut out, blocks.len() as u64);
    let (first, last) = if blocks.is_empty() {
        ((0, 0), (0, 0))
    } else {
        span(0..blocks.len())
    };
    for v in [first.0, first.1, last.0, last.1] {
        put_u32(&mut out, v);
    }
    put_u64(&mut out, full_index_offset);
    put_u32(&mut out, spec.items_per_block as u32);
    put_u32(&mut out, 0);
    let bounds_item = |out: &mut Vec<u8>, ((a, b), (c, d)): ((u32, u32), (u32, u32))| {
        for v in [a, b, c, d] {
            put_u32(out, v);
        }
    };
    write_tree(
        &mut out,
        blocks.len(),
        spec.index_fanout,
        32,
        24,
        &|out, i| {
            bounds_item(out, bounds[i]);
            put_u64(out, refs[i].0);
            put_u64(out, refs[i].1);
        },
        &|out, range| bounds_item(out, span(range)),
    );
    put_u32(&mut out, BIGBED_MAGIC);

    let mut h = Vec::with_capacity(HEADER_LEN);
    put_u32(&mut h, BIGBED_MAGIC);
    put_u16(&mut h, 4);
    put_u16(&mut h, 0);
    put_u64(&mut h, chrom_tree_offset);
    put_u64(&mut h, full_data_offset);
    put_u64(&mut h, full_index_offset);
    put_u16(&mut h, spec.field_count);
    put_u16(&mut h, spec.defined_field_count);
    put_u64(&mut h, auto_sql_offset.unwrap_or(0));
    put_u64(&mut h, total_summary_offset);
    put_u32(&mut h, if spec.compress { max_raw as u32 } else { 0 });
    put_u64(&mut h, 0);
    out[..HEADER_LEN].copy_from_slice(&h);
    out
}
