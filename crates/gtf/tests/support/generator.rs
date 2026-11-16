//! Deterministic annotation generator for differential tests.

/// A tiny deterministic generator (SplitMix64).
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

/// A random GTF with no gene / transcript records for some genes, exons
/// in random order, overlapping genes on both strands; coordinates below
/// `len` on chromosomes `c0`.. so per-base bitmaps can check views.
pub fn random_gtf(rng: &mut Rng, chroms: usize, genes_per_chrom: usize, len: u32) -> String {
    let mut out = String::from("##format: gtf\n");
    let mut lines = Vec::new();
    for c in 0..chroms {
        for g in 0..genes_per_chrom {
            let gid = format!("G{c}_{g}");
            let strand = if rng.below(2) == 0 { '+' } else { '-' };
            let n_tx = 1 + rng.below(3) as usize;
            let base = rng.below(len as u64 - 200) as u32;
            let mut gs = u32::MAX;
            let mut ge = 0;
            let mut tx_lines = Vec::new();
            for t in 0..n_tx {
                let tid = format!("T{c}_{g}_{t}");
                let n_ex = 1 + rng.below(4) as usize;
                let mut at = base + rng.below(20) as u32;
                let mut exons = Vec::new();
                for _ in 0..n_ex {
                    let l = 1 + rng.below(30) as u32;
                    let s = at;
                    let e = (s + l).min(len);
                    if s >= e {
                        break;
                    }
                    exons.push((s, e));
                    at = e + rng.below(40) as u32;
                }
                let ts = exons[0].0;
                let te = exons.last().unwrap().1;
                gs = gs.min(ts);
                ge = ge.max(te);
                let attrs = format!(
                    "gene_id \"{gid}\"; transcript_id \"{tid}\"; transcript_type \"tt{}\";",
                    t % 2
                );
                if rng.below(2) == 0 {
                    tx_lines.push(format!(
                        "c{c}\tsrc\ttranscript\t{}\t{te}\t.\t{strand}\t.\t{attrs}",
                        ts + 1
                    ));
                }
                for &(s, e) in &exons {
                    tx_lines.push(format!(
                        "c{c}\tsrc\texon\t{}\t{e}\t.\t{strand}\t.\t{attrs}",
                        s + 1
                    ));
                }
                let (cs, ce) = exons[0];
                if rng.below(2) == 0 {
                    tx_lines.push(format!(
                        "c{c}\tsrc\tCDS\t{}\t{ce}\t.\t{strand}\t0\t{attrs}",
                        cs + 1
                    ));
                }
            }
            if rng.below(2) == 0 {
                lines.push(format!(
                    "c{c}\tsrc\tgene\t{}\t{ge}\t.\t{strand}\t.\tgene_id \"{gid}\"; gene_type \"bt{}\";",
                    gs + 1,
                    g % 3
                ));
            }
            lines.extend(tx_lines);
        }
    }
    // Shuffle: unsorted input must be fine.
    for i in (1..lines.len()).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        lines.swap(i, j);
    }
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out
}
