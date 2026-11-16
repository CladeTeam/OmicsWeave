//! Independent index structure and BGZF block inspection.

use std::collections::BTreeMap;
use std::path::Path;

/// One reference of an index: bins (loff, chunks) and the linear index.
#[derive(Debug, PartialEq, Eq, Default)]
pub struct IdxRef {
    pub bins: BTreeMap<u32, (u64, Vec<(u64, u64)>)>,
    pub linear: Vec<u64>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Idx {
    pub head: Vec<u8>,
    pub refs: Vec<IdxRef>,
    pub n_no_coor: Option<u64>,
}

/// The text of a BGZF file.
pub fn inflate(raw: &[u8]) -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    omicsweave_bgzf::Reader::new(raw, Box::new(omicsweave_core::exec::Inline::new()))
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// Parses BAI or CSI bytes (bins in any order).
pub fn parse_index(raw: &[u8]) -> Idx {
    let bytes = if raw.starts_with(b"BAI\x01") {
        raw.to_vec()
    } else {
        inflate(raw)
    };
    let mut at = 0usize;
    let mut take = |n: usize| {
        let s = bytes.get(at..at + n).map(<[u8]>::to_vec);
        at += n;
        s.unwrap_or_default()
    };
    let u32_ = |b: Vec<u8>| u32::from_le_bytes(b.try_into().unwrap());
    let u64_ = |b: Vec<u8>| u64::from_le_bytes(b.try_into().unwrap());
    let magic = take(4);
    let csi = magic == b"CSI\x01";
    let mut head = magic.clone();
    if csi {
        head.extend(take(8));
        let l_aux = u32_(take(4)) as usize;
        head.extend(take(l_aux));
    }
    let n_ref = u32_(take(4));
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let mut r = IdxRef::default();
        for _ in 0..u32_(take(4)) {
            let bin = u32_(take(4));
            let loff = if csi { u64_(take(8)) } else { 0 };
            let n = u32_(take(4));
            let chunks = (0..n).map(|_| (u64_(take(8)), u64_(take(8)))).collect();
            r.bins.insert(bin, (loff, chunks));
        }
        if !csi {
            for _ in 0..u32_(take(4)) {
                r.linear.push(u64_(take(8)));
            }
        }
        refs.push(r);
    }
    let tail = take(8);
    let n_no_coor = (tail.len() == 8).then(|| u64_(tail));
    Idx {
        head,
        refs,
        n_no_coor,
    }
}

/// Splits a BGZF file into its blocks.
pub fn blocks(raw: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut at = 0;
    while at < raw.len() {
        let len = u16::from_le_bytes([raw[at + 16], raw[at + 17]]) as usize + 1;
        out.push(&raw[at..at + len]);
        at += len;
    }
    out
}

/// A copy of `bam` with an empty BGZF block after every `every`th block.
pub fn with_empty_blocks(bam: &Path, out: &Path, every: usize) {
    let raw = std::fs::read(bam).unwrap();
    let mut spliced = Vec::new();
    for (i, b) in blocks(&raw).into_iter().enumerate() {
        spliced.extend_from_slice(b);
        if i % every == 0 {
            spliced.extend_from_slice(&omicsweave_bgzf::EOF_BLOCK);
        }
    }
    std::fs::write(out, spliced).unwrap();
}
