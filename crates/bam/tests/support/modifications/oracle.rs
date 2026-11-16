use htslib_sys as h;
use omicsweave_bam::mods::{ModDecoder, ModSeq, ModTags, ModType};
use omicsweave_bam::{Record, Scanner};
use omicsweave_core::{ModificationEntry, exec::Inline};
use std::ffi::CString;
use std::path::Path;

/// A modification as htslib reports it: position, code, canonical base,
/// strand (1 for `-`), ML value (-1 without ML).
pub type HtsMod = (i32, i32, u8, i32, i32);

/// Codes as `bam_mods_recorded` and `bam_mods_queryi` give them: code,
/// canonical base, strand, implicit.
type HtsType = (i32, u8, i32, bool);

fn c(p: &Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

/// htslib's decoding of every record of `bam`, in file order: `None` when
/// parsing or iterating fails.
pub fn htslib_mods(bam: &Path) -> Vec<Option<(Vec<HtsMod>, Vec<HtsType>)>> {
    let mut out = Vec::new();
    unsafe {
        let fp = h::hts_open(c(bam).as_ptr(), c"r".as_ptr());
        assert!(!fp.is_null());
        let hdr = h::sam_hdr_read(fp);
        let b = h::bam_init1();
        let state = h::hts_base_mod_state_alloc();
        while h::sam_read1(fp, hdr, b) >= 0 {
            if h::bam_parse_basemod2(b, state, 0) < 0 {
                out.push(None);
                continue;
            }
            let mut ntype = 0;
            let codes = h::bam_mods_recorded(state, &mut ntype);
            let mut types = Vec::new();
            for i in 0..ntype {
                let (mut strand, mut implicit, mut canonical) = (0, 0, 0 as std::ffi::c_char);
                assert_eq!(
                    h::bam_mods_queryi(state, i, &mut strand, &mut implicit, &mut canonical),
                    0
                );
                types.push((
                    *codes.add(i as usize),
                    canonical as u8,
                    strand,
                    implicit != 0,
                ));
            }
            let mut mods = Vec::new();
            let mut buf = [h::hts_base_mod::default(); 64];
            let mut pos = 0;
            let ok = loop {
                let n = h::bam_next_basemod(b, state, buf.as_mut_ptr(), 64, &mut pos);
                if n < 0 {
                    break false;
                }
                if n == 0 {
                    break true;
                }
                assert!(n <= 64);
                for m in &buf[..n as usize] {
                    mods.push((
                        pos,
                        m.modified_base,
                        m.canonical_base as u8,
                        m.strand,
                        m.qual,
                    ));
                }
            };
            out.push(ok.then_some((mods, types)));
        }
        h::hts_base_mod_state_free(state);
        h::bam_destroy1(b);
        h::sam_hdr_destroy(hdr);
        h::hts_close(fp);
    }
    out
}

pub fn as_hts(rows: &[ModificationEntry]) -> Vec<HtsMod> {
    rows.iter()
        .map(|m| {
            (
                m.qpos,
                m.code,
                m.canonical,
                (m.flags & ModificationEntry::MINUS != 0) as i32,
                if m.flags & ModificationEntry::UNSCORED != 0 {
                    h::HTS_MOD_UNKNOWN
                } else {
                    m.probability as i32
                },
            )
        })
        .collect()
}

pub fn types_as_hts(types: &[ModType]) -> Vec<HtsType> {
    types
        .iter()
        .map(|t| (t.code, t.canonical, t.minus as i32, !t.explicit))
        .collect()
}

/// Our decoding of every record of `bam`, in file order.
pub fn our_mods(bam: &Path) -> Vec<Option<(Vec<ModificationEntry>, Vec<ModType>)>> {
    let raw = std::fs::read(bam).unwrap();
    let mut scan = Scanner::new(&raw[..], Box::new(Inline::new()));
    let header = scan.header().unwrap();
    let mut decoder = ModDecoder::new();
    let mut out = Vec::new();
    while let Some(range) = scan.next_record().unwrap() {
        let rec = Record::parse(scan.slice(range), header.references().len()).unwrap();
        let (a, b) = rec.aux();
        let mut rows = Vec::new();
        let ok = match ModTags::find(&[a, b]) {
            Err(_) => false,
            Ok(None) => {
                decoder = ModDecoder::new();
                true
            }
            Ok(Some(tags)) => decoder
                .decode(
                    &tags,
                    ModSeq::Packed {
                        bytes: rec.packed_seq(),
                        len: rec.l_seq(),
                    },
                    rec.is_reverse(),
                    &[],
                    &mut rows,
                )
                .is_ok(),
        };
        out.push(ok.then(|| (rows, decoder.types().to_vec())));
    }
    out
}
