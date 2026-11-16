use super::alignment::Read;
use super::fixtures::Gene;
use super::output::{gunzip_bytes, parse_tbi};
use super::support::cstr;
use std::collections::BTreeMap;
use std::ffi::{CString, c_char, c_int};
use std::path::Path;
use std::sync::Mutex;

/// C command entry points share getopt and other process-global state.
static TOOLS: Mutex<()> = Mutex::new(());

fn args(v: &[String]) -> (Vec<CString>, Vec<*mut c_char>) {
    let owned: Vec<CString> = v
        .iter()
        .map(|a| CString::new(a.as_str()).unwrap())
        .collect();
    let mut ptrs: Vec<*mut c_char> = owned.iter().map(|a| a.as_ptr() as *mut c_char).collect();
    ptrs.push(std::ptr::null_mut());
    (owned, ptrs)
}

unsafe extern "C" {
    static mut optind: c_int;
}

/// `samtools <args>`.
pub fn samtools(a: &[String]) -> i32 {
    let _g = TOOLS.lock().unwrap_or_else(|e| e.into_inner());
    let (_o, mut p) = args(a);
    unsafe {
        optind = 0;
        samtools_sys::ow_samtools(a.len() as c_int, p.as_mut_ptr())
    }
}

/// `featureCounts <args>`.
pub fn feature_counts(a: &[String]) -> i32 {
    let _g = TOOLS.lock().unwrap_or_else(|e| e.into_inner());
    let mut full = vec!["featureCounts".to_string()];
    full.extend_from_slice(a);
    let (_o, mut p) = args(&full);
    unsafe {
        optind = 0;
        subread_sys::ow_featurecounts_main(full.len() as c_int, p.as_mut_ptr())
    }
}

/// `samtools depth -aa <extra> bam`: per reference, the depth of every base.
pub fn samtools_depth(bam: &Path, extra: &[&str], refs: &[(String, u32)]) -> Vec<Vec<u32>> {
    let out = bam.with_extension("depth.txt");
    let mut a: Vec<String> = vec!["depth".into(), "-aa".into(), "-o".into()];
    a.push(out.to_str().unwrap().into());
    a.extend(extra.iter().map(|s| s.to_string()));
    a.push(bam.to_str().unwrap().into());
    assert_eq!(samtools(&a), 0, "samtools depth failed");
    let mut depth: Vec<Vec<u32>> = refs.iter().map(|r| vec![0; r.1 as usize]).collect();
    for line in std::fs::read_to_string(&out).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let t = refs.iter().position(|r| r.0 == f[0]).unwrap();
        let p: usize = f[1].parse().unwrap();
        depth[t][p - 1] = f[2].parse().unwrap();
    }
    depth
}

/// Whether `[s, e)` overlaps any of `ivs`.
pub fn overlaps(ivs: &[(u32, u32)], s: u32, e: u32) -> bool {
    ivs.iter().any(|&(a, b)| a < e && s < b)
}

/// The genes (indices) whose exons the blocks overlap, on `strand` when
/// given (`Some(true)`: reverse).
pub fn gene_hits(
    genes: &[Gene],
    tid: usize,
    blocks: &[(u32, u32)],
    strand: Option<bool>,
) -> Vec<usize> {
    genes
        .iter()
        .enumerate()
        .filter(|(_, g)| g.tid == tid)
        .filter(|(_, g)| strand.is_none_or(|rev| (g.strand == '-') == rev))
        .filter(|(_, g)| {
            let ex = g.exonic();
            blocks.iter().any(|&(s, e)| overlaps(&ex, s, e))
        })
        .map(|(i, _)| i)
        .collect()
}

/// The fragment strand of a read under a library type (0 / 1 / 2).
pub fn reverse_for(r: &Read, lib: u8) -> bool {
    let read2 = r.flag & 0x81 == 0x81;
    (r.reverse() ^ read2) ^ (lib == 2)
}

/// htslib's `tbx_conf_t`.
#[repr(C)]
pub struct TbxConf {
    preset: i32,
    sc: i32,
    bc: i32,
    ec: i32,
    meta_char: i32,
    line_skip: i32,
}

unsafe extern "C" {
    static tbx_conf_bed: TbxConf;
    fn tbx_index_build(fname: *const c_char, min_shift: c_int, conf: *const TbxConf) -> c_int;
}

/// Checks that `<gz>.tbi` is the index htslib's `tbx_index_build` (BED
/// preset) writes for a copy of `gz` in `scratch` (bins in any order).
pub fn assert_tabix_as_htslib(gz: &Path, scratch: &Path) {
    let copy = scratch.join(format!(
        "copy-{}",
        gz.file_name().unwrap().to_str().unwrap()
    ));
    std::fs::copy(gz, &copy).unwrap();
    {
        let _g = TOOLS.lock().unwrap_or_else(|e| e.into_inner());
        assert_eq!(
            unsafe { tbx_index_build(cstr(&copy).as_ptr(), 0, &tbx_conf_bed) },
            0
        );
    }
    let mut ours = gz.as_os_str().to_owned();
    ours.push(".tbi");
    let mut theirs = copy.as_os_str().to_owned();
    theirs.push(".tbi");
    assert_eq!(
        parse_tbi(&gunzip_bytes(Path::new(&ours))),
        parse_tbi(&gunzip_bytes(Path::new(&theirs)))
    );
}
