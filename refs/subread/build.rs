//! Builds featureCounts and the Subread library files it links, without
//! Subread's Makefile.

use std::env;
use std::path::PathBuf;

/// The release vendored, as `makefile.version` names it.
const VERSION: &str = "2.2.1";

/// `ALL_LIBS` of `Makefile.Linux`: what every Subread program links.
const LIBRARY: &[&str] = &[
    "core.c",
    "core-junction.c",
    "core-indel.c",
    "sambam-file.c",
    "sublog.c",
    "gene-algorithms.c",
    "hashtable.c",
    "input-files.c",
    "sorted-hashtable.c",
    "gene-value-index.c",
    "exon-algorithms.c",
    "HelperFunctions.c",
    "interval_merge.c",
    "long-hashtable.c",
    "core-bigtable.c",
    "seek-zlib.c",
    "input-blc.c",
];

fn build(vendor: &PathBuf) -> cc::Build {
    let mut build = cc::Build::new();
    // The flags of Makefile.Linux; -fcommon because headers define globals
    // that several files include.
    build
        .include(vendor)
        .define("MAKE_FOR_EXON", None)
        .define("MAKE_STANDALONE", None)
        .define("SUBREAD_VERSION", Some(format!("\"{VERSION}\"").as_str()))
        .define("_FILE_OFFSET_BITS", Some("64"))
        .flag_if_supported("-fcommon")
        // htslib (linked into the same test binaries) defines a global
        // `find_path` too.
        .define("find_path", Some("ow_subread_find_path"))
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false);
    if let Ok(include) = env::var("DEP_Z_INCLUDE") {
        build.include(include);
    }
    build
}

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = root.join("vendor/subread");
    // featureCounts' main() becomes an ordinary function; its library is
    // listed first so that the linker resolves it against the others.
    build(&vendor)
        .define("main", Some("ow_featurecounts_main"))
        .file(vendor.join("readSummary.c"))
        .compile("featurecounts");
    let mut lib = build(&vendor);
    for file in LIBRARY {
        lib.file(vendor.join(file));
    }
    lib.compile("subread");
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rerun-if-changed=vendor/subread");
}
