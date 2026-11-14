//! Builds the vendored htslib without running its configure script. The
//! generated headers mirror the defaults htslib's Makefile writes when it is
//! not configured, minus the optional bzip2, lzma and libcurl backends.

use std::env;
use std::fs;
use std::path::PathBuf;

const VERSION: &str = "1.25";

const SOURCES: &[&str] = &[
    "kfunc.c",
    "kstring.c",
    "bcf_sr_sort.c",
    "bgzf.c",
    "errmod.c",
    "faidx.c",
    "header.c",
    "hfile.c",
    "hts.c",
    "hts_expr.c",
    "hts_os.c",
    "md5.c",
    "multipart.c",
    "probaln.c",
    "realn.c",
    "regidx.c",
    "region.c",
    "sam.c",
    "sam_mods.c",
    "simd.c",
    "synced_bcf_reader.c",
    "vcf_sweep.c",
    "tbx.c",
    "textutils.c",
    "thread_pool.c",
    "vcf.c",
    "vcfutils.c",
    "cram/cram_codecs.c",
    "cram/cram_decode.c",
    "cram/cram_encode.c",
    "cram/cram_external.c",
    "cram/cram_index.c",
    "cram/cram_io.c",
    "cram/cram_stats.c",
    "cram/mFILE.c",
    "cram/open_trace_file.c",
    "cram/pooled_alloc.c",
    "cram/string_alloc.c",
    "htscodecs/htscodecs/arith_dynamic.c",
    "htscodecs/htscodecs/fqzcomp_qual.c",
    "htscodecs/htscodecs/htscodecs.c",
    "htscodecs/htscodecs/pack.c",
    "htscodecs/htscodecs/rANS_static4x16pr.c",
    "htscodecs/htscodecs/rANS_static32x16pr.c",
    "htscodecs/htscodecs/rANS_static.c",
    "htscodecs/htscodecs/rle.c",
    "htscodecs/htscodecs/tokenise_name3.c",
    "htscodecs/htscodecs/utils.c",
];

const CONFIG_H: &str = "\
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#define HAVE_DRAND48 1
#if defined __x86_64__ || defined __arm__ || defined __aarch64__
#define HAVE_ATTRIBUTE_CONSTRUCTOR 1
#endif
#if defined __linux__
#define HAVE_GETAUXVAL
#endif
";

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = root.join("vendor/htslib");
    let generated = PathBuf::from(env::var("OUT_DIR").unwrap()).join("include");
    fs::create_dir_all(&generated).unwrap();
    fs::write(generated.join("config.h"), CONFIG_H).unwrap();
    fs::write(
        generated.join("version.h"),
        format!("#define HTS_VERSION_TEXT \"{VERSION}\"\n"),
    )
    .unwrap();
    fs::write(
        generated.join("config_vars.h"),
        "#define HTS_CC \"cc\"\n#define HTS_CPPFLAGS \"\"\n#define HTS_CFLAGS \"\"\n\
         #define HTS_LDFLAGS \"\"\n#define HTS_LIBS \"-lz -lm -lpthread\"\n",
    )
    .unwrap();

    let mut build = cc::Build::new();
    build
        .include(&generated)
        .include(&vendor)
        .include(vendor.join("htscodecs"))
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false);
    if let Ok(include) = env::var("DEP_Z_INCLUDE") {
        build.include(include);
    }
    for file in SOURCES {
        build.file(vendor.join(file));
    }
    build.file(root.join("src/fields.c"));
    build.compile("hts");
    // Dependent reference crates compile against these headers.
    println!("cargo:include={}", vendor.display());
    println!("cargo:generated={}", generated.display());
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rerun-if-changed=vendor/htslib");
    println!("cargo:rerun-if-changed=src/fields.c");
}
