//! Builds the vendored htscodecs without its autotools configure step. The
//! generated `config.h` states what configure finds on a current GCC or
//! Clang: the x86-64 SIMD rANS decoders (SSE4.1, AVX2, AVX-512) are each
//! compiled with their own instruction-set flags and chosen at run time
//! by CPUID, as upstream does; on AArch64 the NEON decoder is used. bzip2
//! (for the arithmetic coder's external mode) comes from `bzip2-sys`.

use std::env;
use std::fs;
use std::path::PathBuf;

const VERSION: &str = "1.6.7";

const SOURCES: &[&str] = &[
    "arith_dynamic.c",
    "fqzcomp_qual.c",
    "htscodecs.c",
    "pack.c",
    "rANS_static4x16pr.c",
    "rANS_static32x16pr.c",
    "rANS_static.c",
    "rle.c",
    "tokenise_name3.c",
    "utils.c",
];

/// x86-64 decoders and the flags each needs.
const X86_SIMD: &[(&str, &[&str])] = &[
    (
        "rANS_static32x16pr_sse4.c",
        &["-msse4.1", "-mssse3", "-mpopcnt"],
    ),
    ("rANS_static32x16pr_avx2.c", &["-mavx2", "-mpopcnt"]),
    ("rANS_static32x16pr_avx512.c", &["-mavx512f", "-mpopcnt"]),
];

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let src = root.join("vendor/htscodecs/htscodecs");
    let generated = PathBuf::from(env::var("OUT_DIR").unwrap()).join("include");
    fs::create_dir_all(&generated).unwrap();
    let arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let unix = env::var("CARGO_CFG_TARGET_FAMILY").is_ok_and(|f| f == "unix");
    let x86 = arch == "x86_64";

    let mut config = String::from(
        "#define HAVE_LIBBZ2 1\n#define HAVE_BUILTIN_PREFETCH 1\n\
         #define HAVE_STDINT_H 1\n#define HAVE_STDLIB_H 1\n#define HAVE_STRING_H 1\n\
         #define HAVE_LIMITS_H 1\n#define HAVE_INTTYPES_H 1\n",
    );
    if unix {
        config.push_str("#define HAVE_UNISTD_H 1\n#define HAVE_FCNTL_H 1\n");
    }
    if x86 {
        config.push_str(
            "#define HAVE_SSE4_1 1\n#define HAVE_SSSE3 1\n#define HAVE_POPCNT 1\n\
             #define HAVE_AVX2 1\n#define HAVE_AVX512 1\n\
             #define HAVE_DECL___CPUID_COUNT 1\n#define HAVE_DECL___GET_CPUID_MAX 1\n",
        );
    }
    fs::write(generated.join("config.h"), config).unwrap();
    fs::write(
        generated.join("version.h"),
        format!("#define HTSCODECS_VERSION_TEXT \"{VERSION}\"\n"),
    )
    .unwrap();

    let base = || {
        let mut build = cc::Build::new();
        build
            .include(&generated)
            .include(&src)
            // Vendored upstream code is compiled as is; its warnings are not ours.
            .warnings(false);
        if let Ok(include) = env::var("DEP_BZIP2_INCLUDE") {
            build.include(include);
        }
        build
    };

    let mut build = base();
    for file in SOURCES {
        build.file(src.join(file));
    }
    if arch == "aarch64" {
        build.file(src.join("rANS_static32x16pr_neon.c"));
    }
    if x86 {
        // One archive: the SIMD objects and the dispatcher refer to each other.
        for (file, flags) in X86_SIMD {
            let mut simd = base();
            for flag in *flags {
                simd.flag(flag);
            }
            simd.file(src.join(file));
            for object in simd.compile_intermediates() {
                build.object(object);
            }
        }
    }
    build.compile("htscodecs");
    if unix {
        // The codec archive uses libm and pthread_once / thread-local pools.
        println!("cargo:rustc-link-lib=m");
        println!("cargo:rustc-link-lib=pthread");
    }
    println!("cargo:include={}", src.display());
    println!("cargo:rerun-if-changed=vendor/htscodecs");
}
