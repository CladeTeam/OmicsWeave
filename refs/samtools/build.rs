//! Builds the vendored samtools commands against `htslib-sys`, without
//! samtools' configure script.

use std::env;
use std::path::PathBuf;

/// The release vendored; `samtools_version()` reports it.
const VERSION: &str = "1.25";

/// The command implementations and the samtools utilities they call.
const SOURCES: &[&str] = &[
    "bam_sort.c",
    "bam2depth.c",
    "bam.c",
    "bedidx.c",
    "sam_opts.c",
    "sam_utils.c",
];

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = root.join("vendor/samtools");
    // The sources include the configure-generated config.h; none of the
    // settings it would hold (curses for tview) concern these files.
    let generated = PathBuf::from(env::var("OUT_DIR").unwrap()).join("include");
    std::fs::create_dir_all(&generated).unwrap();
    std::fs::write(
        generated.join("config.h"),
        "/* samtools without configure: no optional features. */\n",
    )
    .unwrap();
    let mut build = cc::Build::new();
    build
        .include(&generated)
        .include(&vendor)
        .include(env::var("DEP_HTS_INCLUDE").unwrap())
        .include(env::var("DEP_HTS_GENERATED").unwrap())
        .define(
            "OW_SAMTOOLS_VERSION",
            Some(format!("\"{VERSION}\"").as_str()),
        )
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false);
    if let Ok(include) = env::var("DEP_Z_INCLUDE") {
        build.include(include);
    }
    for file in SOURCES {
        build.file(vendor.join(file));
    }
    build.file(root.join("src/shim.c"));
    build.compile("samtools");
    println!("cargo:rerun-if-changed=vendor/samtools");
    println!("cargo:rerun-if-changed=src/shim.c");
}
