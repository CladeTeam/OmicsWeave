use std::env;
use std::path::PathBuf;

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = root.join("vendor/libBigWig");

    let mut build = cc::Build::new();
    build
        .include(&vendor)
        .define("NOCURL", None)
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false);
    for file in ["bwRead.c", "bwStats.c", "bwValues.c", "bwWrite.c", "io.c"] {
        build.file(vendor.join(file));
    }
    if let Ok(include) = env::var("DEP_Z_INCLUDE") {
        build.include(include);
    }
    build.compile("bigwig");
    println!("cargo:rustc-link-lib=m");

    println!("cargo:rerun-if-changed=vendor/libBigWig");
}
