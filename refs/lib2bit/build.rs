use std::env;
use std::path::PathBuf;

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = root.join("vendor/lib2bit");
    cc::Build::new()
        .include(&vendor)
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false)
        .file(vendor.join("2bit.c"))
        .compile("2bit");
    println!("cargo:rerun-if-changed=vendor/lib2bit");
}
