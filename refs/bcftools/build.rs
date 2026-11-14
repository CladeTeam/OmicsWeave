use std::env;
use std::path::PathBuf;

/// bcftools carries its own regidx, whose ABI differs from the htslib one of
/// the same names (32-bit coordinates); its symbols are renamed so both link.
const REGIDX: &[&str] = &[
    "regidx_destroy",
    "regidx_init",
    "regidx_init_string",
    "regidx_insert",
    "regidx_insert_list",
    "regidx_nregs",
    "regidx_overlap",
    "regidx_parse_bed",
    "regidx_parse_reg",
    "regidx_parse_tab",
    "regidx_parse_vcf",
    "regidx_push",
    "regidx_seq_names",
    "regidx_seq_nregs",
    "regitr_copy",
    "regitr_destroy",
    "regitr_init",
    "regitr_loop",
    "regitr_overlap",
    "regitr_reset",
];

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = root.join("vendor/bcftools");
    // filter.c includes the configure-generated config.h; the only setting
    // it reads, ENABLE_PERL_FILTERS, stays undefined (no Perl filters).
    let generated = PathBuf::from(env::var("OUT_DIR").unwrap()).join("include");
    std::fs::create_dir_all(&generated).unwrap();
    std::fs::write(
        generated.join("config.h"),
        "/* bcftools configure defaults: no Perl filters. */\n",
    )
    .unwrap();
    let mut build = cc::Build::new();
    build
        .include(&generated)
        .include(&vendor)
        .include(env::var("DEP_HTS_INCLUDE").unwrap())
        .include(env::var("DEP_HTS_GENERATED").unwrap())
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false);
    for name in REGIDX {
        build.define(name, Some(format!("bcftools_{name}").as_str()));
    }
    for file in ["consensus.c", "filter.c", "regidx.c", "smpl_ilist.c"] {
        build.file(vendor.join(file));
    }
    build.file(root.join("src/shim.c"));
    build.compile("bcftools_consensus");
    println!("cargo:rerun-if-changed=vendor/bcftools");
    println!("cargo:rerun-if-changed=src/shim.c");
}
