use std::env;
use std::fs;
use std::path::PathBuf;

/// straw looks up block index entries with `std::map::operator[]` inside
/// tasks running on its thread pool. For a candidate block the file does not
/// store, that inserts into the shared map from several threads at once and
/// crashes. The build compiles a copy whose lookup only reads (an absent
/// block reads as empty, which is what the inserted default meant).
const RACE: &str = "fileName, blockMap[blockNumber], version,";
const READ_ONLY: &str = "fileName, (blockMap.count(blockNumber) ? blockMap.at(blockNumber) : indexEntry{0, 0}), version,";

fn main() {
    let root = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let source = fs::read_to_string(root.join("vendor/straw/straw.cpp")).unwrap();
    assert_eq!(
        source.matches(RACE).count(),
        1,
        "straw.cpp changed; revisit the lookup patch"
    );
    let patched = out.join("straw.cpp");
    fs::write(&patched, source.replace(RACE, READ_ONLY)).unwrap();

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++14")
        // The stand-in curl header comes first; straw only reads local files here.
        .include(root.join("src/curl-stub"))
        .include(root.join("vendor/straw"))
        // Vendored upstream code is compiled as is; its warnings are not ours.
        .warnings(false)
        .file(&patched)
        .file(root.join("src/shim.cpp"));
    if let Ok(include) = env::var("DEP_Z_INCLUDE") {
        build.include(include);
    }
    build.compile("straw");
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rerun-if-changed=vendor/straw");
    println!("cargo:rerun-if-changed=src");
}
