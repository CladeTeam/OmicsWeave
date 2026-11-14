# straw (reference implementation)

- Upstream: https://github.com/aidenlab/straw
- Version: commit `82fba9cf7f323a432bad4ac9a33d1ee65de00a47` (2025-04-05)
- License: MIT (`vendor/straw/LICENSE`)
- Local changes: sources unmodified. Only `C++/straw.cpp`, `C++/straw.h`, `C++/hic_slice.h`,
  `LICENSE` and `README.md` were kept; the R, MATLAB and Python bindings,
  examples and build scripts were removed.

This crate (`straw-sys`) is a dev-dependency of `omicsweave-hic` only. straw
is the Aiden Lab's reader, the one behind `hicstraw` and the Juicer tools'
dump path. It provides:

- the source and its build (`build.rs`, C++14);
- `src/curl-stub/curl/curl.h`, a stand-in for libcurl, which straw includes
  unconditionally for remote files: every call fails, so only local files
  are read and libcurl is not needed;
- `src/shim.cpp`, a C ABI over straw's C++ entry point (C strings in, the
  record vector out), and raw bindings to it in `src/lib.rs`.

There is no test logic here; tests live in `crates/hic/tests/`.

Build-time patch (the vendored file stays unmodified; `build.rs` compiles a
patched copy from `OUT_DIR`): `MatrixZoomData::getRecords` looks up block
index entries with `std::map::operator[]` inside thread-pool tasks. For a
candidate block that the file does not store, this inserts into the shared map
from several threads at once and crashes (SIGSEGV). The patch replaces it with
a read-only lookup that treats an absent block as empty, which is what the
inserted default meant.

Known behaviour left as is: blocks are inflated into a buffer of ten times
their compressed size, so a block that compresses better than 10:1 is
truncated. Test inputs compared against straw stay below that ratio.
