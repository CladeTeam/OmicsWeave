# lib2bit (reference implementation)

- Upstream: https://github.com/dpryan79/lib2bit
- Version: tag `1.0.1` (commit `cbd7e79174e3d0944ca7ea7984efe7675a45762f`)
- Archive: `https://github.com/dpryan79/lib2bit/archive/refs/tags/1.0.1.tar.gz`,
  sha256 `d526163459c43976e5e02d0fd9c8e9c89e9fcbe7fb4f24fa77b0920a6a017aae`
- License: MIT (`vendor/lib2bit/LICENSE`)
- Local changes: sources unmodified. Only `2bit.c`, `2bit.h`, `LICENSE` and
  `README.md` were kept; docs, tests and the Makefile were removed.

This crate (`lib2bit-sys`) is a dev-dependency only: the source, its build and
raw bindings, without test logic. lib2bit is the library under py2bit and was
written by the author of libBigWig. It reads 2bit files but cannot write them;
2bit tests generate their inputs in test code.
