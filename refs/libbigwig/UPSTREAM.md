# libBigWig (reference implementation)

- Upstream: https://github.com/dpryan79/libBigWig
- Version: tag `0.4.10` (commit `b5cc843188c93eeef842aa208d13732bc424043f`)
- Archive: `https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.10.tar.gz`,
  sha256 `d7126ee8dc9732f36cd958efccff46e6eafdc1ee7d8da6ff5d43faef962e4d98`
- License: MIT (`vendor/libBigWig/LICENSE`)
- Local changes: sources unmodified. Only the files needed to build and to
  keep provenance were kept from the archive (`*.c`, `*.h`, `LICENSE`,
  `README.md`); `docs/` (generated Doxygen HTML), `test/`, `Doxyfile`,
  `CMakeLists.txt`, `Makefile` and the upstream CI and git files (`.github/`,
  `.travis.yml`, `.environmentLinux.yaml`, `.gitignore`) were removed.

This crate (`libbigwig-sys`) is a dev-dependency of `omicsweave-bbi` only. It
provides the source, its build (`build.rs`, without libcurl via `-DNOCURL`) and
raw bindings mirroring `bigWig.h` (`src/lib.rs`); it contains no test logic.
The differential tests and their wrapper over these bindings live in
`crates/bbi/tests/`. libBigWig is the library under pyBigWig, so matching it
keeps results identical to existing Python pipelines. It cannot write bigBed
files; bigBed tests generate their inputs in test code.
