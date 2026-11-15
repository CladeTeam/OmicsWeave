# htscodecs

- Upstream: https://github.com/samtools/htscodecs
- Version: 1.6.7 (release archive), the version bundled with htslib 1.25
- Archive: `https://github.com/samtools/htscodecs/releases/download/v1.6.7/htscodecs-1.6.7.tar.gz`,
  sha256 `5100e4b27646a27042a00b07cada35ce9fa3c7fc5aacd0de9874cd5d16c35fbf`
- License: BSD-3-Clause (`vendor/htscodecs/LICENSE.md`; `c_range_coder.h`,
  `rANS_byte.h` and `rANS_word.h` are public domain, as stated there)
- Local changes: sources unmodified. Copied from the archive: `LICENSE.md`,
  `README.md`, `NEWS.md` and every `.c` / `.h` file of `htscodecs/`. Left
  out: the autotools build system (`configure`, `Makefile.*`, `m4/`, ...)
  and `tests/`.

`build.rs` compiles the sources with `cc`, writing in `OUT_DIR` the
`config.h` that configure produces with a current compiler (SIMD rANS
decoders on x86-64 with their own `-msse4.1` / `-mavx2` / `-mavx512f`
flags and run-time CPU dispatch, NEON on AArch64, bzip2 from `bzip2-sys`)
and the `version.h` the upstream Makefile generates. `src/lib.rs` declares
the decoders `omicsweave-cram` uses; nothing else is wrapped here.
The native archive declares its Unix `libm` and pthread link dependencies
explicitly. Runtime dispatch and thread-local scratch pools remain the
unmodified upstream implementations.
