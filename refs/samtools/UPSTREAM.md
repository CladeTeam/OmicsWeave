# samtools (reference implementation)

- Upstream: https://github.com/samtools/samtools
- Version: 1.25 (release archive), built against the vendored htslib 1.25
- Archive: `https://github.com/samtools/samtools/releases/download/1.25/samtools-1.25.tar.bz2`,
  sha256 `8bbd2b3cb447953149a74ab0b00e3ec4e88a789a27b8373a31fa6b5b2dee02b1`
- License: MIT/Expat (`vendor/samtools/LICENSE`)
- Files: `bam_sort.c` (`samtools sort`), `bam2depth.c` (`samtools depth`)
  and what they need — `bam.c/.h`, `bedidx.c/.h`, `sam_opts.c/.h`,
  `sam_utils.c/.h`, `samtools.h` — unmodified. `bam2depth.c`'s `-f` file
  list reader lives in `bam_plcmd.c` (mpileup), which is not vendored; the
  shim's stand-in rejects the option.

This crate (`samtools-sys`) is a dev-dependency only; `samtools sort` defines
the record order, header and BGZF framing the BAM sort of
`omicsweave-convert` matches, and `samtools depth` the per-base coverage of
`omicsweave-quant`. `build.rs` compiles the files against
`htslib-sys` with an empty `config.h` (no configure). `src/shim.c` carries
`samtools_version()` from `bamtk.c` (which also holds `main`) and exposes
`ow_samtools(argc, argv)`, which runs a command by name as samtools' `main`
does; `src/lib.rs` declares it and contains no test logic. Another command
is added by vendoring its file (and what it needs), listing it in
`build.rs` and adding a row to the command table of `src/shim.c`.
