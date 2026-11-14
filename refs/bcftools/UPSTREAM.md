# bcftools (reference implementation)

- Upstream: https://github.com/samtools/bcftools
- Version: 1.23.1-91-gd4254b11 (commit `d4254b11ec201fa531fe809a9c8c748e2707c515`)
- License: MIT/Expat (`vendor/bcftools/LICENSE`; the GPL option applies only
  to builds with the GNU Scientific Library, which is not used)
- Files: `consensus.c` and `filter.c` (the `-i`/`-e` expression engine)
  and what they need — `regidx.c/.h`, `smpl_ilist.c/.h`, `rbuf.h`,
  `bcftools.h`, `filter.h` — unmodified.

This crate (`bcftools-sys`) is a dev-dependency only; it defines the
behaviour `omicsweave-consensus` and the filter of `omicsweave-vcf` match
byte for byte. `build.rs` compiles the files against `htslib-sys`, writes
the `config.h` filter.c includes (empty: `ENABLE_PERL_FILTERS` undefined,
as bcftools' configure leaves it without `--enable-perl-filters`) and
renames bcftools' own `regidx_*`/`regitr_*` symbols (their ABI differs from
htslib's functions of the same names). `src/shim.c` replaces bcftools'
`error()` (which exits) with a return to the caller, carries
`apply_verbosity` and `parse_overlap_option` verbatim from `version.c`
(which also holds the exiting `error()`), and exposes
`ow_bcftools_consensus(argc, argv)`, `ow_filter_init(hdr, expr)` and
`ow_filter_test(filter, rec)` (the latter two return NULL / -1 where
bcftools reports an error). `src/lib.rs` declares these and the raw
`filter_init`/`filter_test`/`filter_destroy`; it contains no test logic.
