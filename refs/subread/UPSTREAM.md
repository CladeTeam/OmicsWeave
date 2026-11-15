# Subread featureCounts (reference implementation)

- Upstream: https://github.com/ShiLab-Bioinformatics/subread
- Version: 2.2.1 (`src/makefile.version`), commit
  `eb1cec17054f3acdc72631becaf39d581b40ae98`
- Archive: `https://codeload.github.com/ShiLab-Bioinformatics/subread/tar.gz/eb1cec17054f3acdc72631becaf39d581b40ae98`,
  sha256 `e2ba7a064af23359c11f5f9b9f126686398d5a95d5e8c4b4e85292bc17a4c8fc`
- License: GPL-3.0-or-later (the notice at the head of every source file).
  The crate is a test-only dev-dependency (`publish = false`) and is never
  linked into anything that is distributed.
- Files: `readSummary.c` (featureCounts) and the library files every Subread
  program links (`ALL_LIBS` of `src/Makefile.Linux`) with all headers of
  `src/`, unmodified.

This crate (`subread-sys`) defines the bulk gene counts `omicsweave-quant`
matches. `build.rs` compiles the files with the flags of `Makefile.Linux`
(`-DMAKE_FOR_EXON -DMAKE_STANDALONE`) plus `-fcommon`, renaming
featureCounts' `main` to `ow_featurecounts_main`, which `src/lib.rs`
declares, and Subread's global `find_path` (htslib, linked into the same
test binaries, has one too); it contains no test logic.
