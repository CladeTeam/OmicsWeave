# htslib (reference implementation)

- Upstream: https://github.com/samtools/htslib
- Version: 1.25 (release archive)
- Archive: `https://github.com/samtools/htslib/releases/download/1.25/htslib-1.25.tar.bz2`,
  sha256 `158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730`
- License: MIT/Expat (`vendor/htslib/LICENSE`; bundled htscodecs:
  `vendor/htslib/htscodecs/LICENSE.md`)
- Local changes: sources unmodified. Removed from the archive: `test/`,
  `ref_cache/`, `samples/`, `m4/`, `htscodecs/tests/` and htscodecs docs, man
  pages, the configure/Makefile build system, the command-line programs
  (`bgzip.c`, `tabix.c`, `htsfile.c`, `annot-tsv.c`), and the plugin and
  remote-file backends (`plugin.c`, `hfile_libcurl.c`, `hfile_gcs.c`,
  `hfile_s3.c`) and the x86/ARM SIMD variants of the htscodecs rANS codec.

This crate (`htslib-sys`) is a dev-dependency only. `build.rs` compiles the
library without configure, generating `config.h`, `version.h` and
`config_vars.h` in `OUT_DIR` from the defaults htslib's own Makefile uses when
unconfigured, plus bzip2 and lzma (linked from `bzip2-sys` and `lzma-sys`, so
CRAM blocks of every method can be written and read) and without libcurl.
`src/lib.rs` declares raw bindings for the FASTA index, BGZF, VCF/BCF,
SAM/BAM/CRAM, format options and logging API;
`src/fields.c` reads `bcf1_t`/`bcf_hdr_t` fields whose layout is not stable,
and `src/bam.c` exposes `bam1_t` fields and the inline region iterator of
`sam.h`; none contains test logic. `samtools faidx`, `bgzip`, `samtools
index`, bcftools and pysam all run on this code, so it defines the behaviour
`omicsweave-fasta` matches for `.fa`, `.fai`, `.fa.gz` and `.gzi`,
`omicsweave-vcf` for VCF and BCF, `omicsweave-bam` for BAM, BAI and CSI, and
`omicsweave-cram` for CRAM and CRAI.
`build.rs` exports the header directories (`DEP_HTS_INCLUDE`,
`DEP_HTS_GENERATED`) for `bcftools-sys`.
