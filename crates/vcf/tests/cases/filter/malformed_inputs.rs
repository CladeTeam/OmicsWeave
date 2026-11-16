//! Broad malformed inputs exercise error recovery and independent Filter clones.
//! C's observations after uninitialized reads are not expected boolean values.

use crate::support::{rich_variants, variants::Rng};
use omicsweave_core::Error;
use omicsweave_vcf::{Filter, Record, VariantReader};
use proptest::prelude::*;
use proptest::test_runner::{Config, TestCaseError, TestRunner};

#[derive(Debug, Default)]
struct Tally {
    compile_errors: usize,
    records: usize,
    rejected: usize,
    runtime_errors: usize,
}

fn check(bytes: &[u8], expression: &str, tally: &mut Tally) -> Result<(), TestCaseError> {
    let mut reader = VariantReader::new(bytes).unwrap();
    let mut filter = match Filter::new(reader.header(), expression) {
        Ok(filter) => filter,
        Err(Error::Invalid { .. } | Error::Unsupported(_)) => {
            tally.compile_errors += 1;
            return Ok(());
        }
        Err(error) => {
            return Err(TestCaseError::fail(format!(
                "unexpected compile error: {error:?}"
            )));
        }
    };
    let mut copy = filter.clone();
    let mut record = Record::default();
    let mut row = 0;
    loop {
        match reader.read(&mut record) {
            Ok(false) => break,
            Err(_) => {
                prop_assert_eq!(&record, &Record::default());
                tally.rejected += 1;
            }
            Ok(true) => {
                let left = filter.test(reader.header(), &record);
                let right = copy.test(reader.header(), &record);
                prop_assert_eq!(
                    format!("{left:?}"),
                    format!("{right:?}"),
                    "expression={}, row={}",
                    expression,
                    row
                );
                match left {
                    Ok(_) => tally.records += 1,
                    Err(Error::Invalid { .. } | Error::Unsupported(_)) => tally.runtime_errors += 1,
                    Err(error) => {
                        return Err(TestCaseError::fail(format!(
                            "unexpected evaluation error: {error:?}"
                        )));
                    }
                }
                row += 1;
                if row % 3 == 0 {
                    // A snapshot taken after failures must own independent
                    // mutable scratch, just like the original compile clone.
                    copy = filter.clone();
                }
            }
        }
    }
    Ok(())
}

fn run(seed: u64, tally: &mut Tally) -> Result<(), TestCaseError> {
    let mut rng = Rng(seed);
    let (text, spec) = rich_variants::random_vcf(&mut rng);
    let directory = tempfile::tempdir().unwrap();
    let vcf = directory.path().join("malformed.vcf");
    let bcf = directory.path().join("malformed.bcf");
    std::fs::write(&vcf, &text).unwrap();
    rich_variants::htslib_to_bcf(&vcf, &bcf);
    let bcf = std::fs::read(&bcf).unwrap();
    for expression in crate::filter_expressions::generate(&mut rng, &spec, directory.path()) {
        for bytes in [text.as_bytes(), &bcf] {
            check(bytes, &expression, tally).map_err(|error| {
                TestCaseError::fail(format!(
                    "seed={seed}, expression={expression}: {error}\n{text}"
                ))
            })?;
        }
    }
    Ok(())
}

#[test]
fn malformed_records_and_expressions_preserve_independent_filter_states() {
    let mut runner = TestRunner::new(Config {
        cases: 400,
        source_file: Some(file!()),
        failure_persistence: Some(Box::new(
            proptest::test_runner::FileFailurePersistence::Direct(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/cases/filter/malformed_inputs.proptest-regressions"
            )),
        )),
        ..Config::default()
    });
    let tally = std::cell::RefCell::new(Tally::default());
    runner
        .run(&any::<u64>(), |seed| run(seed, &mut tally.borrow_mut()))
        .unwrap();
    let tally = tally.into_inner();
    assert!(
        tally.records > 1000 && tally.compile_errors > 50 && tally.rejected > 50,
        "{tally:?}"
    );
}

#[test]
fn initialized_buffer_boundaries_preserve_error_recovery_and_clone_state() {
    let mut tally = Tally::default();
    for seed in [4_796_652_997_821_466_650, 18_387_510_760_104_561_416] {
        run(seed, &mut tally).unwrap();
    }
    assert!(tally.runtime_errors > 0, "{tally:?}");
}
