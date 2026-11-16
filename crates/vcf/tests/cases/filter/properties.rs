//! Filter expressions against bcftools' filter.c: random expressions over
//! the tags of type-correct VCFs with initialized string layouts (and the
//! BCF htslib writes from them) must
//! compile, fail to compile, pass, fail and report errors exactly where
//! `filter_init` and `filter_test(filter, rec, NULL)` do. Outcomes filter.c
//! leaves undefined (process exits, assertion failures, crashes) end a
//! comparison without deciding it.
//! Implementation failures against a defined oracle always fail.

use crate::filter_oracle;
use crate::support::rich_variants;
use crate::support::variants;

use std::path::Path;

use filter_oracle::{Tally, c_filter, compare, ours};
use proptest::prelude::*;
use proptest::test_runner::{Config, TestRunner};
use variants::Rng;

fn check(
    path: &Path,
    bytes: &[u8],
    expr: &str,
    what: &str,
    text: &str,
    tally: &mut Tally,
    compile_error: bool,
) -> Result<(), TestCaseError> {
    let c = c_filter(path, expr);
    match &c {
        filter_oracle::Run::InitError if compile_error => {}
        filter_oracle::Run::Records(rows)
            if !compile_error
                && rows.len() == 24
                && rows
                    .iter()
                    .all(|row| matches!(row, filter_oracle::Outcome::Pass(_))) => {}
        _ => {
            return Err(TestCaseError::fail(format!(
                "typed oracle must finish every valid record (compile_error={compile_error}): {c:?}\nexpression: {expr}\n{text}"
            )));
        }
    }
    let o = ours(bytes, expr);
    compare(&c, &o, tally)
        .map(|_| ())
        .map_err(|e| TestCaseError::fail(format!("{what}: {e}\nexpression: {expr}\n{text}")))
}

fn run(seed: u64, tally: &mut Tally) -> Result<(), TestCaseError> {
    let mut rng = Rng(seed);
    let (text, _) = crate::filter_inputs::defined_vcf(&mut rng);
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("t.vcf");
    std::fs::write(&vcf, &text).unwrap();
    let bcf = dir.path().join("t.bcf");
    let kept = rich_variants::htslib_to_bcf(&vcf, &bcf);
    prop_assert!(
        kept.len() == 24 && kept.iter().all(|accepted| *accepted),
        "typed corpus must be accepted completely: {:?}\n{}",
        kept,
        text
    );
    let kept_text = rich_variants::kept_lines(&text, &kept);
    let bcf_bytes = std::fs::read(&bcf).unwrap();

    let expressions = crate::filter_inputs::defined_expressions(&mut rng, dir.path());
    for (i, expr) in expressions.iter().enumerate() {
        let compile_error = i == expressions.len() - 1;
        check(
            &vcf,
            text.as_bytes(),
            expr,
            "vcf",
            &text,
            tally,
            compile_error,
        )?;
        check(
            &bcf,
            &bcf_bytes,
            expr,
            "bcf",
            &kept_text,
            tally,
            compile_error,
        )?;
    }
    Ok(())
}

#[test]
fn filter_matches_bcftools() {
    let mut runner = TestRunner::new(Config {
        cases: 400,
        source_file: Some(file!()),
        failure_persistence: Some(Box::new(
            proptest::test_runner::FileFailurePersistence::Direct(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/cases/filter/properties.proptest-regressions"
            )),
        )),
        ..Config::default()
    });
    let tally = std::cell::RefCell::new(Tally::default());
    let result = runner.run(&any::<u64>(), |seed| {
        run(seed, &mut tally.borrow_mut())
            .map_err(|error| TestCaseError::fail(format!("seed={seed}: {error}")))
    });
    let tally = tally.into_inner();
    eprintln!("{tally:?}");
    result.unwrap();
    // The comparisons must not be vacuous.
    assert!(tally.records > 2000, "{tally:?}");
    assert!(tally.passed > 200, "{tally:?}");
    assert!(tally.init_errors > 50, "{tally:?}");
}
