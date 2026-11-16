//! Oracle incomparability must never conceal a defined implementation failure.

use crate::filter_oracle::{Comparison, Outcome, Run, Tally, compare};

#[test]
fn unwritten_string_errors_do_not_skip_a_defined_oracle() {
    let oracle = Run::Records(vec![Outcome::Pass(true)]);
    let implementation = Run::Records(vec![Outcome::ImplementationFailure(
        "Unsupported(\"read of a filter string byte never written\")".into(),
    )]);
    assert!(compare(&oracle, &implementation, &mut Tally::default()).is_err());
}

#[test]
fn defined_oracle_rejects_unsupported_or_missing_implementation_outcomes() {
    let expected = Run::Records(vec![Outcome::Pass(true), Outcome::Pass(false)]);
    for implementation in [
        Run::InitUndefined,
        Run::InitError,
        Run::InitFailure("unsupported compile".into()),
        Run::Records(vec![Outcome::ImplementationFailure(
            "unsupported evaluation".into(),
        )]),
        Run::Records(vec![Outcome::Undefined]),
        Run::Records(vec![Outcome::Pass(true), Outcome::Undefined]),
        Run::Records(vec![Outcome::Pass(true)]),
        Run::Records(vec![Outcome::Error]),
    ] {
        let mut tally = Tally::default();
        assert!(
            compare(&expected, &implementation, &mut tally).is_err(),
            "{implementation:?}"
        );
        assert_eq!(
            tally.oracle_init_undefined + tally.oracle_runtime_undefined,
            0
        );
    }
    assert!(compare(&Run::InitError, &Run::InitUndefined, &mut Tally::default()).is_err());
}

#[test]
fn only_oracle_undefined_is_incomparable_and_defined_prefix_still_must_match() {
    let expected = Run::Records(vec![Outcome::Pass(true), Outcome::Undefined]);
    let actual = Run::Records(vec![Outcome::Pass(true)]);
    let mut tally = Tally::default();
    assert_eq!(
        compare(&expected, &actual, &mut tally).unwrap(),
        Comparison::OracleUndefined
    );
    assert_eq!((tally.records, tally.oracle_runtime_undefined), (1, 1));
    assert!(
        compare(
            &expected,
            &Run::Records(vec![Outcome::Pass(false)]),
            &mut Tally::default()
        )
        .is_err()
    );
    assert_eq!(
        compare(&Run::InitUndefined, &actual, &mut tally).unwrap(),
        Comparison::OracleUndefined
    );
    assert_eq!(tally.oracle_init_undefined, 1);
    assert_eq!(
        compare(&actual, &actual, &mut tally).unwrap(),
        Comparison::Compared
    );
}
