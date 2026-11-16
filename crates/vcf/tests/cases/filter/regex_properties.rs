//! Regular expressions of `~` against glibc's `regcomp`/`regexec` as
//! bcftools' filter.c uses them (`REG_NOSUB`, with `REG_ICASE` for `/i`):
//! random basic regular expressions, tried as `ID~"..."` on records with
//! random IDs, must compile, fail to compile and match alike.

use crate::filter_oracle;
use crate::support::variants;

use filter_oracle::{Tally, c_filter, compare, ours};
use proptest::prelude::*;
use proptest::test_runner::{Config, TestRunner};
use variants::Rng;

/// A random basic regular expression (no `"` and no back-references,
/// whose matching POSIX leaves to the implementation).
fn pattern(rng: &mut Rng) -> String {
    let mut s = String::new();
    if rng.chance(25) {
        s.push('^');
    }
    for _ in 0..rng.range(1, 5) {
        let piece = match rng.range(0, 12) {
            0..=3 => rng
                .pick(&["a", "b", "A", "x", "1", "_", "-", ".", ",", " "])
                .to_string(),
            4 => rng
                .pick(&[
                    "[ab]",
                    "[^a]",
                    "[a-c]",
                    "[]a]",
                    "[^]a]",
                    "[a-]",
                    "[-a]",
                    "[[:alpha:]]",
                    "[[:digit:]]",
                    "[[:upper:]]",
                    "[[:lower:]]",
                    "[[:space:]]",
                    "[[:punct:]]",
                    "[[.a.]]",
                    "[[=b=]]",
                    "[A-Z]",
                    "[0-9_]",
                    "[[:alnum:]_]",
                ])
                .to_string(),
            5 => rng
                .pick(&[
                    "\\w", "\\W", "\\s", "\\S", "\\b", "\\B", "\\<", "\\>", "\\`", "\\'",
                ])
                .to_string(),
            6 => rng
                .pick(&["\\.", "\\*", "\\[", "\\]", "\\a", "\\^", "\\$", "\\-"])
                .to_string(),
            7 => {
                let inner = rng
                    .pick(&["a", "ab", "a\\|b", "[0-9]", "x*", ""])
                    .to_string();
                format!("\\({inner}\\)")
            }
            8 => rng.pick(&["a\\|b", "x\\|", "\\|a"]).to_string(),
            9 => rng
                .pick(&[
                    "[",
                    "\\(",
                    "\\)",
                    "a\\{",
                    "\\{1\\}",
                    "[z-a]",
                    "[[:nope:]]",
                    "[[.ab.]]",
                    "\\",
                    "a**",
                    "a*\\{2\\}",
                    "\\1",
                    "\\(a\\)\\2",
                    "a\\{3,1\\}",
                    "a\\{,2\\}",
                    "a\\{99999\\}",
                    "[a",
                    "*a",
                    "^*",
                    "a\\{1,\\}",
                ])
                .to_string(),
            _ => "$".to_string(),
        };
        s.push_str(&piece);
        if rng.chance(30) {
            let rep: &&str = rng.pick(&[
                "*",
                "\\+",
                "\\?",
                "\\{2\\}",
                "\\{0,1\\}",
                "\\{1,\\}",
                "\\{0\\}",
            ]);
            s.push_str(rep);
        }
    }
    if rng.chance(15) {
        s.push('$');
    }
    if rng.chance(25) {
        s.push_str("/i");
    }
    s
}

fn subject(rng: &mut Rng) -> String {
    let n = rng.range(1, 8);
    (0..n)
        .map(|_| {
            *rng.pick(&[
                'a', 'b', 'A', 'B', 'x', 'X', '1', '9', '_', '-', '.', ',', '*', '[', ']', 'c',
            ])
        })
        .collect()
}

fn run(seed: u64, tally: &mut Tally) -> Result<(), TestCaseError> {
    let mut rng = Rng(seed);
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 0..10 {
        text += &format!("chr1\t{}\t{}\tA\tC\t.\t.\t.\n", i + 1, subject(&mut rng));
    }
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("t.vcf");
    std::fs::write(&vcf, &text).unwrap();
    for _ in 0..4 {
        let expr = format!("ID~\"{}\"", pattern(&mut rng));
        let c = c_filter(&vcf, &expr);
        let o = ours(text.as_bytes(), &expr);
        compare(&c, &o, tally)
            .map_err(|e| TestCaseError::fail(format!("{e}\nexpression: {expr}\n{text}")))?;
    }
    Ok(())
}

#[test]
fn regex_matches_glibc() {
    let mut runner = TestRunner::new(Config {
        cases: 400,
        source_file: Some(file!()),
        failure_persistence: Some(Box::new(
            proptest::test_runner::FileFailurePersistence::Direct(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/cases/filter/regex_properties.proptest-regressions"
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
    assert!(
        tally.records > 5000 && tally.passed > 500 && tally.init_errors > 100,
        "{tally:?}"
    );
}
