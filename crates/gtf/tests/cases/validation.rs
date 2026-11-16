//! Recoverable hierarchy diagnostics and rejected annotation syntax.
use crate::support::{gtf_text, read};
use omicsweave_gtf::{Annotation, GtfLayout, ProblemKind, ReadOptions};

#[test]
fn problems_are_reported_or_raised() {
    let text = "\
c\ts\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";
c\ts\ttranscript\t1\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";
c\ts\texon\t1\t60\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";
c\ts\texon\t70\t80\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";
c\ts\texon\t70\t80\t.\t+\t.\ttranscript_id \"t9\";
c\ts\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";
c\ts\texon\t1\t10\t.\t+\t.\tgene_id \"g2\"; transcript_id \"t1\";
";
    let ann = read(text);
    let r = ann.report();
    for kind in [
        ProblemKind::OutsideParent,
        ProblemKind::StrandConflict,
        ProblemKind::MissingGeneId,
        ProblemKind::DuplicateId,
        ProblemKind::GeneConflict,
    ] {
        assert!(r.count(kind) > 0, "{kind:?}: {:?}", r.problems);
    }
    // The orphan and the duplicate gene line are kept as others.
    assert_eq!(ann.others().len(), 2);
    assert_eq!(gtf_text(&ann, GtfLayout::Original), text);
    let strict = ReadOptions {
        strict: true,
        ..ReadOptions::default()
    };
    let err = Annotation::parse(text.as_bytes(), &strict).unwrap_err();
    assert!(
        matches!(err, omicsweave_core::Error::Format { line: 3, .. }),
        "{err}"
    );

    let gff = "##gff-version 3\nc\ts\tgene\t1\t100\t.\t+\t.\tID=g\nc\ts\texon\t1\t10\t.\t+\t.\tParent=nope\nc\ts\texon\t1\t10\t.\t+\t.\tID=e\nc\ts\tmRNA\t1\t100\t.\t+\t.\tID=g\n";
    let ann = read(gff);
    let r = ann.report();
    assert_eq!(r.count(ProblemKind::UnknownParent), 1);
    assert_eq!(r.count(ProblemKind::MissingParent), 1);
    assert_eq!(r.count(ProblemKind::DuplicateId), 1);
    assert!(Annotation::parse(gff.as_bytes(), &strict).is_err());
}

#[test]
fn malformed_records_are_rejected_in_report_and_strict_modes() {
    for bad in [
        "c\ts\texon\t0\t10\t.\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t10\t9\t.\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\t.\tx\t.\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\t.\t+\t3\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\tabc\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t1\t9\n",
        "c\ts\texon\t1\t9\t.\t+\t.\tgene_id \"g;\n",
        "c\ts\texon\t1\t4294967296\t.\t+\t.\tgene_id \"g\";\n",
        "c\ts\texon\t-1\t9\t.\t+\t.\tgene_id \"g\";\n",
    ] {
        let text = format!("# annotation\n\n{bad}");
        for strict in [false, true] {
            for (threads, batch_bytes) in [(1, 1), (4, 32)] {
                let opts = ReadOptions {
                    strict,
                    threads,
                    batch_bytes,
                    ..ReadOptions::default()
                };
                let err = Annotation::parse(text.as_bytes(), &opts).unwrap_err();
                assert!(
                    matches!(err, omicsweave_core::Error::Format { line: 3, .. }),
                    "input={text:?}, strict={strict}, threads={threads}, batch={batch_bytes}: {err}"
                );
            }
        }
    }
}

#[test]
fn isolated_hierarchy_problems_keep_their_line_and_kind_across_batches() {
    let gene = "c\t.\tgene\t1\t20\t.\t+\t.\tgene_id \"g\";\n";
    for (row, kind) in [
        (
            "c\t.\ttranscript\t1\t21\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n",
            ProblemKind::OutsideParent,
        ),
        (
            "c\t.\ttranscript\t1\t20\t.\t-\t.\tgene_id \"g\"; transcript_id \"t\";\n",
            ProblemKind::StrandConflict,
        ),
        (
            "c\t.\ttranscript\t1\t20\t.\t+\t.\tgene_id \"g\";\n",
            ProblemKind::MissingTranscriptId,
        ),
    ] {
        let text = format!("# annotation\n{gene}{row}");
        for (threads, batch_bytes) in [(1, 1), (4, 32), (2, 1 << 20)] {
            let opts = ReadOptions {
                threads,
                batch_bytes,
                ..ReadOptions::default()
            };
            let ann = Annotation::parse(text.as_bytes(), &opts).unwrap();
            assert_eq!(
                ann.report().count(kind),
                1,
                "kind={kind:?}, input={text:?}, options={opts:?}"
            );
            assert_eq!(ann.report().problems[0].line, 3, "{:?}", ann.report());
            let err = Annotation::parse(
                text.as_bytes(),
                &ReadOptions {
                    strict: true,
                    ..opts
                },
            )
            .unwrap_err();
            match err {
                omicsweave_core::Error::Format { line, detail } => {
                    assert_eq!(line, 3, "input={text:?}");
                    assert!(
                        detail.contains(kind.describe()),
                        "kind={kind:?}, detail={detail}"
                    );
                }
                other => panic!("expected annotation problem for {kind:?}: {other}"),
            }
        }
    }
}

#[test]
fn gff3_parent_cycles_are_diagnosed_without_recursion() {
    for body in [
        "c\t.\tregion\t1\t10\t.\t+\t.\tID=a;Parent=a\n",
        "c\t.\tregion\t1\t10\t.\t+\t.\tID=a;Parent=b\nc\t.\tregion\t1\t10\t.\t+\t.\tID=b;Parent=a\n",
        "c\t.\tregion\t1\t10\t.\t+\t.\tID=root\nc\t.\tregion\t1\t10\t.\t+\t.\tID=a;Parent=root,b\nc\t.\tregion\t1\t10\t.\t+\t.\tID=b;Parent=a\n",
    ] {
        let text = format!("##gff-version 3\n{body}");
        let mut expected = None;
        for (threads, batch_bytes) in [(1, 1), (4, 32), (2, 1 << 20)] {
            let opts = ReadOptions {
                threads,
                batch_bytes,
                ..ReadOptions::default()
            };
            let ann = Annotation::parse(text.as_bytes(), &opts).unwrap();
            assert!(ann.report().count(ProblemKind::BadNesting) > 0);
            let strict = ReadOptions {
                strict: true,
                ..opts
            };
            let err = Annotation::parse(text.as_bytes(), &strict)
                .unwrap_err()
                .to_string();
            assert_eq!(expected.get_or_insert_with(|| err.clone()), &err);
        }
    }
    // A long parent chain must neither exhaust the stack nor look cyclic.
    let mut text = "##gff-version 3\n".to_owned();
    for i in 0..16_384 {
        use std::fmt::Write;
        writeln!(
            text,
            "c\t.\tregion\t1\t10\t.\t+\t.\tID=r{i};Parent=r{}",
            i + 1
        )
        .unwrap();
    }
    text.push_str("c\t.\tregion\t1\t10\t.\t+\t.\tID=r16384\n");
    let ann = read(&text);
    assert_eq!(ann.report().count(ProblemKind::BadNesting), 0);

    let diamond = b"##gff-version 3\nc\t.\tregion\t1\t10\t.\t+\t.\tID=leaf;Parent=left,right\nc\t.\tregion\t1\t10\t.\t+\t.\tID=left;Parent=root\nc\t.\tregion\t1\t10\t.\t+\t.\tID=right;Parent=root\nc\t.\tregion\t1\t10\t.\t+\t.\tID=root\n";
    let ann = Annotation::parse(diamond, &ReadOptions::default()).unwrap();
    assert_eq!(ann.report().count(ProblemKind::BadNesting), 0);
}
