mod common;

use common::{Rng, htslib, ours};

#[test]
fn floating_literals_match_htslib() {
    let mut numbers: Vec<String> = [
        "0",
        "-0",
        "+0",
        "0.",
        ".0",
        "+.5",
        "  -1.5",
        "1.e2",
        ".5e-2",
        "2.2250738585072014e-308",
        "2.2250738585072012e-308",
        "4.9406564584124654e-324",
        "2.4703282292062327e-324",
        "1.7976931348623157e308",
        "1.7976931348623159e308",
        "1e999999999999999999999999",
        "-1e-999999999999999999999999",
        "nan",
        "-NAN(0x123456)",
        "nan(0777)",
        "nan(999999999999999999999999999999)",
        "nan(not_a_number)",
        "nan(0x)",
        "nan(09)",
        "nan()",
        "nan(123)",
        "infinity",
        "-INF",
        "1e",
        "0x",
        "0x.p1",
        "0x0p9999999999",
        "0x1.fffffffffffffp1023",
        "0x1.fffffffffffff8p1023",
        "0x1.fffffffffffff7p1023",
        "0x1p-1074",
        "0x1p-1075",
        "0x1.00000000000001p-1075",
        "0x1.fffffffffffffp-1023",
        "0x0.fffffffffffff8p-1022",
        "0x1.00000000000018p0",
        "abc",
        "+",
        "-",
        ".",
        "0x1p+",
        "0x1p-",
        "1junk",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect();
    let mut rng = Rng(0x11_5a_4d);
    for _ in 0..1024 {
        let sign = if rng.chance(50) { "-" } else { "" };
        let mantissa = rng.next();
        let exponent = rng.range(0, 2400) as i64 - 1200;
        numbers.push(format!("{sign}{mantissa}.{:019}e{exponent}", rng.next()));
        numbers.push(format!(
            "{sign}0x{mantissa:x}.{:016x}p{exponent}",
            rng.next()
        ));
    }
    let mut text = b"@SQ\tSN:chr1\tLN:1000\n".to_vec();
    for number in &numbers {
        text.extend_from_slice(
            format!(
                "r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tDD:d:{number}\tFF:f:{number}\tBB:B:f,{number}\n"
            )
            .as_bytes(),
        );
    }
    let expected = htslib(&text, true);
    let got = ours(&text);
    assert_eq!(got.header, expected.header);
    assert_eq!(got.records.len(), numbers.len());
    assert_eq!(expected.records.len(), numbers.len());
    for (i, number) in numbers.iter().enumerate() {
        assert_eq!(got.records[i], expected.records[i], "literal {number}");
    }
}

#[test]
fn byte_values_and_signed_char_aux_rules_match_htslib() {
    let mut text = b"@SQ\tSN:chr1\tLN:1000\n".to_vec();
    for value in 1..=255 {
        // LF ends a record before it reaches the encoder.
        if value == b'\n' {
            continue;
        }
        for prefix in [
            b"XX:A:".as_slice(),
            b"XX:i:1",
            b"XX:f:1",
            b"XX:Z:",
            b"XX:B:C,-",
        ] {
            text.extend_from_slice(b"r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\t");
            text.extend_from_slice(prefix);
            text.push(value);
            text.extend_from_slice(b"\tYY:i:2\n");
        }
    }
    let expected = htslib(&text, true);
    let got = ours(&text);
    assert_eq!(got.header, expected.header);
    assert_eq!(got.records, expected.records);
}
