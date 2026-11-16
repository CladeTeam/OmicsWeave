//! Compact diagnostics for potentially large encoded outputs.

#[track_caller]
pub fn bytes_eq(actual: &[u8], expected: &[u8]) {
    if actual == expected {
        return;
    }
    let at = actual
        .iter()
        .zip(expected)
        .position(|(a, b)| a != b)
        .unwrap_or(actual.len().min(expected.len()));
    panic!(
        "output differs at byte {at}: actual {:?} ({} bytes), expected {:?} ({} bytes)",
        actual.get(at),
        actual.len(),
        expected.get(at),
        expected.len()
    );
}
