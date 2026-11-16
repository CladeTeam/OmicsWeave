use super::*;

/// Table lines `(tid, start, end, context, strand, meth, unmeth)`, the
/// CpGs of both strands merged with `merge`.
pub(super) type Line = (usize, u32, u32, usize, char, u32, u32);

pub(super) fn lines(c: &BTreeMap<(usize, u32, u8), (u32, u32)>, merge: bool) -> Vec<Line> {
    // Lines at one position come forward strand first, then by context.
    type Key = (usize, u32, u32, bool, usize, char);
    let mut m: BTreeMap<Key, (u32, u32)> = BTreeMap::new();
    for (&(t, p, code), &(a, b)) in c {
        let ctx = (code & 3) as usize;
        let rev = code & 4 != 0;
        let key = if merge && ctx == 0 {
            let s = if rev { p - 1 } else { p };
            (t, s, s + 2, false, 0, '.')
        } else {
            (t, p, p + 1, rev, ctx, if rev { '-' } else { '+' })
        };
        let e = m.entry(key).or_default();
        e.0 += a;
        e.1 += b;
    }
    m.into_iter()
        .map(|((t, s, e, _, c, st), (a, b))| (t, s, e, c, st, a, b))
        .collect()
}

pub(super) fn text(f: &Fixture, lines: &[Line]) -> String {
    let mut s = String::new();
    for &(t, a, b, c, st, m, u) in lines {
        s.push_str(&format!(
            "{}\t{a}\t{b}\t{}\t{st}\t{m}\t{u}\n",
            f.refs[t].0, CONTEXTS[c]
        ));
    }
    s
}

/// Equal texts, or a panic naming the first line that differs.
pub(super) fn same_text(got: &str, want: &str) {
    if got != want {
        let (g, w): (Vec<&str>, Vec<&str>) = (got.lines().collect(), want.lines().collect());
        let i = g
            .iter()
            .zip(&w)
            .position(|(a, b)| a != b)
            .unwrap_or(g.len().min(w.len()));
        panic!(
            "{} vs {} lines; first difference at {i}: got {:?}, want {:?}",
            g.len(),
            w.len(),
            g.get(i.saturating_sub(1)..(i + 2).min(g.len())),
            w.get(i.saturating_sub(1)..(i + 2).min(w.len()))
        );
    }
}
