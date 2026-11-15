//! A small JSON value for `summary.json`.

use std::fmt::Write as _;

/// A JSON value; objects keep their key order.
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i64),
    UInt(u64),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn obj<K: Into<String>>(fields: impl IntoIterator<Item = (K, Json)>) -> Json {
        Json::Obj(fields.into_iter().map(|(k, v)| (k.into(), v)).collect())
    }

    /// The value of `key` in an object.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(f) => f.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_u64(&self) -> Option<u64> {
        match *self {
            Json::UInt(v) => Some(v),
            Json::Int(v) if v >= 0 => Some(v as u64),
            _ => None,
        }
    }

    /// Indented text.
    pub fn pretty(&self) -> String {
        let mut s = String::new();
        self.write(&mut s, 0);
        s.push('\n');
        s
    }

    fn write(&self, s: &mut String, depth: usize) {
        let pad = |s: &mut String, d: usize| {
            for _ in 0..d {
                s.push_str("  ");
            }
        };
        match self {
            Json::Null => s.push_str("null"),
            Json::Bool(b) => s.push_str(if *b { "true" } else { "false" }),
            Json::Int(v) => {
                let _ = write!(s, "{v}");
            }
            Json::UInt(v) => {
                let _ = write!(s, "{v}");
            }
            Json::Num(v) if v.is_finite() => {
                let _ = write!(s, "{v}");
            }
            Json::Num(_) => s.push_str("null"),
            Json::Str(t) => quote(s, t),
            Json::Arr(v) if v.is_empty() => s.push_str("[]"),
            Json::Arr(v) => {
                s.push_str("[\n");
                for (i, x) in v.iter().enumerate() {
                    pad(s, depth + 1);
                    x.write(s, depth + 1);
                    if i + 1 < v.len() {
                        s.push(',');
                    }
                    s.push('\n');
                }
                pad(s, depth);
                s.push(']');
            }
            Json::Obj(f) if f.is_empty() => s.push_str("{}"),
            Json::Obj(f) => {
                s.push_str("{\n");
                for (i, (k, x)) in f.iter().enumerate() {
                    pad(s, depth + 1);
                    quote(s, k);
                    s.push_str(": ");
                    x.write(s, depth + 1);
                    if i + 1 < f.len() {
                        s.push(',');
                    }
                    s.push('\n');
                }
                pad(s, depth);
                s.push('}');
            }
        }
    }
}

fn quote(s: &mut String, t: &str) {
    s.push('"');
    for c in t.chars() {
        match c {
            '"' => s.push_str("\\\""),
            '\\' => s.push_str("\\\\"),
            '\n' => s.push_str("\\n"),
            '\r' => s.push_str("\\r"),
            '\t' => s.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(s, "\\u{:04x}", c as u32);
            }
            c => s.push(c),
        }
    }
    s.push('"');
}

impl From<u64> for Json {
    fn from(v: u64) -> Json {
        Json::UInt(v)
    }
}

impl From<usize> for Json {
    fn from(v: usize) -> Json {
        Json::UInt(v as u64)
    }
}

impl From<f64> for Json {
    fn from(v: f64) -> Json {
        Json::Num(v)
    }
}

impl From<bool> for Json {
    fn from(v: bool) -> Json {
        Json::Bool(v)
    }
}

impl From<&str> for Json {
    fn from(v: &str) -> Json {
        Json::Str(v.to_owned())
    }
}

impl From<String> for Json {
    fn from(v: String) -> Json {
        Json::Str(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_nested_values() {
        let j = Json::obj([
            ("a", Json::from(1u64)),
            ("b", Json::Arr(vec![Json::from("x\"y"), Json::Null])),
            ("c", Json::obj(Vec::<(String, Json)>::new())),
        ]);
        assert_eq!(
            j.pretty(),
            "{\n  \"a\": 1,\n  \"b\": [\n    \"x\\\"y\",\n    null\n  ],\n  \"c\": {}\n}\n"
        );
    }
}
