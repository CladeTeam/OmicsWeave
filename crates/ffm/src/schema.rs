//! Annotation keys: declared once per file; each record stores only values.

use crate::error::{Result, invalid};
use crate::layout::Cursor;
use std::collections::HashSet;

/// How a key's non-empty values must parse.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum KeyType {
    String,
    /// Decimal `i64`.
    Integer,
    /// Decimal `f64`.
    Float,
}

/// A non-missing annotation interpreted according to its schema.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AnnotationValue<'a> {
    String(&'a str),
    Integer(i64),
    Float(f64),
}

impl KeyType {
    fn id(self) -> u8 {
        match self {
            KeyType::String => 0,
            KeyType::Integer => 1,
            KeyType::Float => 2,
        }
    }

    fn from_id(id: u8) -> Option<KeyType> {
        match id {
            0 => Some(KeyType::String),
            1 => Some(KeyType::Integer),
            2 => Some(KeyType::Float),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Key {
    pub name: String,
    pub kind: KeyType,
}

impl Key {
    pub fn new(name: impl Into<String>, kind: KeyType) -> Key {
        Key {
            name: name.into(),
            kind,
        }
    }
}

/// The ordered annotation keys of a file. Motif IDs are not keys.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Schema {
    keys: Vec<Key>,
}

/// Bytes no stored name or value may contain.
pub(crate) fn forbidden(b: u8) -> bool {
    matches!(b, b'\t' | b'\n' | b'\r' | 0)
}

impl Schema {
    pub fn new(keys: Vec<Key>) -> Result<Schema> {
        if keys.len() > u16::MAX as usize {
            return Err(invalid("FFM schema", "more than 65535 keys"));
        }
        let mut names = HashSet::with_capacity(keys.len());
        for (i, key) in keys.iter().enumerate() {
            let name = key.name.as_bytes();
            if name.is_empty() || name.len() > u16::MAX as usize {
                return Err(invalid("FFM schema", format!("key {i} name length")));
            }
            if name.iter().copied().any(forbidden) {
                return Err(invalid(
                    "FFM schema",
                    format!("key {:?} contains a tab, newline or NUL", key.name),
                ));
            }
            if !names.insert(key.name.as_str()) {
                return Err(invalid("FFM schema", format!("key {:?} repeats", key.name)));
            }
        }
        Ok(Schema { keys })
    }

    pub fn keys(&self) -> &[Key] {
        &self.keys
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// The position of a key by name.
    pub fn index(&self, name: &str) -> Option<usize> {
        self.keys.iter().position(|k| k.name == name)
    }

    pub(crate) fn value<'a>(&self, i: usize, text: &'a str) -> Option<AnnotationValue<'a>> {
        if text.is_empty() {
            return None;
        }
        Some(match self.keys.get(i)?.kind {
            KeyType::String => AnnotationValue::String(text),
            KeyType::Integer => AnnotationValue::Integer(text.parse().expect("validated integer")),
            KeyType::Float => AnnotationValue::Float(text.parse().expect("validated float")),
        })
    }

    pub(crate) fn heap_bytes(&self) -> usize {
        self.keys.capacity() * std::mem::size_of::<Key>()
            + self.keys.iter().map(|k| k.name.capacity()).sum::<usize>()
    }

    /// Checks one value of key `i`; an empty value is a missing one.
    pub(crate) fn check_value(&self, i: usize, value: &str) -> Result<()> {
        let key = &self.keys[i];
        if value.bytes().any(|b| b == b'\t' || b == 0) {
            return Err(invalid(
                "FFM annotation",
                format!("value of {:?} contains a tab or NUL", key.name),
            ));
        }
        if value.is_empty() {
            return Ok(());
        }
        let ok = match key.kind {
            KeyType::String => true,
            KeyType::Integer => value.parse::<i64>().is_ok(),
            KeyType::Float => value.parse::<f64>().is_ok(),
        };
        if !ok {
            return Err(invalid(
                "FFM annotation",
                format!("value {value:?} of {:?} is not a {:?}", key.name, key.kind),
            ));
        }
        Ok(())
    }

    pub(crate) fn write_to(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&(self.keys.len() as u16).to_le_bytes());
        for key in &self.keys {
            out.push(key.kind.id());
            out.extend_from_slice(&(key.name.len() as u16).to_le_bytes());
            out.extend_from_slice(key.name.as_bytes());
        }
    }

    pub(crate) fn parse(payload: &[u8]) -> Result<Schema> {
        let mut r = Cursor::new(payload, "FFM schema");
        let count = r.u16()? as usize;
        let mut keys = Vec::with_capacity(count);
        for _ in 0..count {
            let id = r.u8()?;
            let kind = KeyType::from_id(id)
                .ok_or_else(|| invalid("FFM schema", format!("unknown key type {id}")))?;
            let len = r.u16()? as usize;
            let name = std::str::from_utf8(r.bytes(len)?)
                .map_err(|_| invalid("FFM schema", "key name is not UTF-8"))?;
            keys.push(Key::new(name, kind));
        }
        r.finish()?;
        Schema::new(keys)
    }
}
