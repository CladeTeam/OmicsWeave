//! File-level alphabets: symbols, global codes, the byte lookup table and an
//! optional complement.

use crate::error::{Result, invalid};

/// The lookup-table value of a byte outside the alphabet. It stands for an
/// unknown input, not for any symbol, and never indexes a table.
pub const INVALID: u8 = 0xFF;

/// Symbols an alphabet may hold; `0xFF` is reserved for [`INVALID`].
pub const MAX_SYMBOLS: usize = 255;

/// An ordered set of single-byte symbols. Symbol `c` has global code `c`;
/// the lookup table maps every input byte to a global code or [`INVALID`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alphabet {
    name: String,
    symbols: Vec<u8>,
    complement: Option<Vec<u8>>,
    lut: [u8; 256],
}

impl Alphabet {
    /// An alphabet whose lookup table maps each symbol byte to its own code
    /// and every other byte to [`INVALID`].
    pub fn new(name: impl Into<String>, symbols: &[u8]) -> Result<Alphabet> {
        let mut lut = [INVALID; 256];
        for (code, &s) in symbols.iter().enumerate() {
            if code < MAX_SYMBOLS {
                lut[s as usize] = code as u8;
            }
        }
        Alphabet::from_parts(name.into(), symbols.to_vec(), None, lut)
    }

    /// Adds `byte` as another spelling of the symbol `symbol`, e.g. `a` for
    /// `A` or `U` for `T`.
    pub fn with_alias(mut self, byte: u8, symbol: u8) -> Result<Alphabet> {
        let Some(code) = self.code_of_symbol(symbol) else {
            return Err(invalid(
                "FFM alphabet",
                format!("alias target {:?} is not a symbol", symbol as char),
            ));
        };
        if self.symbols.contains(&byte) && byte != symbol {
            return Err(invalid(
                "FFM alphabet",
                format!("{:?} is already a symbol", byte as char),
            ));
        }
        self.lut[byte as usize] = code;
        Ok(self)
    }

    /// Maps lowercase ASCII letters to their uppercase symbols.
    pub fn with_lowercase(mut self) -> Result<Alphabet> {
        for &s in &self.symbols.clone() {
            if s.is_ascii_uppercase() && !self.symbols.contains(&s.to_ascii_lowercase()) {
                self = self.with_alias(s.to_ascii_lowercase(), s)?;
            }
        }
        Ok(self)
    }

    /// Sets the complement: `partners[c]` is the symbol paired with symbol
    /// `c`. The pairing must be an involution.
    pub fn with_complement(mut self, partners: &[u8]) -> Result<Alphabet> {
        if partners.len() != self.symbols.len() {
            return Err(invalid(
                "FFM alphabet",
                format!(
                    "complement has {} entries for {} symbols",
                    partners.len(),
                    self.symbols.len()
                ),
            ));
        }
        let mut complement = Vec::with_capacity(partners.len());
        for &p in partners {
            let Some(code) = self.code_of_symbol(p) else {
                return Err(invalid(
                    "FFM alphabet",
                    format!("complement {:?} is not a symbol", p as char),
                ));
            };
            complement.push(code);
        }
        self.complement = Some(complement);
        self.validate()?;
        Ok(self)
    }

    /// `ACGT`, lowercase accepted, with the Watson-Crick complement.
    pub fn dna() -> Alphabet {
        Alphabet::new("DNA", b"ACGT")
            .and_then(Alphabet::with_lowercase)
            .and_then(|a| a.with_complement(b"TGCA"))
            .expect("valid preset")
    }

    /// `ACGU`, lowercase accepted, with the Watson-Crick complement.
    pub fn rna() -> Alphabet {
        Alphabet::new("RNA", b"ACGU")
            .and_then(Alphabet::with_lowercase)
            .and_then(|a| a.with_complement(b"UGCA"))
            .expect("valid preset")
    }

    /// The twenty standard amino acids, lowercase accepted, no complement.
    pub fn protein() -> Alphabet {
        Alphabet::new("protein", b"ACDEFGHIKLMNPQRSTVWY")
            .and_then(Alphabet::with_lowercase)
            .expect("valid preset")
    }

    /// Builds an alphabet from its stored parts, checking every constraint.
    pub fn from_parts(
        name: String,
        symbols: Vec<u8>,
        complement: Option<Vec<u8>>,
        lut: [u8; 256],
    ) -> Result<Alphabet> {
        let alphabet = Alphabet {
            name,
            symbols,
            complement,
            lut,
        };
        alphabet.validate()?;
        Ok(alphabet)
    }

    fn validate(&self) -> Result<()> {
        let n = self.symbols.len();
        if n == 0 || n > MAX_SYMBOLS {
            return Err(invalid(
                "FFM alphabet",
                format!("{n} symbols, need 1..={MAX_SYMBOLS}"),
            ));
        }
        if self.name.len() > u16::MAX as usize {
            return Err(invalid("FFM alphabet", "name longer than 65535 bytes"));
        }
        let mut seen = [false; 256];
        for (code, &s) in self.symbols.iter().enumerate() {
            if std::mem::replace(&mut seen[s as usize], true) {
                return Err(invalid(
                    "FFM alphabet",
                    format!("symbol {:?} repeats", s as char),
                ));
            }
            if self.lut[s as usize] as usize != code {
                return Err(invalid(
                    "FFM alphabet",
                    format!("lookup of symbol {:?} is not its code {code}", s as char),
                ));
            }
        }
        if let Some(byte) = self
            .lut
            .iter()
            .position(|&c| c != INVALID && c as usize >= n)
        {
            return Err(invalid(
                "FFM alphabet",
                format!("lookup of byte {byte} is outside {n} codes"),
            ));
        }
        if let Some(comp) = &self.complement {
            if comp.len() != n {
                return Err(invalid("FFM alphabet", "complement length differs"));
            }
            for (c, &p) in comp.iter().enumerate() {
                if p as usize >= n || comp[p as usize] as usize != c {
                    return Err(invalid(
                        "FFM alphabet",
                        "complement is not an involution on the codes",
                    ));
                }
            }
        }
        Ok(())
    }

    fn code_of_symbol(&self, symbol: u8) -> Option<u8> {
        self.symbols
            .iter()
            .position(|&s| s == symbol)
            .map(|c| c as u8)
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Number of symbols, `n`.
    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    /// The canonical byte of each global code.
    pub fn symbols(&self) -> &[u8] {
        &self.symbols
    }

    /// The complement code of each global code, if the alphabet has one.
    pub fn complement(&self) -> Option<&[u8]> {
        self.complement.as_deref()
    }

    /// Byte → global code or [`INVALID`].
    pub fn lut(&self) -> &[u8; 256] {
        &self.lut
    }

    pub(crate) fn heap_bytes(&self) -> usize {
        self.name.capacity()
            + self.symbols.capacity()
            + self.complement.as_ref().map_or(0, Vec::capacity)
    }

    /// The global code of an input byte, `None` when it is invalid.
    pub fn code(&self, byte: u8) -> Option<u8> {
        let c = self.lut[byte as usize];
        (c != INVALID).then_some(c)
    }

    /// Maps `input` to global codes (or [`INVALID`]) in `out`.
    pub fn encode(&self, input: &[u8], out: &mut [u8]) -> Result<()> {
        omicsweave_core::check_output("encoded symbols", out.len(), input.len())?;
        for (o, &b) in out.iter_mut().zip(input) {
            // SAFETY: every u8 indexes the complete 256-entry byte table.
            *o = unsafe { *self.lut.get_unchecked(b as usize) };
        }
        Ok(())
    }

    pub(crate) fn write_to(&self, out: &mut Vec<u8>) {
        out.push(self.symbols.len() as u8);
        out.push(u8::from(self.complement.is_some()));
        out.extend_from_slice(&(self.name.len() as u16).to_le_bytes());
        out.extend_from_slice(self.name.as_bytes());
        out.extend_from_slice(&self.symbols);
        if let Some(comp) = &self.complement {
            out.extend_from_slice(comp);
        }
        out.extend_from_slice(&self.lut);
    }

    pub(crate) fn parse(payload: &[u8]) -> Result<Alphabet> {
        let mut r = crate::layout::Cursor::new(payload, "FFM alphabet");
        let n = r.u8()? as usize;
        let flags = r.u8()?;
        if flags & !1 != 0 {
            return Err(invalid(
                "FFM alphabet",
                format!("unknown flags {flags:#04x}"),
            ));
        }
        let name_len = r.u16()? as usize;
        let name = std::str::from_utf8(r.bytes(name_len)?)
            .map_err(|_| invalid("FFM alphabet", "name is not UTF-8"))?
            .to_owned();
        let symbols = r.bytes(n)?.to_vec();
        let complement = if flags & 1 != 0 {
            Some(r.bytes(n)?.to_vec())
        } else {
            None
        };
        let lut: [u8; 256] = r.bytes(256)?.try_into().expect("256 bytes");
        r.finish()?;
        Alphabet::from_parts(name, symbols, complement, lut)
    }
}
