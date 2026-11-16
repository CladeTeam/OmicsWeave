//! The records of a region as columns: QUAL, FILTER and chosen INFO and
//! FORMAT fields beside the [`GenotypeMatrix`], optionally only those
//! passing a bcftools filter expression ([`BcfReader::records`]).
//!
//! Numeric fields are laid out densely at the widest record of the region
//! (and, for FORMAT, the widest sample), padded with htslib's vector-end
//! sentinel; a record without the field holds only padding. Missing values
//! keep htslib's missing sentinel.
//!
//! [`BcfReader::records`]: crate::BcfReader::records

use omicsweave_core::{Error, Result};

use crate::fields::{FLOAT_VECTOR_END, GetError, INT32_VECTOR_END, Span, format_view, info_view};
use crate::filter::Filter;
use crate::header::{Header, ValueType};
use crate::indexed::{GenotypeMatrix, Strings};
use crate::record::Record;

/// Which fields a records query decodes besides the genotypes.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RecordQuery {
    /// INFO keys, in the order their columns come.
    pub info: Vec<String>,
    /// FORMAT keys (other than `GT`, which is the genotype matrix).
    pub format: Vec<String>,
}

/// Records kept by a bcftools expression (`-i`) or left out by one (`-e`).
pub struct RecordFilter {
    filter: Filter,
    include: bool,
}

impl RecordFilter {
    pub fn include(header: &Header, expr: &str) -> Result<RecordFilter> {
        Ok(RecordFilter {
            filter: Filter::new(header, expr)?,
            include: true,
        })
    }

    pub fn exclude(header: &Header, expr: &str) -> Result<RecordFilter> {
        Ok(RecordFilter {
            filter: Filter::new(header, expr)?,
            include: false,
        })
    }

    pub fn keeps(&mut self, header: &Header, rec: &Record) -> Result<bool> {
        Ok(self.filter.test(header, rec)? == self.include)
    }
}

/// The values of one field over a region's records.
#[derive(Debug, Clone, PartialEq)]
pub enum Values {
    /// An INFO flag per record.
    Flag(Vec<bool>),
    /// `[record][(sample)][width]`, missing values [`crate::INT32_MISSING`],
    /// padding [`INT32_VECTOR_END`].
    Int { width: usize, values: Vec<i32> },
    /// As `Int`, with the float sentinels' bit patterns
    /// ([`crate::FLOAT_MISSING`], [`FLOAT_VECTOR_END`]).
    Float { width: usize, values: Vec<f32> },
    /// A string per record (INFO) or per record and sample (FORMAT), as
    /// stored; empty when absent.
    Str(Strings),
}

/// A field's name and values.
#[derive(Debug, Clone, PartialEq)]
pub struct Column {
    pub name: String,
    pub values: Values,
}

/// The records of a region.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantTable {
    /// Positions, ID, alleles and the genotypes of the chosen samples.
    pub genotypes: GenotypeMatrix,
    /// QUAL; the missing bit pattern ([`crate::FLOAT_MISSING`]) when `.`.
    pub qual: Vec<f32>,
    /// FILTER names joined by `;`; empty for `.`.
    pub filters: Strings,
    pub info: Vec<Column>,
    pub format: Vec<Column>,
}

/// One field being collected: per record, its values and how many.
enum Collect {
    Flag(Vec<bool>),
    Int(Vec<i32>, Vec<usize>),
    Float(Vec<f32>, Vec<usize>),
    Str(Strings),
}

impl Collect {
    fn of(ty: ValueType) -> Collect {
        match ty {
            ValueType::Flag => Collect::Flag(Vec::new()),
            ValueType::Integer => Collect::Int(Vec::new(), Vec::new()),
            ValueType::Float => Collect::Float(Vec::new(), Vec::new()),
            ValueType::String => Collect::Str(Strings::default()),
        }
    }
}

struct Field {
    name: String,
    key: usize,
    collect: Collect,
}

/// Resolves requested dictionary IDs in one pass over a record's fields.
/// Only requested slots need resetting; no per-record header-wide clearing.
fn positions(slots: &mut [usize], requested: &[Field], fields: &[Span]) {
    const WANTED: usize = usize::MAX - 1;
    let mut left = 0usize;
    for f in requested {
        if slots[f.key] != WANTED {
            slots[f.key] = WANTED;
            left += 1;
        }
    }
    // Keep the first occurrence, and stop once every requested key is found.
    for (i, f) in fields.iter().enumerate() {
        if let Some(slot) = slots.get_mut(f.key as usize)
            && *slot == WANTED
        {
            *slot = i;
            left -= 1;
            if left == 0 {
                return;
            }
        }
    }
    // No WANTED markers may survive into the next INFO or FORMAT lookup.
    for f in requested {
        if slots[f.key] == WANTED {
            slots[f.key] = usize::MAX;
        }
    }
}

/// `values` of `counts[i]` each (times `groups` per record) laid out at
/// the widest, padded with `end`. Uniform widths already have this layout.
fn dense<T: Copy>(
    values: Vec<T>,
    counts: &[usize],
    groups: usize,
    end: T,
) -> Result<(usize, Vec<T>)> {
    let width = counts.iter().copied().max().unwrap_or(0);
    let len = counts
        .len()
        .checked_mul(groups)
        .and_then(|n| n.checked_mul(width))
        .ok_or_else(|| Error::invalid("variants", "column length overflows usize"))?;
    if counts.iter().all(|&n| n == width) {
        debug_assert_eq!(values.len(), len);
        return Ok((width, values));
    }
    let mut out = Vec::new();
    out.try_reserve_exact(len)
        .map_err(|e| Error::invalid("variants", e.to_string()))?;
    let mut source = values.as_slice();
    for &n in counts {
        for _ in 0..groups {
            let (row, rest) = source.split_at(n);
            out.extend_from_slice(row);
            out.resize(out.len() + width - n, end);
            source = rest;
        }
    }
    debug_assert!(source.is_empty());
    Ok((width, out))
}

fn finish(f: Field, groups: usize) -> Result<Column> {
    let values = match f.collect {
        Collect::Flag(v) => Values::Flag(v),
        Collect::Int(v, n) => {
            let (width, values) = dense(v, &n, groups, INT32_VECTOR_END)?;
            Values::Int { width, values }
        }
        Collect::Float(v, n) => {
            let (width, values) = dense(v, &n, groups, f32::from_bits(FLOAT_VECTOR_END))?;
            Values::Float { width, values }
        }
        Collect::Str(s) => Values::Str(s),
    };
    Ok(Column {
        name: f.name,
        values,
    })
}

fn get_error(e: GetError, kind: &str, tag: &str) -> Error {
    Error::invalid(
        "variants",
        format!(
            "{kind}/{tag} cannot be read as declared (htslib code {})",
            e.code()
        ),
    )
}

/// Collects the fields of a [`RecordQuery`], one record at a time.
pub(crate) struct TableBuilder<'a> {
    samples: &'a [usize],
    qual: Vec<f32>,
    filters: Strings,
    info: Vec<Field>,
    format: Vec<Field>,
    slots: Vec<usize>,
}

impl<'a> TableBuilder<'a> {
    pub(crate) fn new(
        header: &Header,
        query: &RecordQuery,
        samples: &'a [usize],
    ) -> Result<TableBuilder<'a>> {
        let declared = |name: &String, ty: Option<ValueType>, kind: &str| {
            let ty = ty.ok_or_else(|| {
                Error::invalid(
                    "variants",
                    format!("{kind}/{name} is not declared in the header"),
                )
            })?;
            Ok(Field {
                name: name.clone(),
                key: header
                    .key_id(name)
                    .expect("declared field has a dictionary ID") as usize,
                collect: Collect::of(ty),
            })
        };
        let info = query
            .info
            .iter()
            .map(|n| declared(n, header.info_type(n), "INFO"))
            .collect::<Result<Vec<_>>>()?;
        let mut format = Vec::with_capacity(query.format.len());
        for n in &query.format {
            if n == "GT" {
                return Err(Error::invalid(
                    "variants",
                    "FORMAT/GT is the genotype matrix; ask for other FORMAT fields",
                ));
            }
            let field = declared(n, header.format_type(n), "FORMAT")?;
            if matches!(field.collect, Collect::Flag(_)) {
                return Err(Error::invalid(
                    "variants",
                    format!("FORMAT/{n} is declared as a Flag, which FORMAT fields cannot be"),
                ));
            }
            format.push(field);
        }
        let slots = if info.len() > 2 || format.len() > 2 {
            vec![usize::MAX; header.keys.len()]
        } else {
            Vec::new()
        };
        Ok(TableBuilder {
            samples,
            qual: Vec::new(),
            filters: Strings::default(),
            info,
            format,
            slots,
        })
    }

    pub(crate) fn push(&mut self, header: &Header, rec: &Record) -> Result<()> {
        self.qual.push(f32::from_bits(rec.qual_bits()));
        for (i, &f) in rec.filters().iter().enumerate() {
            if i > 0 {
                self.filters.data.push(b';');
            }
            let name = header.key(f as u32).map_or("", |k| k.name.as_str());
            self.filters.data.extend_from_slice(name.as_bytes());
        }
        self.filters.end();
        // Metadata-only queries need no field lookup or collection state.
        if self.info.is_empty() && self.format.is_empty() {
            return Ok(());
        }
        // Up to two early-stop searches have linear cost and avoid the
        // position table's setup for the common small-field queries.
        let multiple = self.info.len() > 2;
        if multiple {
            positions(&mut self.slots, &self.info, &rec.data.info);
        }
        for f in &mut self.info {
            let field = if multiple {
                rec.data
                    .info
                    .get(self.slots[f.key])
                    .map(|s| info_view(&rec.data, s))
            } else {
                rec.info_field(f.key as u32)
            };
            match &mut f.collect {
                Collect::Flag(v) => v.push(field.is_some()),
                Collect::Int(v, n) => {
                    let k = field
                        .map_or(Ok(0), |field| field.numeric(|b| v.push(b as i32)))
                        .map_err(|e| get_error(e, "INFO", &f.name))?;
                    n.push(k);
                }
                Collect::Float(v, n) => {
                    let k = field
                        .map_or(Ok(0), |field| field.numeric(|b| v.push(f32::from_bits(b))))
                        .map_err(|e| get_error(e, "INFO", &f.name))?;
                    n.push(k);
                }
                Collect::Str(out) => {
                    if let Some(field) = field {
                        // The string getter copies the value count in bytes,
                        // including fields stored with a wider numeric type.
                        out.data
                            .extend_from_slice(&field.values.bytes[..field.values.len()]);
                    }
                    out.end();
                }
            }
        }
        let multiple = self.format.len() > 2;
        if multiple {
            positions(&mut self.slots, &self.format, &rec.data.fmt);
        }
        for f in &mut self.format {
            let field = if multiple {
                rec.data
                    .fmt
                    .get(self.slots[f.key])
                    .map(|s| format_view(&rec.data, s))
            } else {
                rec.format_field(f.key as u32)
            };
            if field.is_some() && rec.n_sample() != header.samples.len() {
                return Err(get_error(GetError::Unsupported, "FORMAT", &f.name));
            }
            let per = field.map_or(0, |field| field.n);
            match &mut f.collect {
                Collect::Flag(_) => unreachable!("FORMAT flags are refused"),
                Collect::Int(v, n) => {
                    if let Some(field) = field {
                        field
                            .numeric(self.samples.iter().copied(), |b| v.push(b as i32))
                            .map_err(|e| get_error(e, "FORMAT", &f.name))?;
                    }
                    n.push(per);
                }
                Collect::Float(v, n) => {
                    if let Some(field) = field {
                        field
                            .numeric(self.samples.iter().copied(), |b| v.push(f32::from_bits(b)))
                            .map_err(|e| get_error(e, "FORMAT", &f.name))?;
                    }
                    n.push(per);
                }
                Collect::Str(out) => {
                    for &s in self.samples {
                        if let Some(field) = field {
                            // bcf_get_format_char treats the payload as n
                            // bytes per sample regardless of its stored type;
                            // missing tail bytes are zero padding.
                            let at = s * per;
                            let bytes = field.values.bytes.get(at..).unwrap_or(&[]);
                            let bytes = &bytes[..per.min(bytes.len())];
                            let len = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
                            out.data.extend_from_slice(&bytes[..len]);
                        }
                        out.end();
                    }
                }
            }
        }
        Ok(())
    }

    pub(crate) fn finish(self, genotypes: GenotypeMatrix) -> Result<VariantTable> {
        let groups = self.samples.len();
        Ok(VariantTable {
            genotypes,
            qual: self.qual,
            filters: self.filters,
            info: self
                .info
                .into_iter()
                .map(|f| finish(f, 1))
                .collect::<Result<_>>()?,
            format: self
                .format
                .into_iter()
                .map(|f| finish(f, groups))
                .collect::<Result<_>>()?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VariantReader;

    const VCF: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=c,length=100>\n\
        ##INFO=<ID=I,Number=.,Type=Integer,Description=\"i\",IDX=20>\n\
        ##INFO=<ID=F,Number=.,Type=Float,Description=\"f\">\n\
        ##INFO=<ID=S,Number=1,Type=String,Description=\"s\">\n\
        ##INFO=<ID=B,Number=0,Type=Flag,Description=\"b\">\n\
        ##FORMAT=<ID=I,Number=.,Type=Integer,Description=\"i\">\n\
        ##FORMAT=<ID=F,Number=.,Type=Float,Description=\"f\">\n\
        ##FORMAT=<ID=S,Number=1,Type=String,Description=\"s\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n\
        c\t1\t.\tA\tC\t.\t.\tI=1,.;F=.,0.5;S=x;B\tI:F:S\t1,.:0.1,.:ok\t300:.:.\t70000,2,3:0.2,0.3,0.4:long\n\
        c\t2\t.\tA\tC\t.\t.\t.\tS\ta\t.\tc\n\
        c\t3\t.\tA\tC\t.\t.\tF=2;I=-200,30000;S=y\tF:I\t.:1\t1,2:30000,4\t3:70000\n";

    #[test]
    fn selected_columns_match_public_getters_for_sparse_ids_and_missing_fields() {
        for samples in [vec![], vec![2, 0, 2], vec![0, 1, 2]] {
            for names in [vec![], vec!["I"], vec!["F", "I"], vec!["S", "F", "I", "I"]] {
                let mut reader = VariantReader::new(VCF.as_bytes()).unwrap();
                let query = RecordQuery {
                    info: names
                        .iter()
                        .map(|s| (*s).to_owned())
                        .chain((names.len() > 2).then(|| "B".into()))
                        .collect(),
                    format: names.iter().map(|s| (*s).to_owned()).collect(),
                };
                let mut builder = TableBuilder::new(reader.header(), &query, &samples).unwrap();
                let mut rec = Record::default();
                let mut row = 0;
                while reader.read(&mut rec).unwrap() {
                    let header = reader.header();
                    builder.push(header, &rec).unwrap();
                    for (fields, is_format) in [(&builder.info, false), (&builder.format, true)] {
                        for f in fields {
                            match &f.collect {
                                Collect::Flag(values) => {
                                    assert_eq!(values[row], rec.info_flag(header, &f.name).unwrap())
                                }
                                Collect::Int(values, counts) => {
                                    let mut want = Vec::new();
                                    if is_format {
                                        let n = rec
                                            .format_int32(header, &f.name, &mut want)
                                            .unwrap_or(0);
                                        let width = n / header.samples.len();
                                        let selected: Vec<_> = samples
                                            .iter()
                                            .flat_map(|&s| {
                                                want[s * width..(s + 1) * width].iter().copied()
                                            })
                                            .collect();
                                        assert_eq!(counts[row], width);
                                        assert!(values.ends_with(&selected));
                                    } else {
                                        let n =
                                            rec.info_int32(header, &f.name, &mut want).unwrap_or(0);
                                        assert_eq!(counts[row], n);
                                        assert!(values.ends_with(&want));
                                    }
                                }
                                Collect::Float(values, counts) => {
                                    let mut want = Vec::new();
                                    let expected: Vec<_> = if is_format {
                                        let n = rec
                                            .format_float(header, &f.name, &mut want)
                                            .unwrap_or(0);
                                        let width = n / header.samples.len();
                                        assert_eq!(counts[row], width);
                                        samples
                                            .iter()
                                            .flat_map(|&s| {
                                                want[s * width..(s + 1) * width]
                                                    .iter()
                                                    .map(|v| v.to_bits())
                                            })
                                            .collect()
                                    } else {
                                        let n =
                                            rec.info_float(header, &f.name, &mut want).unwrap_or(0);
                                        assert_eq!(counts[row], n);
                                        want.iter().map(|v| v.to_bits()).collect()
                                    };
                                    let bits: Vec<_> = values.iter().map(|v| v.to_bits()).collect();
                                    assert!(bits.ends_with(&expected));
                                }
                                Collect::Str(values) => {
                                    let mut want = Vec::new();
                                    if is_format {
                                        let n = rec
                                            .format_char(header, &f.name, &mut want)
                                            .unwrap_or(0);
                                        let width = n / header.samples.len();
                                        for (j, &s) in samples.iter().enumerate() {
                                            let bytes = &want[s * width..(s + 1) * width];
                                            let end = bytes
                                                .iter()
                                                .position(|&b| b == 0)
                                                .unwrap_or(bytes.len());
                                            assert_eq!(
                                                values.get(row * samples.len() + j),
                                                &bytes[..end]
                                            );
                                        }
                                    } else {
                                        rec.info_string(header, &f.name, &mut want).unwrap_or(0);
                                        assert_eq!(values.get(row), want);
                                    }
                                }
                            }
                        }
                    }
                    row += 1;
                }
                assert_eq!(row, 3);
                assert_eq!(builder.qual.len(), row);
                assert_eq!(builder.filters.len(), row);
            }
        }
    }

    #[test]
    fn dense_columns_pad_each_sample_and_check_output_shape() {
        assert_eq!(
            dense(vec![1, 2, 3, 4, 5, 6], &[1, 2], 2, -1).unwrap(),
            (2, vec![1, -1, 2, -1, 3, 4, 5, 6])
        );
        assert_eq!(
            dense(vec![1, 2, 3, 4], &[1, 1], 2, -1).unwrap(),
            (1, vec![1, 2, 3, 4])
        );
        assert_eq!(dense::<i32>(vec![], &[1, 2], 0, -1).unwrap(), (2, vec![]));
        assert!(dense(vec![1], &[usize::MAX], 2, -1).is_err());
    }
}
