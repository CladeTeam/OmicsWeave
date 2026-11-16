use std::hint::black_box;
use std::sync::Arc;
use std::time::Instant;

use omicsweave_core::cache::BlockCache;
use omicsweave_ffm::{
    Alphabet, FfmReader, FfmWriter, Key, KeyType, Motif, MotifSet, Schema, Support, WriterOptions,
};

struct OwnedRecord {
    id: String,
    model: Motif,
    annotation: Vec<String>,
}

fn fixture() -> Vec<u8> {
    let schema = Schema::new(vec![
        Key::new("name", KeyType::String),
        Key::new("n", KeyType::Integer),
    ])
    .unwrap();
    let mut writer = FfmWriter::new(
        Vec::new(),
        Alphabet::dna(),
        schema,
        WriterOptions::default(),
    )
    .unwrap();
    for i in 0..2000 {
        let rows: Vec<Vec<f64>> = (0..10 + i % 20)
            .map(|p| {
                (0..4)
                    .map(|b| ((i * 7 + p * 3 + b) % 17 + 1) as f64)
                    .collect()
            })
            .collect();
        writer
            .push(
                &format!("id{i}"),
                &Motif::pfm(&rows, Some(100.0)),
                &[&format!("TF{i}"), "100"],
            )
            .unwrap();
    }
    writer.finish().unwrap()
}

fn load_owned(bytes: &[u8]) -> Vec<OwnedRecord> {
    let mut reader = FfmReader::open(bytes).unwrap();
    reader.set_cache(Some(Arc::new(BlockCache::new(64 << 20, 256))));
    (0..reader.len())
        .map(|i| OwnedRecord {
            id: reader.id(i).unwrap().to_owned(),
            model: reader.motif(i).unwrap().view().to_motif(),
            annotation: reader
                .annotation(i)
                .unwrap()
                .values()
                .map(str::to_owned)
                .collect(),
        })
        .collect()
}

fn owned_heap(records: &[OwnedRecord]) -> usize {
    std::mem::size_of_val(records)
        + records
            .iter()
            .map(|r| {
                r.id.capacity()
                    + r.model.factors.capacity() * std::mem::size_of::<omicsweave_ffm::Factor>()
                    + r.model
                        .factors
                        .iter()
                        .map(|f| f.positions.capacity() * 4 + f.table.capacity() * 8)
                        .sum::<usize>()
                    + match &r.model.support {
                        Support::All => 0,
                        Support::Codes(c) => c.capacity(),
                    }
                    + r.annotation.capacity() * std::mem::size_of::<String>()
                    + r.annotation.iter().map(String::capacity).sum::<usize>()
            })
            .sum::<usize>()
}

fn measure<T>(name: &str, repeats: usize, operations: usize, mut f: impl FnMut() -> T) {
    let mut samples = Vec::new();
    black_box(f());
    for _ in 0..7 {
        let start = Instant::now();
        for _ in 0..repeats {
            black_box(f());
        }
        samples.push(start.elapsed().as_nanos() as f64 / (repeats * operations) as f64);
    }
    samples.sort_by(f64::total_cmp);
    println!(
        "{name}: {:.2} ns/op (median; min {:.2}, max {:.2})",
        samples[3], samples[0], samples[6]
    );
}

fn resident_bench(bytes: &[u8]) {
    let set = MotifSet::load(bytes).unwrap();
    println!("resident heap: {} bytes", set.heap_bytes());
    measure("resident load per motif", 4, set.len(), || {
        MotifSet::load(black_box(bytes)).unwrap()
    });
    measure("resident random motif", 100, 2000, || {
        for j in 0..2000 {
            black_box(set.motif(black_box(j * 997 % 2000)).unwrap());
        }
    });
    measure("resident random tail factor", 100, 2000, || {
        for j in 0..2000 {
            let model = set.motif(black_box(j * 997 % 2000)).unwrap();
            black_box(model.factor(model.factor_count() - 1).unwrap());
        }
    });
    measure(
        "resident factor iteration per cell",
        20,
        set.value_count(),
        || {
            let mut sum = 0.0;
            for entry in set.iter() {
                for factor in entry.model().factors() {
                    for &v in factor.table {
                        sum += black_box(v);
                    }
                }
            }
            sum
        },
    );
    measure(
        "resident PFM iteration per cell",
        20,
        set.value_count(),
        || {
            let mut sum = 0.0;
            for entry in set.iter() {
                for &v in entry.model().pfm().unwrap().values() {
                    sum += black_box(v);
                }
            }
            sum
        },
    );
    measure("resident random annotation field", 100, 2000, || {
        for j in 0..2000 {
            black_box(
                set.annotation(black_box(j * 997 % 2000))
                    .unwrap()
                    .get(1)
                    .unwrap(),
            );
        }
    });
}

fn main() {
    let bytes = fixture();
    let reader = FfmReader::open(&bytes).unwrap();
    let n = reader.len();
    println!(
        "fixture: {n} motifs, {} bytes, {} blocks",
        bytes.len(),
        reader.block_count()
    );
    measure("open", 50, 1, || {
        FfmReader::open(black_box(&bytes)).unwrap()
    });
    measure("writer per motif", 1, n, fixture);
    measure("sequential decode per motif", 4, n, || {
        let reader = FfmReader::open(black_box(&bytes)).unwrap();
        for model in reader.motifs() {
            black_box(model.unwrap().view());
        }
    });
    measure("nested load per motif", 4, n, || {
        load_owned(black_box(&bytes))
    });
    println!(
        "nested model/ID/annotation heap: {} bytes",
        owned_heap(&load_owned(&bytes))
    );
    measure("uncached block decode", 100, 1, || {
        reader.motif(black_box(0)).unwrap()
    });
    let held = reader.motif(0).unwrap();
    measure("held block PFM view", 200_000, 1, || {
        held.view().pfm().unwrap()
    });
    measure("held block table slice", 200_000, 1, || {
        held.view().table_values()
    });
    measure("held block tail factor", 200_000, 1, || {
        let model = held.view();
        model.factor(model.factor_count() - 1).unwrap()
    });
    let mut cached = FfmReader::open(&bytes).unwrap();
    cached.set_cache(Some(Arc::new(BlockCache::new(64 << 20, 256))));
    for model in cached.motifs() {
        black_box(model.unwrap());
    }
    for i in 0..n {
        cached.annotation(i).unwrap();
    }
    measure("cached random motif", 100, n, || {
        for j in 0..n {
            black_box(cached.motif(black_box(j * 997 % n)).unwrap());
        }
    });
    measure("cached random tail factor", 100, n, || {
        for j in 0..n {
            let held = cached.motif(black_box(j * 997 % n)).unwrap();
            let model = held.view();
            black_box(model.factors().nth(model.factor_count() - 1).unwrap());
        }
    });
    measure("cached random annotation field", 100, n, || {
        for j in 0..n {
            black_box(cached.annotation(black_box(j * 997 % n)).unwrap().get(1));
        }
    });
    resident_bench(&bytes);
}
