use std::sync::atomic::{AtomicUsize, Ordering};

use crate::support::{Storage, expected, spec, write_and_open};
use omicsweave_cool::{Cool, PixelBatch, PixelScratch};

#[test]
fn shared_reader_matches_serial_across_native_and_library_paths() {
    fn send_sync<T: Send + Sync>() {}
    fn send<T: Send>() {}
    send_sync::<Cool>();
    send::<PixelScratch>();

    for storage in [
        [Storage::gzip_shuffle(64); 3],
        [Storage::Contiguous; 3],
        [
            Storage::gzip_shuffle(64),
            Storage::Contiguous,
            Storage::ScaleOffset { chunk: 77 },
        ],
    ] {
        let dir = tempfile::tempdir().unwrap();
        let mut s = spec(20_000, 14);
        s.storage = storage;
        let cool = write_and_open(dir.path(), &s);
        let spans = cool.pixel_spans(700);
        assert!(spans.len() > 8);
        let next = AtomicUsize::new(0);
        let results = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(|| {
                        let mut scratch = PixelScratch::new();
                        let mut batch = PixelBatch::new();
                        let mut results = Vec::new();
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            let Some(span) = spans.get(i) else { break };
                            cool.read_pixels(span.clone(), &mut scratch, &mut batch)
                                .unwrap_or_else(|err| {
                                    panic!("storage {storage:?}, span {span:?}: {err}")
                                });
                            results.push((i, batch.clone()));
                        }
                        results
                    })
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        let mut seen = vec![false; spans.len()];
        for (i, batch) in results {
            assert!(!std::mem::replace(&mut seen[i], true), "duplicate span {i}");
            assert_eq!(
                batch,
                expected(&s, spans[i].clone()),
                "storage {storage:?}, span {:?}",
                spans[i]
            );
        }
        assert!(
            seen.into_iter().all(|read| read),
            "storage {storage:?}: missing span"
        );
    }
}
