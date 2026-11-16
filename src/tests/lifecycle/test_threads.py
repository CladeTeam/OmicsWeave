from concurrent.futures import ThreadPoolExecutor

import numpy as np

from omicsweave.formats import FastaWriter
from omicsweave.query import SequenceReader


def test_threads_share_a_reader(tmp_path):
    rng = np.random.default_rng(0)
    seq = np.frombuffer(b"ACGT", np.uint8)[rng.integers(0, 4, 200_000)].tobytes()
    with FastaWriter(tmp_path / "g.fa") as writer:
        writer.add("c", seq)
    with SequenceReader(tmp_path / "g.fa") as reader:

        def work(seed):
            random = np.random.default_rng(seed)
            for _ in range(200):
                start = int(random.integers(0, len(seq)))
                end = int(random.integers(start, len(seq) + 1))
                assert reader.fetch(("c", start, end)).bases.tobytes() == seq[start:end], (
                    seed,
                    start,
                    end,
                )

        with ThreadPoolExecutor(max_workers=8) as pool:
            list(pool.map(work, range(8)))
