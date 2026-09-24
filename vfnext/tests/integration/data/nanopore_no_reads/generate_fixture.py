#!/usr/bin/env python3
"""Generate the zero-coverage NANOPORE fixture.

unmapped.fastq.gz holds reads of seeded random sequence, none of which align to
SARS-CoV-2. This is what a negative control looks like: reads come off the
flowcell, but nothing covers the reference, so the consensus must come out
fully masked rather than as a copy of the reference.

Regenerate from the `vfnext` directory with:

    python3 tests/integration/data/nanopore_no_reads/generate_fixture.py
"""

import gzip
import hashlib
import random
from pathlib import Path

SEED = 20260923
READ_COUNT = 20
READ_LENGTH = 1000
READS_FILE = "unmapped.fastq.gz"


def unmapped_reads():
    rng = random.Random(SEED)
    records = []
    for index in range(READ_COUNT):
        sequence = "".join(rng.choice("ACGT") for _ in range(READ_LENGTH))
        records.append(f"@unmapped_{index}\n{sequence}\n+\n{'5' * READ_LENGTH}\n")
    return "".join(records)


def main():
    fixture_dir = Path(__file__).resolve().parent
    reads_path = fixture_dir / READS_FILE

    # mtime=0 and no filename keep the archive byte-identical across runs.
    with reads_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(unmapped_reads().encode())

    digest = hashlib.sha256(reads_path.read_bytes()).hexdigest()
    (fixture_dir / "SHA256SUMS").write_text(f"{digest}  {READS_FILE}\n")


if __name__ == "__main__":
    main()
