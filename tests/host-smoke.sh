#!/usr/bin/env bash
# Run the unchanged pipeline on random, non-biological fixtures.
set -euo pipefail

RUN_DIR="${1:-$(mktemp -d "${TMPDIR:-/tmp}/viralflow-host-smoke.XXXXXX")}"
VIRALFLOW_COMMAND="${VIRALFLOW_COMMAND:-viralflow}"
mkdir -p "${RUN_DIR}"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"
[[ ! -e "${RUN_DIR}/analysis.params" ]] || {
  printf 'Choose a fresh validation directory.\n' >&2
  exit 1
}

python3 - "${RUN_DIR}" <<'PY'
import gzip
from pathlib import Path
import random
import sys

root = Path(sys.argv[1])
reads = root / "input"
reads.mkdir()
rng = random.Random(2604)
reference = "".join(rng.choices("ACGT", k=1800))
(root / "reference.fasta").write_text(
    ">synthetic Non-biological random runtime fixture\n" + reference + "\n"
)
(root / "reference.gff3").write_text(
    "##gff-version 3\n"
    "synthetic\tsynthetic\tregion\t1\t1800\t.\t+\t.\tID=synthetic_region\n"
)
complement = str.maketrans("ACGT", "TGCA")
with gzip.open(reads / "paired_R1.fq.gz", "wt") as r1, \
        gzip.open(reads / "paired_R2.fq.gz", "wt") as r2, \
        gzip.open(reads / "single.fq.gz", "wt") as single:
    for index in range(2000):
        start = rng.randrange(0, len(reference) - 350 + 1)
        forward = reference[start:start + 150]
        reverse = reference[start + 200:start + 350].translate(complement)[::-1]
        r1.write(f"@paired_{index}/1\n{forward}\n+\n{'I' * 150}\n")
        r2.write(f"@paired_{index}/2\n{reverse}\n+\n{'I' * 150}\n")
        single.write(f"@single_{index}\n{forward}\n+\n{'I' * 150}\n")
(root / "analysis.params").write_text(
    f"virus custom\ninDir {reads}\noutDir {root / 'output'}\n"
    f"referenceGenome {root / 'reference.fasta'}\n"
    f"referenceGFF {root / 'reference.gff3'}\n"
    "runSnpEff false\nnextflowSimCalls 2\n"
)
PY

cd "${RUN_DIR}"
printf 'Running host validation in %s\n' "${RUN_DIR}"
"${VIRALFLOW_COMMAND}" run --params-file "${RUN_DIR}/analysis.params" > pipeline.log 2>&1
[[ -s output/paired_results/paired.depth25.fa ]]
[[ -s output/single_results/single.depth25.fa ]]
[[ -d output/COMPILED_OUTPUT ]]
printf 'Host validation passed; log: %s/pipeline.log\n' "${RUN_DIR}"
