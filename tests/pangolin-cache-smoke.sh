#!/usr/bin/env bash
# Exercise Snakemake's cache through the production Singularity configuration.
# No sequence input or biological analysis is used.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONTAINER="${1:-${REPO_DIR}/vfnext/containers/pangolin:4.4.sif}"
RUN_DIR="${2:-$(mktemp -d "${TMPDIR:-/tmp}/viralflow-cache-smoke.XXXXXX")}"
NEXTFLOW_COMMAND="${NEXTFLOW_COMMAND:-nextflow}"
[[ -d "${CONTAINER}" ]] || { printf 'Pangolin sandbox not found.\n' >&2; exit 1; }
CONTAINER="$(cd "${CONTAINER}" && pwd)"
mkdir -p "${RUN_DIR}"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"
[[ ! -e "${RUN_DIR}/main.nf" ]] || { printf 'Choose a fresh validation directory.\n' >&2; exit 1; }

cat > "${RUN_DIR}/cache_probe.py" <<'PY'
import json
import os
from pathlib import Path
import sys

from snakemake.common import get_appdirs
from snakemake.sourcecache import SourceCache

work_dir = Path.cwd()
cache_home = Path(os.environ["XDG_CACHE_HOME"])
assert cache_home == work_dir / ".cache", (cache_home, work_dir)
source_path = Path(get_appdirs().user_cache_dir) / "snakemake/source-cache"
assert source_path.is_relative_to(cache_home), source_path
cache = SourceCache(source_path)
for directory in (cache.cache_path, Path(cache.runtime_cache_path)):
    probe = directory / "write-probe.txt"
    probe.write_text(sys.argv[1])
    assert probe.read_text() == sys.argv[1]
Path("cache-probe.json").write_text(json.dumps({
    "task": sys.argv[1], "cache_home": str(cache_home),
    "source_cache": str(cache.cache_path), "runtime_cache": cache.runtime_cache_path,
}, indent=2) + "\n")
PY

cat > "${RUN_DIR}/main.nf" <<'NF'
nextflow.enable.dsl = 2
process runPangolin {
    tag "$caseId"
    input:
    val caseId
    path probe
    output:
    path 'cache-probe.json'
    script:
    "python ${probe} ${caseId}"
}
workflow {
    runPangolin(Channel.of('first', 'second'), file("${projectDir}/cache_probe.py"))
}
NF

# Override only the image path; containerOptions comes from the real config.
python3 - "${REPO_DIR}" "${CONTAINER}" "${RUN_DIR}" <<'PY'
import json
from pathlib import Path
import sys
repo, container, root = map(Path, sys.argv[1:])
(root / "nextflow.config").write_text(
    f"includeConfig {json.dumps(str(repo / 'vfnext/configs/containers.config'))}\n"
    f"process {{ withName:runPangolin {{ container = {json.dumps(str(container))} }} }}\n"
    "process.executor = 'local'\nprocess.maxForks = 2\n"
)
PY

cd "${RUN_DIR}"
"${NEXTFLOW_COMMAND}" run main.nf -ansi-log false -with-trace trace.txt > pipeline.log 2>&1
python3 - <<'PY'
import json
from pathlib import Path
reports = [json.loads(p.read_text()) for p in Path("work").glob("*/*/cache-probe.json")]
assert len(reports) == 2, reports
assert {r["task"] for r in reports} == {"first", "second"}, reports
assert len({r["cache_home"] for r in reports}) == 2, reports
print("Snakemake cache passed: two tasks, two writable independent caches.")
PY
printf 'Validation artifacts: %s\n' "${RUN_DIR}"
