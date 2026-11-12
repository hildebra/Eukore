#!/bin/bash
# Eukore-benchmark (since 2026-11-12): a finished Eukore package against the other references and the read benchmark,
# into OUT (never the package's own folder), with the package's own benchmarks beside them.
#
#   bash benchmarks/eukore-benchmark.sh BENCH PACKAGE OUT [-t N] [--tmp DIR]
#
#   BENCH     eukore-bench-prep.sh's folder
#   PACKAGE   the package (OUT/Eukore_<version>)
#   OUT       the results: the tables and summary.json/.md; eukore-figures.sh draws them
#
# Threads: -t, else SLURM's --cpus-per-task, else 8. The read benchmark aligns reads of ~40 test genomes to the whole
# package and each other reference (minimap2): a node with tens of GB of memory; --tmp a local disk for its files.
set -Eeuo pipefail
source "$(dirname "$0")/bench_env.sh"
BENCH=${1:?give BENCH, PACKAGE and OUT}
PACKAGE=${2:?give PACKAGE after BENCH}
OUT=${3:?give OUT (the results folder) after PACKAGE}
shift 3
bench_env "$BENCH"
python3 -m eukore_bench.compare --package "$PACKAGE" --bench "$BENCH" --out "$OUT" -t "${SLURM_CPUS_PER_TASK:-8}" "$@"
