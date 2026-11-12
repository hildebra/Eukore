#!/bin/bash
# Eukore-figures (since 2026-11-12): every figure that needs benchmarking (Figure 1 and the supplementary ones), from
# eukore-benchmark.sh's results, into OUT/figures (with index.md, their captions and tables, and summary.md).
#
#   bash benchmarks/eukore-figures.sh BENCH PACKAGE OUT [--format pdf] [--out DIR]
#
#   BENCH     eukore-bench-prep.sh's folder (its environment)
#   PACKAGE   the package (OUT/Eukore_<version>)
#   OUT       eukore-benchmark.sh's results
set -Eeuo pipefail
source "$(dirname "$0")/bench_env.sh"
BENCH=${1:?give BENCH, PACKAGE and OUT}
PACKAGE=${2:?give PACKAGE after BENCH}
OUT=${3:?give OUT, the results of eukore-benchmark.sh, after PACKAGE}
shift 3
bench_env "$BENCH"
python3 -m eukore_bench.figures "$PACKAGE" --benchmarks "$OUT" "$@"
