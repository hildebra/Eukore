#!/bin/bash
# Eukore-bench-prep (since 2026-11-12): every download the benchmarks need, into BENCH (a folder for large files:
# ~9 GB), and the other references made comparable; the benchmarks' environment made at BENCH/env first.
#
#   bash benchmarks/eukore-bench-prep.sh BENCH PACKAGE [--genomes DIR] [--only pr2,metaphlan,...]
#
#   BENCH     the benchmarks' folder: downloads/ (PR2, MetaPhlAn 4's database, EukDetect2's database),
#             read_background/, read_genomes/ (the read benchmark's test genomes), compare/, env/, pkgs/
#   PACKAGE   a finished Eukore package (OUT/Eukore_<version>): its taxonomy labels the references, its genomes
#             give the test genomes
#   --genomes the build's genomes folder: test genomes there are linked, not downloaded
#   the rest: python3 -m eukore_bench.prep's options (--only, --read_genomes, --mpa_version)
#
# Needs the internet (a download node on a cluster). Run again to finish what failed: whole files are kept.
set -Eeuo pipefail
source "$(dirname "$0")/bench_env.sh"
BENCH=${1:?give BENCH (the benchmarks folder) and PACKAGE (OUT/Eukore_<version>)}
PACKAGE=${2:?give PACKAGE (OUT/Eukore_<version>) after BENCH}
shift 2
bench_env "$BENCH"
python3 -m eukore_bench.prep "$BENCH" --package "$PACKAGE" "$@"
