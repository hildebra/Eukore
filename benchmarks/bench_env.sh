# shellcheck shell=bash
# Sourced by the benchmarks' scripts (since 2026-11-12): bench_env BENCH makes the benchmarks' environment at
# BENCH/env where it is not there (benchmarks/env.yaml, by micromamba, mamba or conda; their package cache at
# BENCH/pkgs, so that nothing goes to a home folder's conda directories), then puts it first on PATH and this folder
# on PYTHONPATH (python3 -m eukore_bench.<module>). No activation needed.
BENCH_HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

bench_env() {
    local bench=$1 tool
    mkdir -p "$bench"
    local env=$bench/env
    if [[ ! -x $env/bin/python3 ]]; then
        if ! tool=$(command -v micromamba || command -v mamba || command -v conda); then
            echo "no micromamba, mamba or conda to make the benchmarks' environment at $env" >&2
            exit 1
        fi
        export CONDA_PKGS_DIRS=$bench/pkgs
        echo "the benchmarks' environment: $env (packages cached in $CONDA_PKGS_DIRS)" >&2
        if [[ ${tool##*/} == micromamba ]]; then
            "$tool" create -y -p "$env" -f "$BENCH_HERE/env.yaml"
        else
            "$tool" env create -p "$env" -f "$BENCH_HERE/env.yaml"
        fi
    fi
    export PATH=$env/bin:$PATH
    export PYTHONPATH=$BENCH_HERE${PYTHONPATH:+:$PYTHONPATH}
}
