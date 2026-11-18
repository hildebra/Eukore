#!/bin/bash
# Eukore's build on a SLURM cluster: the download as one job, everything after it as one chain (since 2026-11-10).
#
#   bash tools/eukore.sh fetch        # the download job, on the node with internet
#   bash tools/eukore.sh build [N]    # N process jobs (an array; default TASKS), then the final job after them
#   AFTER=<job id> bash tools/eukore.sh build   # the chain begun only once that job (the fetch) has ended well
#   REVERSION=1 bash tools/eukore.sh build      # once, on a build whose genomes an older Eukore ran: --reversion
#   bash tools/eukore.sh final        # the final job alone, no array: the genomes done, the package made again
#   bash tools/eukore.sh qc [N]      # N helpers of the package stage's vsearch checks (the final job submits them)
#   bash tools/eukore.sh bench_downloads  # benchmarks/' downloads into BENCH (a download node; AFTER=<final job>)
#   bash tools/eukore.sh benchmark   # benchmarks/ on the finished package: the comparisons, the read benchmark and
#   AFTER=<bench_downloads job> bash tools/eukore.sh benchmark   # the figures (a compute node), after the downloads
#   bash tools/eukore.sh bundle      # the benchmark run in one tar.gz (OUT/Eukore_<V>_benchmarks_evaluation.tar.gz: `benchmark` makes it too)
#
# The comparisons with other references, the read benchmark and the figures are not the build's (since
# 2026-11-12; until then the `compare` step and the final job's --read_compare): benchmarks/README.md runs them on
# the finished package, `bench_downloads` and `benchmark` here (since 2026-11-17; two jobs, each on its own
# partition) with the settings' folders: the other references'
# downloads and their comparable forms in BENCH (kept between builds), the results in BENCH_RESULTS.
#
# `build` submits both jobs at once. The process jobs share the genomes out by their claims; each, once no genome
# is left for it, goes on with the steps that take one done genome at a time under its claim (the backfill with its
# clade hold-outs, the assembly metadata, the packing of older done folders, since 2026-11-11 the organelle genes
# and the gene models' comparison), and stops
# (--stage process --genomes_only). The final job starts once every one of them has ended, however it ended
# (afterany: a task that timed out or crashed does not leave the chain pending for good), runs the genomes and
# steps they left, waits for any still claimed, makes the passes over all genomes and writes the package
# (--stage process,package, one run and one plan). Every job is resumable: `build` again goes on where the last
# chain stopped.
#
# Settings (the folders, the partitions, the build's options): tools/eukore.conf.example, copied to
# local/eukore.conf of this checkout, or to the file $EUKORE_CONF names.
#SBATCH -N 1
#SBATCH --export=ALL
set -Eeuo pipefail
trap 'rc=$?; echo "line $LINENO: \"$BASH_COMMAND\" exited $rc" >&2' ERR

STEP=${1:?give fetch, build, final, qc, bench_downloads, benchmark or bundle}
# sbatch runs a copy of this script from its spool folder: a job gets the script's own path from the submitting call
if [[ -n ${SLURM_JOB_ID:-} && -n ${EUKORE_SH:-} ]]; then SELF=$EUKORE_SH; else SELF=$(readlink -f "$0"); fi
CONF=${EUKORE_CONF:-$(dirname "$(dirname "$SELF")")/local/eukore.conf}
if [[ ! -f $CONF ]]; then
    echo "no settings at $CONF: copy tools/eukore.conf.example there, or name the file with EUKORE_CONF" >&2
    exit 1
fi
TASKS=8; QC_TASKS=0; QC_IDLE=1800; MAMBA_ENV=eukore; GEN=; SOURCES=; SCR=; LOGS=
BUILD_OPTS=(); FETCH_OPTS=(); FINAL_OPTS=(); SBATCH_ALL=(); SBATCH_FETCH=(); SBATCH_PROCESS=(); SBATCH_FINAL=(); SBATCH_QC=()
BENCH=; BENCH_RESULTS=; BENCH_DOWNLOADS_OPTS=(); BENCH_OPTS=(); FIGURE_OPTS=(); SBATCH_BENCH_DOWNLOADS=()
SBATCH_BENCH=()
source "$CONF"
: "${PDIR:?$CONF: PDIR}" "${WORK:?$CONF: WORK}" "${OUT:?$CONF: OUT}" "${VERSION:?$CONF: VERSION}"
PACKAGE=$OUT/Eukore_$VERSION
BENCH_RESULTS=${BENCH_RESULTS:-$OUT/Eukore_${VERSION}_benchmarks}
# each benchmarks' job on its own partition: the downloads where the fetch runs, the benchmarks where the processes do
(( ${#SBATCH_BENCH_DOWNLOADS[@]} )) || SBATCH_BENCH_DOWNLOADS=(${SBATCH_FETCH[@]+"${SBATCH_FETCH[@]}"})
(( ${#SBATCH_BENCH[@]} )) || SBATCH_BENCH=(${SBATCH_PROCESS[@]+"${SBATCH_PROCESS[@]}"})

submit() {  # submit NAME [sbatch's options] -- STEP: the job's ID; fails when sbatch refuses the job
    local name=$1 opts=() call jid
    shift
    while [[ $1 != -- ]]; do opts+=("$1"); shift; done
    # the logs named by the job's step (since 2026-11-16; until then eukore.<job>.*txt for every step): fetch.<job>,
    # final.<job>, and an array's process.<array>_<task>, qc.<array>_<task>
    local step=${name#Eukore-} id=%j
    step=${step,,}
    [[ " ${opts[*]} " != *" --array="* ]] || id=%A_%a
    [[ -z $LOGS ]] || opts+=(-o "$LOGS/$step.$id.otxt" -e "$LOGS/$step.$id.etxt")
    call=(sbatch --parsable -J "$name" ${SBATCH_ALL[@]+"${SBATCH_ALL[@]}"} "${opts[@]}" "$SELF" "$2")
    # checked here: submit runs inside $(...), where bash drops set -e (until 2026-11-10 a refused job gave an empty
    # ID, and `build` went on to a final job with a dependency on no job)
    if ! jid=$(EUKORE_SH=$SELF EUKORE_CONF=$CONF "${call[@]}") || [[ -z $jid ]]; then
        {
            echo "sbatch refused $name, submitted as:"
            printf ' %q' "${call[@]}"
            echo
            # 2026-11-10 on the EI cluster: not the job count, as first thought; the list qib-medium,ei-medium,
            # ei-largemem refused even a single job, since ei-largemem alone refused it
            printf '%s\n' "SLURM gives this one message for limits on the number of jobs, their size and their time;" \
                "with several partitions (-p a,b,c), one partition that refuses refuses the whole job. The same" \
                "call with --test-only, each partition alone, then one option changed at a time, finds which." \
                "The limits: sacctmgr show qos format=name%20,mintres%30,maxwall,maxsubmitpu,maxjobspu,maxtrespu%30"
        } >&2
        exit 1  # leaves the $(...) only: the caller decides
    fi
    echo "${jid%%;*}"
}

run_stage() {  # run_stage STAGE [the stage's options]: Eukore's build with the settings' folders, inside a job
    local stage=$1 cmd tiara
    shift
    echo "SLURM job ${SLURM_JOB_ID:--} (array ${SLURM_ARRAY_JOB_ID:--} task ${SLURM_ARRAY_TASK_ID:--}) on" \
        "$HOSTNAME: --stage $stage"
    activate
    cmd=(python3 "$PDIR/eukore_build.py" --build_dir "$WORK" --out "$OUT" --version "$VERSION" --stage "$stage")
    [[ -z $GEN ]] || cmd+=(--genomes "$GEN")
    [[ -z $SOURCES ]] || cmd+=(--sources "$SOURCES")
    [[ -z ${REVERSION:-} ]] || cmd+=(--reversion)
    tiara=$(tiara_of)
    [[ -z $tiara ]] || cmd+=(--tiara "$tiara")
    "${cmd[@]}" ${BUILD_OPTS[@]+"${BUILD_OPTS[@]}"} -t "${SLURM_CPUS_PER_TASK:-2}" "$@"
}

tiara_of() {  # tiara's own environment's program (it needs a Python of 3.9 or older: envs/tiara.yaml), or nothing
    # since 2026-11-14: the 2026-11 build ran without tiara ("no tiara: windows not classified"), as the eukore
    # environment holds none. TIARA_ENV in the settings names its environment (default eukore-tiara, envs/tiara.yaml's
    # name), found under micromamba's root; TIARA a program itself
    local env=${TIARA_ENV:-eukore-tiara} root=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}
    if [[ -n ${TIARA:-} ]]; then
        echo "$TIARA"
    elif [[ -x $root/envs/$env/bin/tiara ]]; then
        echo "$root/envs/$env/bin/tiara"
    elif [[ -n ${TIARA_ENV:-} ]] && command -v micromamba > /dev/null \
            && micromamba env list 2> /dev/null | awk -v e="$env" '$1 == e {print $NF}' | grep -q .; then
        echo "$(micromamba env list | awk -v e="$env" '$1 == e {print $NF}')/bin/tiara"
    fi
}

activate() {  # the settings' micromamba environment, inside a job
    set +eu  # an environment's activation scripts are not written for these
    eval "$(micromamba shell hook -s bash)"
    [[ ${CONDA_DEFAULT_ENV:-} == "$MAMBA_ENV" ]] || micromamba activate "$MAMBA_ENV"
    set -eu
}

need_bench() {  # the settings' BENCH: the benchmarks' steps need it
    if [[ -z $BENCH ]]; then
        echo "$CONF: BENCH, the benchmarks' folder (the other references, kept between builds), is not set" >&2
        exit 1
    fi
}

newest_log() {  # NEWEST_LOG: the most recently written LOGS/NAME.<job>.otxt, or empty where there is none
    local f
    NEWEST_LOG=
    for f in "$LOGS/$1".*.otxt; do  # a pattern that matches nothing stays as it is: not a file
        if [[ -f $f && ( -z $NEWEST_LOG || $f -nt $NEWEST_LOG ) ]]; then NEWEST_LOG=$f; fi
    done
}

bench_bundle() {  # the benchmark run's evaluation bundle (since 2026-11-18): OUT/Eukore_<version>_benchmarks_evaluation.tar.gz
    # the results' tables and figures, what attracts the background reads, the newest benchmark job's logs (the running
    # job's own, or, for the step `bundle` alone, the run's that made the results) and the newest bench_downloads job's:
    # a file to copy into local/ of the checkout to have the run looked at. Never fails the job.
    local logs=() name
    if [[ -n $LOGS ]]; then
        for name in benchmark bench_downloads; do
            newest_log "$name"
            [[ -z $NEWEST_LOG ]] || logs+=("$NEWEST_LOG" "${NEWEST_LOG%.otxt}.etxt")
        done
    fi
    local IFS=:
    PYTHONPATH=$PDIR/benchmarks${PYTHONPATH:+:$PYTHONPATH} python3 -m eukore_bench.bundle "$PACKAGE" "$BENCH_RESULTS" \
        --bench "$BENCH" --logs "${logs[*]-}" || echo "no benchmark bundle (eukore_bench.bundle failed)" >&2
}

need_local_tmp() {  # the node's local disk: the process and the final job need it (the shared file system's I/O)
    if [[ -z $SCR ]]; then
        echo "no local temporary folder on $HOSTNAME (SCR of $CONF is empty: SLURM_LOCAL_SCRATCH not set?)" >&2
        exit 1
    fi
}

case $STEP in
fetch)
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    jid=$(submit Eukore-fetch ${SBATCH_FETCH[@]+"${SBATCH_FETCH[@]}"} -- _fetch)
    echo "fetch: job $jid. The chain after it: AFTER=$jid bash $SELF build"
    ;;
build)
    tasks=${2:-$TASKS}
    # no chain on sources the fetch stage never filled, unless it waits for a fetch job (2026-11-11: on a new WORK,
    # every task of the array crashed in seconds, NCBI's taxonomy missing, and the final job would have too)
    if [[ -z ${AFTER:-} && ! -s ${SOURCES:-$WORK/sources}/taxdump.tar.gz ]]; then
        echo "no sources in ${SOURCES:-$WORK/sources} (its taxdump.tar.gz missing): the fetch stage first, the chain" \
            "after it: bash $SELF fetch, then AFTER=<its job id> bash $SELF build" >&2
        exit 1
    fi
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    # a fetch job that failed cancels the chain (both jobs wait for it), rather than leaving it pending
    after=${AFTER:+afterok:$AFTER}
    process=$(submit Eukore-process ${SBATCH_PROCESS[@]+"${SBATCH_PROCESS[@]}"} --array=1-"$tasks" \
        ${after:+--dependency=$after --kill-on-invalid-dep=yes} -- _process)
    # the chain whole or not at all: an array without its final job would never be packaged
    if ! final=$(submit Eukore-final ${SBATCH_FINAL[@]+"${SBATCH_FINAL[@]}"} \
        --dependency="afterany:$process${after:+,$after}" --kill-on-invalid-dep=yes -- _final); then
        scancel "$process"
        echo "array job $process cancelled: no final job after it" >&2
        exit 1
    fi
    echo "process: array job $process, $tasks tasks; final (the genomes left, the passes, the package): job $final"
    echo "to stop the chain: scancel $final $process (the final job first: it starts when the array ends)"
    ;;
final)  # the final job alone (since 2026-11-12): its genomes all done, a package stage to run again (no array)
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    after=${AFTER:+afterany:$AFTER}
    final=$(submit Eukore-final ${SBATCH_FINAL[@]+"${SBATCH_FINAL[@]}"} \
        ${after:+--dependency=$after --kill-on-invalid-dep=yes} -- _final)
    echo "final (the genomes left, the passes, the package): job $final"
    ;;
_fetch)  # unfinished downloads on the node's local disk where it has one, else beside the genomes
    run_stage fetch ${SCR:+--local_tmp "$SCR"} ${FETCH_OPTS[@]+"${FETCH_OPTS[@]}"}
    ;;
_process)
    need_local_tmp
    run_stage process --genomes_only --local_tmp "$SCR"
    ;;
_final)
    need_local_tmp
    # the vsearch checks across species shared with QC_TASKS helper jobs, which the package stage submits itself once
    # their inputs are written (since 2026-11-11)
    helpers=()
    if (( ${QC_TASKS:-0} > 0 )); then helpers=(--qc_helpers "bash $(printf %q "$SELF") qc" --qc_tasks "$QC_TASKS"); fi
    # this job's own logs, in the evaluation bundle the package stage writes (OUT/Eukore_<version>_evaluation.tar.gz)
    [[ -z $LOGS ]] || export EUKORE_LOGS="$LOGS/final.$SLURM_JOB_ID.otxt:$LOGS/final.$SLURM_JOB_ID.etxt"
    run_stage process,package --local_tmp "$SCR" ${helpers[@]+"${helpers[@]}"} \
        ${FINAL_OPTS[@]+"${FINAL_OPTS[@]}"}
    ;;
qc)  # the helpers of the package stage's searches (the loci's vsearch, since 2026-11-13 the rRNA checks', the SSU
    # models' cmsearch and the 18S cmalign parts too): QC_TASKS (or N) jobs, each searching loci till none is left
    tasks=${2:-$QC_TASKS}
    if (( tasks < 1 )); then
        echo "no helpers: QC_TASKS of $CONF is 0, and no number given (bash $SELF qc N)" >&2
        exit 1
    fi
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    jid=$(submit Eukore-qc ${SBATCH_QC[@]+"${SBATCH_QC[@]}"} --array=1-"$tasks" -- _qc)
    echo "qc: array job $jid, $tasks helpers searching the loci of $WORK/qc" >&2
    echo "$jid"  # the last line: the package stage cancels the helpers when it ends
    ;;
_qc)
    activate
    echo "SLURM job ${SLURM_JOB_ID:--} (array ${SLURM_ARRAY_JOB_ID:--} task ${SLURM_ARRAY_TASK_ID:--}) on $HOSTNAME: qc"
    cd "$PDIR"
    # the searches' temporary files on this node's local disk (SCR), where it has one (since 2026-11-17)
    if [[ -n $SCR ]] && mkdir -p "$SCR/qc"; then export TMPDIR=$SCR/qc; fi
    # waits QC_IDLE seconds for the final job's next search once none is left (since 2026-11-15: one helper array
    # serves the final job's every search), till the final job ends (WORK/qc/helpers.stop)
    python3 -m eukore.qc helper "$WORK/qc" -t "${SLURM_CPUS_PER_TASK:-2}" --idle "${QC_IDLE:-1800}"
    ;;
bench_downloads)  # the benchmarks' downloads into BENCH (since 2026-11-17): a job of its own, on a node with internet
    need_bench
    if [[ -z ${AFTER:-} && ! -f $PACKAGE/manifest.json ]]; then
        echo "no finished package at $PACKAGE: the build first (or AFTER=<its final job> bash $SELF bench_downloads)" >&2
        exit 1
    fi
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    after=${AFTER:+afterok:$AFTER}
    jid=$(submit Eukore-bench_downloads ${SBATCH_BENCH_DOWNLOADS[@]+"${SBATCH_BENCH_DOWNLOADS[@]}"} \
        ${after:+--dependency=$after --kill-on-invalid-dep=yes} -- _bench_downloads)
    echo "bench_downloads: job $jid, into $BENCH. The benchmarks after it: AFTER=$jid bash $SELF benchmark"
    ;;
benchmark)  # the comparisons, the read benchmark and the figures (since 2026-11-17): a compute node's job
    need_bench
    # no benchmarks on a BENCH no downloads job filled, unless it waits for one (as `build` for the fetch)
    if [[ -z ${AFTER:-} && ! -f $BENCH/read_genomes/tests.tsv ]]; then
        echo "no downloads in $BENCH (its read_genomes/tests.tsv missing): bash $SELF bench_downloads first, the" \
            "benchmarks after it: AFTER=<its job id> bash $SELF benchmark" >&2
        exit 1
    fi
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    # afterany: a download that failed (a server down) leaves that reference out, not the benchmarks
    after=${AFTER:+afterany:$AFTER}
    jid=$(submit Eukore-benchmark ${SBATCH_BENCH[@]+"${SBATCH_BENCH[@]}"} \
        ${after:+--dependency=$after --kill-on-invalid-dep=yes} -- _benchmark)
    echo "benchmark: job $jid, the tables and figures into $BENCH_RESULTS"
    ;;
_bench_downloads)  # every download once in BENCH (kept: a rerun fetches only what is missing), the test genomes chosen
    need_bench
    activate
    export BENCH_ENV=active  # Eukore's environment (envs/eukore.yaml): no benchmarks' environment made in BENCH
    src=${SOURCES:-$WORK/sources}
    opts=(--genomes "${GEN:-$OUT/genomes}" --sources "$src")  # SILVA's and UNITE's files of its labels linked
    # the build's own NCBI taxonomy and BUSCO dataset, the package's versions, rather than downloads of other ones
    [[ ! -s $src/taxdump.tar.gz ]] || opts+=(--taxdump "$src/taxdump.tar.gz")
    [[ ! -d $src/eukaryota_odb12/hmms ]] || opts+=(--odb "$src/eukaryota_odb12")
    bash "$PDIR/benchmarks/eukore-bench-prep.sh" "$BENCH" "$PACKAGE" "${opts[@]}" \
        ${BENCH_DOWNLOADS_OPTS[@]+"${BENCH_DOWNLOADS_OPTS[@]}"}
    ;;
_benchmark)
    need_bench
    if [[ ! -f $PACKAGE/manifest.json ]]; then  # the final job it waited for ended without one
        echo "no finished package at $PACKAGE" >&2
        exit 1
    fi
    activate
    export BENCH_ENV=active
    src=${SOURCES:-$WORK/sources}
    opts=()
    [[ ! -s $src/taxdump.tar.gz ]] || opts+=(--taxdump "$src/taxdump.tar.gz")
    [[ ! -d $WORK/align ]] || opts+=(--align_cache "$WORK/align")  # the build's proteins on BUSCO's profiles
    if [[ -n $SCR ]] && mkdir -p "$SCR/bench"; then opts+=(--tmp "$SCR/bench"); fi
    bash "$PDIR/benchmarks/eukore-benchmark.sh" "$BENCH" "$PACKAGE" "$BENCH_RESULTS" ${opts[@]+"${opts[@]}"} \
        ${BENCH_OPTS[@]+"${BENCH_OPTS[@]}"}
    # every figure (since 2026-11-17; until then the paper's alone): Figure 1 and the supplementary ones
    # (RESULTS/figures), as PNG too where FIGURE_OPTS asks another format (the website takes PNG) ...
    bash "$PDIR/benchmarks/eukore-figures.sh" "$BENCH" "$PACKAGE" "$BENCH_RESULTS" \
        ${FIGURE_OPTS[@]+"${FIGURE_OPTS[@]}"}
    if [[ " ${FIGURE_OPTS[*]-} " == *" --format "* ]]; then
        bash "$PDIR/benchmarks/eukore-figures.sh" "$BENCH" "$PACKAGE" "$BENCH_RESULTS"
    fi
    # ... the overview with the package's counts, whole and as a slide (docs/figures/overview.py; the website's Home
    # shows website/images/ copies of them) ...
    for slide in "" --slide; do
        python3 "$PDIR/docs/figures/overview.py" --package "$PACKAGE" --out "$BENCH_RESULTS/figures" --no_marks $slide
    done
    # ... and the website's Benchmarks page: its tables and figures into the build's release (OUT/release/data)
    if [[ -d $OUT/release/data ]]; then
        PYTHONPATH=$PDIR/benchmarks${PYTHONPATH:+:$PYTHONPATH} python3 -m eukore_bench.site "$BENCH_RESULTS" --release "$OUT/release" --package "$PACKAGE"
    else
        echo "no $OUT/release/data (the build's release step): the website's benchmarks not written" >&2
    fi
    # the run in one file to look at elsewhere (since 2026-11-18)
    bench_bundle
    ;;
bundle)  # the benchmark run's evaluation bundle alone, from the results already in BENCH_RESULTS (since 2026-11-18)
    need_bench
    if [[ ! -d $BENCH_RESULTS ]]; then
        echo "no benchmark results at $BENCH_RESULTS: bash $SELF benchmark first" >&2
        exit 1
    fi
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    after=${AFTER:+afterany:$AFTER}
    jid=$(submit Eukore-bundle ${SBATCH_BENCH[@]+"${SBATCH_BENCH[@]}"} \
        ${after:+--dependency=$after --kill-on-invalid-dep=yes} -- _bundle)
    echo "bundle: job $jid, into $OUT/$(basename "$BENCH_RESULTS")_evaluation.tar.gz"
    ;;
_bundle)
    need_bench
    activate
    bench_bundle
    ;;
*)
    echo "unknown step: $STEP (fetch, build, final, qc, bench_downloads, benchmark or bundle)" >&2
    exit 1
    ;;
esac
