#!/bin/bash
# Eukore's build on a SLURM cluster: the download as one job, everything after it as one chain (since 2026-11-10).
#
#   bash tools/eukore.sh fetch        # the download job, on the node with internet
#   bash tools/eukore.sh build [N]    # N process jobs (an array; default TASKS), then the final job after them
#   AFTER=<job id> bash tools/eukore.sh build   # the chain begun only once that job (the fetch) has ended well
#   REVERSION=1 bash tools/eukore.sh build      # once, on a build whose genomes an older Eukore ran: --reversion
#
# `build` submits both jobs at once. The process jobs share the genomes out by their claims; each, once no genome
# is left for it, goes on with the steps that take one done genome at a time under its claim (the backfill with its
# clade hold-outs, the assembly metadata, the packing of older done folders), and stops
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

STEP=${1:?give fetch or build}
# sbatch runs a copy of this script from its spool folder: a job gets the script's own path from the submitting call
if [[ -n ${SLURM_JOB_ID:-} && -n ${EUKORE_SH:-} ]]; then SELF=$EUKORE_SH; else SELF=$(readlink -f "$0"); fi
CONF=${EUKORE_CONF:-$(dirname "$(dirname "$SELF")")/local/eukore.conf}
if [[ ! -f $CONF ]]; then
    echo "no settings at $CONF: copy tools/eukore.conf.example there, or name the file with EUKORE_CONF" >&2
    exit 1
fi
TASKS=8; MAMBA_ENV=eukore; GEN=; SCR=; LOGS=
BUILD_OPTS=(); FETCH_OPTS=(); FINAL_OPTS=(); SBATCH_ALL=(); SBATCH_FETCH=(); SBATCH_PROCESS=(); SBATCH_FINAL=()
source "$CONF"
: "${PDIR:?$CONF: PDIR}" "${WORK:?$CONF: WORK}" "${OUT:?$CONF: OUT}" "${VERSION:?$CONF: VERSION}"

submit() {  # submit NAME [sbatch's options] -- STEP: the job's ID
    local name=$1 opts=() jid
    shift
    while [[ $1 != -- ]]; do opts+=("$1"); shift; done
    [[ -z $LOGS ]] || opts+=(-o "$LOGS/eukore.%j.otxt" -e "$LOGS/eukore.%j.etxt")
    jid=$(EUKORE_SH=$SELF EUKORE_CONF=$CONF sbatch --parsable -J "$name" ${SBATCH_ALL[@]+"${SBATCH_ALL[@]}"} \
        "${opts[@]}" "$SELF" "$2")
    echo "${jid%%;*}"
}

run_stage() {  # run_stage STAGE [the stage's options]: Eukore's build with the settings' folders, inside a job
    local stage=$1 cmd
    shift
    echo "SLURM job ${SLURM_JOB_ID:--} (array ${SLURM_ARRAY_JOB_ID:--} task ${SLURM_ARRAY_TASK_ID:--}) on" \
        "$HOSTNAME: --stage $stage"
    set +eu  # an environment's activation scripts are not written for these
    eval "$(micromamba shell hook -s bash)"
    [[ ${CONDA_DEFAULT_ENV:-} == "$MAMBA_ENV" ]] || micromamba activate "$MAMBA_ENV"
    set -eu
    cmd=(python3 "$PDIR/eukore_build.py" --build_dir "$WORK" --out "$OUT" --version "$VERSION" --stage "$stage")
    [[ -z $GEN ]] || cmd+=(--genomes "$GEN")
    [[ -z ${REVERSION:-} ]] || cmd+=(--reversion)
    "${cmd[@]}" ${BUILD_OPTS[@]+"${BUILD_OPTS[@]}"} -t "${SLURM_CPUS_PER_TASK:-2}" "$@"
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
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    # a fetch job that failed cancels the chain (both jobs wait for it), rather than leaving it pending
    after=${AFTER:+afterok:$AFTER}
    process=$(submit Eukore-process ${SBATCH_PROCESS[@]+"${SBATCH_PROCESS[@]}"} --array=1-"$tasks" \
        ${after:+--dependency=$after --kill-on-invalid-dep=yes} -- _process)
    final=$(submit Eukore-final ${SBATCH_FINAL[@]+"${SBATCH_FINAL[@]}"} \
        --dependency="afterany:$process${after:+,$after}" --kill-on-invalid-dep=yes -- _final)
    echo "process: array job $process, $tasks tasks; final (the genomes left, the passes, the package): job $final"
    echo "to stop the chain: scancel $final $process (the final job first: it starts when the array ends)"
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
    run_stage process,package --local_tmp "$SCR" ${FINAL_OPTS[@]+"${FINAL_OPTS[@]}"}
    ;;
*)
    echo "unknown step: $STEP (fetch or build)" >&2
    exit 1
    ;;
esac
