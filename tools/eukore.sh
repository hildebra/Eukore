#!/bin/bash
# Eukore's build on a SLURM cluster: the download as one job, everything after it as one chain (since 2026-11-10).
#
#   bash tools/eukore.sh fetch        # the download job, on the node with internet
#   bash tools/eukore.sh build [N]    # N process jobs (an array; default TASKS), then the final job after them
#   AFTER=<job id> bash tools/eukore.sh build   # the chain begun only once that job (the fetch) has ended well
#   REVERSION=1 bash tools/eukore.sh build      # once, on a build whose genomes an older Eukore ran: --reversion
#   bash tools/eukore.sh qc [N]       # N helpers of the package stage's vsearch checks (the final job submits them)
#   bash tools/eukore.sh compare      # the read benchmark's other references into COMPARE, on the node with internet
#
# The final job hands the build every other reference whose files are there (since 2026-11-10): EukDetect2's
# (EUKDETECT_FASTA, EUKDETECT_TABLE), MetaPhlAn 4's and PR2's (COMPARE/<name>.fna.gz and <name>_taxids.tsv, made by
# `compare`): the read benchmark and Figure 1 compare Eukore with them; without any, with no other reference.
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

STEP=${1:?give fetch, build or compare}
# sbatch runs a copy of this script from its spool folder: a job gets the script's own path from the submitting call
if [[ -n ${SLURM_JOB_ID:-} && -n ${EUKORE_SH:-} ]]; then SELF=$EUKORE_SH; else SELF=$(readlink -f "$0"); fi
CONF=${EUKORE_CONF:-$(dirname "$(dirname "$SELF")")/local/eukore.conf}
if [[ ! -f $CONF ]]; then
    echo "no settings at $CONF: copy tools/eukore.conf.example there, or name the file with EUKORE_CONF" >&2
    exit 1
fi
TASKS=8; QC_TASKS=0; MAMBA_ENV=eukore; GEN=; SOURCES=; SCR=; LOGS=
COMPARE=; PR2_VERSION=5.1.1; MPA_ARCHIVE=; MPA_INFO=; EUKDETECT_FASTA=; EUKDETECT_TABLE=
BUILD_OPTS=(); FETCH_OPTS=(); FINAL_OPTS=(); SBATCH_ALL=(); SBATCH_FETCH=(); SBATCH_PROCESS=(); SBATCH_FINAL=(); SBATCH_QC=()
source "$CONF"
: "${PDIR:?$CONF: PDIR}" "${WORK:?$CONF: WORK}" "${OUT:?$CONF: OUT}" "${VERSION:?$CONF: VERSION}"

submit() {  # submit NAME [sbatch's options] -- STEP: the job's ID; fails when sbatch refuses the job
    local name=$1 opts=() call jid
    shift
    while [[ $1 != -- ]]; do opts+=("$1"); shift; done
    [[ -z $LOGS ]] || opts+=(-o "$LOGS/eukore.%j.otxt" -e "$LOGS/eukore.%j.etxt")
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
    local stage=$1 cmd
    shift
    echo "SLURM job ${SLURM_JOB_ID:--} (array ${SLURM_ARRAY_JOB_ID:--} task ${SLURM_ARRAY_TASK_ID:--}) on" \
        "$HOSTNAME: --stage $stage"
    activate
    cmd=(python3 "$PDIR/eukore_build.py" --build_dir "$WORK" --out "$OUT" --version "$VERSION" --stage "$stage")
    [[ -z $GEN ]] || cmd+=(--genomes "$GEN")
    [[ -z $SOURCES ]] || cmd+=(--sources "$SOURCES")
    [[ -z ${REVERSION:-} ]] || cmd+=(--reversion)
    "${cmd[@]}" ${BUILD_OPTS[@]+"${BUILD_OPTS[@]}"} -t "${SLURM_CPUS_PER_TASK:-2}" "$@"
}

activate() {  # the settings' micromamba environment, inside a job
    set +eu  # an environment's activation scripts are not written for these
    eval "$(micromamba shell hook -s bash)"
    [[ ${CONDA_DEFAULT_ENV:-} == "$MAMBA_ENV" ]] || micromamba activate "$MAMBA_ENV"
    set -eu
}

compare_opts() {  # the --read_compare options of the other references whose files are there, EukDetect2's first
    local name fasta table
    while read -r name fasta table; do
        if [[ -n $fasta && -s $fasta && -s $table ]]; then
            printf '%s\n' --read_compare "$name=$fasta,$table"
        else
            echo "the read benchmark without $name: ${fasta:-its files not set} ${table} not there" >&2
        fi
    done <<EOF
EukDetect2 $EUKDETECT_FASTA $EUKDETECT_TABLE
MetaPhlAn4 ${COMPARE:+$COMPARE/MetaPhlAn4.fna.gz} ${COMPARE:+$COMPARE/MetaPhlAn4_taxids.tsv}
PR2 ${COMPARE:+$COMPARE/PR2.fna.gz} ${COMPARE:+$COMPARE/PR2_taxids.tsv}
EOF
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
_fetch)  # unfinished downloads on the node's local disk where it has one, else beside the genomes
    run_stage fetch ${SCR:+--local_tmp "$SCR"} ${FETCH_OPTS[@]+"${FETCH_OPTS[@]}"}
    ;;
_process)
    need_local_tmp
    run_stage process --genomes_only --local_tmp "$SCR"
    ;;
_final)
    need_local_tmp
    mapfile -t compared < <(compare_opts)
    # the vsearch checks across species shared with QC_TASKS helper jobs, which the package stage submits itself once
    # their inputs are written (since 2026-11-11)
    helpers=()
    if (( ${QC_TASKS:-0} > 0 )); then helpers=(--qc_helpers "bash $(printf %q "$SELF") qc"); fi
    run_stage process,package --local_tmp "$SCR" ${compared[@]+"${compared[@]}"} ${helpers[@]+"${helpers[@]}"} \
        ${FINAL_OPTS[@]+"${FINAL_OPTS[@]}"}
    ;;
qc)  # the helpers of the package stage's vsearch checks: QC_TASKS (or N) jobs, each searching loci till none is left
    tasks=${2:-$QC_TASKS}
    if (( tasks < 1 )); then
        echo "no helpers: QC_TASKS of $CONF is 0, and no number given (bash $SELF qc N)" >&2
        exit 1
    fi
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    jid=$(submit Eukore-qc ${SBATCH_QC[@]+"${SBATCH_QC[@]}"} --array=1-"$tasks" -- _qc)
    echo "qc: array job $jid, $tasks helpers searching the loci of $WORK/qc" >&2
    echo "$jid"  # the last line: the package stage cancels the helpers once every locus is searched
    ;;
_qc)
    activate
    echo "SLURM job ${SLURM_JOB_ID:--} (array ${SLURM_ARRAY_JOB_ID:--} task ${SLURM_ARRAY_TASK_ID:--}) on $HOSTNAME: qc"
    cd "$PDIR"
    python3 -m eukore.qc helper "$WORK/qc" -t "${SLURM_CPUS_PER_TASK:-2}"
    ;;
compare)
    : "${COMPARE:?$CONF: COMPARE, the folder of the other references of the read benchmark}"
    [[ -z $LOGS ]] || mkdir -p "$LOGS"
    jid=$(submit Eukore-compare ${SBATCH_FETCH[@]+"${SBATCH_FETCH[@]}"} -- _compare)
    echo "compare: job $jid (the 18S of PR2${MPA_ARCHIVE:+ and the markers of MetaPhlAn 4} into $COMPARE)"
    ;;
_compare)  # PR2's 18S (downloaded) and MetaPhlAn 4's eukaryote markers (its archive given): tools/compare_references.py
    activate
    mkdir -p "$COMPARE"
    pr2=$COMPARE/pr2_version_${PR2_VERSION}_SSU_taxo_long.fasta.gz
    if [[ ! -s $pr2 ]]; then
        curl -fsSL -o "$pr2.part" \
            "https://github.com/pr2database/pr2database/releases/download/v$PR2_VERSION/${pr2##*/}"
        mv "$pr2.part" "$pr2"
    fi
    python3 "$PDIR/tools/compare_references.py" pr2 --fasta "$pr2" --taxdump "${SOURCES:-$WORK/sources}/taxdump.tar.gz" \
        --out "$COMPARE"
    if [[ -n $MPA_ARCHIVE && -n $MPA_INFO ]]; then
        python3 "$PDIR/tools/compare_references.py" metaphlan --archive "$MPA_ARCHIVE" --marker_info "$MPA_INFO" \
            --out "$COMPARE"
    else
        echo "no MetaPhlAn 4 reference: MPA_ARCHIVE and MPA_INFO of $CONF not set" >&2
    fi
    ;;
*)
    echo "unknown step: $STEP (fetch, build or compare)" >&2
    exit 1
    ;;
esac
