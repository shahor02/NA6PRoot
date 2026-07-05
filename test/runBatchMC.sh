#!/bin/bash
set -euo pipefail

CMD=""
NEVENTS=100
NJOBS=1
NPAR=1
SEED_START=0
OUTDIR="."
DRYRUN=0

if [[ -z $(which parallel) ]] ; then 
  echo "parallel is not available, install it first"
  exit 1
fi

usage() {
  echo "Usage: to run NJOBS na6psim <params> with at most NPARARALLEL running jobs in parallel, each with"
  echo "       producing   NEVENTS and using seed = SEED_START+jobID and writing to OUTDIR/seed<seed>:"
  echo "  runBatchMC -c 'na6psim <params>' [-n NEVENTS(=100)] [-m NPARARALLEL(=1)] [-j NJOBS(=1)] [-s SEED_START(=0)] [-o OUTDIR=(./)] [--dry]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--cmd) CMD="$2"; shift 2 ;;
    -n|--nevents) NEVENTS="$2"; shift 2 ;;
    -m|--mparallel) NPAR="$2"; shift 2 ;;
    -j|--jobs) NJOBS="$2"; shift 2 ;;
    -s|--seed-start) SEED_START="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    --dry) DRYRUN=1; shift 1 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -z "$CMD" ]] && { echo "Error: -cmd is required"; usage; }

[[ ! -d "$OUTDIR" ]] && [[ "$DRYRUN" == 0 ]] && mkdir -p "$OUTDIR"

export CMD NEVENTS OUTDIR SEED_START
echo "Will run $NJOBS jobs [$CMD] with at most $NPAR jobs in parallel, $NEVENTS each, output to $OUTDIR/job<id>_seed<seed>"
[[ $DRYRUN == 1 ]] && echo "Dry run"

run_one() {
  local jobid="$1"
  local seed=$((SEED_START + jobid))
  local odir="${OUTDIR%/}/job${jobid}_seed${seed}"
  local cmd="$CMD -n $NEVENTS --rnd-seed $seed"
  if [[ $DRYRUN == 1 ]] ; then
    echo "would run: $cmd"
  else 
    [[ ! -d "$odir" ]] && mkdir -p "$odir"
    echo "Starting [$cmd] in directory=$odir"
    (
      cd "$odir"
      eval "$cmd"
    ) > "${odir}/stdout.log" 2> "${odir}/stderr.log"
    echo "Finished job $jobid seed=$seed"
  fi
}

export -f run_one
for ((jobid=0; jobid<NJOBS; jobid++)); do
    run_one "$jobid" &

    while (( $(jobs -rp | wc -l) >= NPAR )); do
        wait -n
    done
done

wait
