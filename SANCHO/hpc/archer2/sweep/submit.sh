#!/bin/bash
# =============================================================================
# submit.sh -- submit runs as one Slurm job array
#
#   ./submit.sh            runs never submitted            (default)
#   ./submit.sh failed     FAILED / TIMEOUT runs, and runs whose job vanished
#                          (UNCONV runs are not included: resubmitting them as
#                          they are gives the same answer. Change their inputs,
#                          FORCE=1 ./generate.sh, then submit them by name.)
#   ./submit.sh all        both of the above
#   ./submit.sh N08_Ut1.0_o2 N12_Ut1.0_o2    just these runs
#   DRY=1 ./submit.sh ...  show what would be submitted, submit nothing
#
# Runs that are DONE, or QUEUED/RUNNING in a live job, are never resubmitted.
# Runs are packed RUNS_PER_NODE to a node: one array task = one node.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
SWEEP_ROOT=$(pwd)
source ./sweep.conf

[[ -f runs.list ]] || { echo "run ./generate.sh first" >&2; exit 1; }
[[ -x $BIN ]] || { echo "executable not found: $BIN" >&2; exit 1; }
mkdir -p "$LOG_DIR" lists
(( RUNS_PER_NODE * CPUS <= 128 )) || { echo "RUNS_PER_NODE*CPUS = $((RUNS_PER_NODE*CPUS)) > 128 cores" >&2; exit 1; }

mode=${1:-new}
# jobs Slurm still knows about (array job ids)
live=" $(squeue -h -u "${USER:-$(whoami)}" -o '%F' 2>/dev/null | sort -u | tr '\n' ' ' || true) "

state_of() {  # prints the effective state of a run
  local d=$RUNS_DIR/$1 s j
  if   [[ -f $d/STATUS ]];    then read -r s j < "$d/STATUS"
  elif [[ -f $d/SUBMITTED ]]; then s=QUEUED; read -r j < "$d/SUBMITTED"
  else echo NEW; return; fi
  if [[ $s == QUEUED || $s == RUNNING ]] && [[ $live != *" ${j%_*} "* ]]; then echo STALE; else echo "$s"; fi
}

selected=()
case $mode in
  new|failed|all)
    while read -r name; do
      st=$(state_of "$name")
      case $mode:$st in
        new:NEW|failed:FAILED|failed:TIMEOUT|failed:STALE|all:NEW|all:FAILED|all:TIMEOUT|all:STALE) selected+=("$name") ;;
      esac
    done < runs.list ;;
  *)
    for name in "$@"; do
      grep -qx "$name" runs.list || { echo "not in runs.list: $name" >&2; exit 1; }
      st=$(state_of "$name")
      [[ $st == QUEUED || $st == RUNNING ]] && { echo "skipping $name: $st" >&2; continue; }
      selected+=("$name")
    done ;;
esac

n=${#selected[@]}
(( n > 0 )) || { echo "nothing to submit ($mode)"; exit 0; }

list=lists/$(date +%Y%m%d_%H%M%S)_$mode.list
printf '%s\n' "${selected[@]}" > "$list"

echo "$n run(s) -> $list"
if [[ ${DRY:-0} == 1 ]]; then cat "$list"; rm -f "$list"; exit 0; fi

# clear the old state first; the array tasks write STATUS themselves
for name in "${selected[@]}"; do rm -f "$RUNS_DIR/$name/STATUS" "$RUNS_DIR/$name/SUBMITTED"; done

ntasks=$(( (n + RUNS_PER_NODE - 1) / RUNS_PER_NODE ))
jobid=$(sbatch --parsable \
  --job-name="$JOB_NAME" --account="$ACCOUNT" --partition="$PARTITION" --qos="$QOS" \
  --time="$TIME" --nodes=1 --exclusive \
  --array=1-"$ntasks"%"$MAX_PARALLEL" \
  --output="$LOG_DIR/%x.%A_%a.log" \
  --export=ALL,SWEEP_ROOT="$SWEEP_ROOT",RUN_LIST="$SWEEP_ROOT/$list" \
  job_array.sh)
jobid=${jobid%%;*}

i=0
for name in "${selected[@]}"; do
  echo "${jobid}_$(( i / RUNS_PER_NODE + 1 ))" > "$RUNS_DIR/$name/SUBMITTED"; i=$((i + 1))
done
echo "submitted array job $jobid: $n run(s) on $ntasks node(s), $RUNS_PER_NODE per node, at most $MAX_PARALLEL nodes at once"
