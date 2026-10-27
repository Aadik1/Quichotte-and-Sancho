#!/bin/bash -l
# =============================================================================
# job_array.sh -- one Slurm array task = one ARCHER2 node running up to
# RUNS_PER_NODE runs side by side, each as its own srun step.
# Don't sbatch this directly: submit.sh sets the resources, the list of runs
# (RUN_LIST) and the folder (SWEEP_ROOT).
#
# Exit codes of the code -> run state:
#   0 DONE      3 UNCONV (finished, some voltages not converged)
#   other FAILED           time limit -> TIMEOUT
# =============================================================================
#SBATCH --signal=B:USR2@300          # warn the script 5 min before the time limit

set -u
cd "$SWEEP_ROOT"
source ./sweep.conf

module load $MODULES

export OMP_NUM_THREADS=$CPUS
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

task=$SLURM_ARRAY_TASK_ID
first=$(( (task - 1) * RUNS_PER_NODE + 1 ))
last=$((  task * RUNS_PER_NODE ))
mapfile -t names < <(sed -n "${first},${last}p" "$RUN_LIST")
(( ${#names[@]} > 0 )) || { echo "no runs for task $task in $RUN_LIST"; exit 1; }

id="${SLURM_ARRAY_JOB_ID}_${task}"
flag=$SWEEP_ROOT/$LOG_DIR/.timeout_$id
rm -f "$flag"

run_one() {   # $1 = run name; runs in the background, writes STATUS and meta.txt
  local name=$1 d=$SWEEP_ROOT/$RUNS_DIR/$1 start end rc state
  cd "$d" || return
  rm -f Volt_Current_*.dat scf_V_*.dat populations.dat spectral_function.dat run.log
  echo "RUNNING $id" > STATUS
  start=$(date +%s)
  srun --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --mem="$MEM_PER_RUN" \
       --exact --hint=nomultithread --distribution=block:block \
       "$BIN" > run.log 2>&1
  rc=$?
  end=$(date +%s)
  if   [[ -f $flag ]]; then state=TIMEOUT
  elif [[ $rc == 0 ]]; then state=DONE
  elif [[ $rc == 3 ]]; then state=UNCONV
  else                      state=FAILED
  fi
  echo "$state $id" > STATUS
  cat > meta.txt <<META
run         $name
state       $state
exit_code   $rc
job         $id
host        $(hostname)
cpus        $CPUS
start       $(date -d @$start '+%F %T')
end         $(date -d @$end '+%F %T')
walltime_s  $((end - start))
binary      $BIN
binary_md5  $(md5sum "$BIN" | cut -d' ' -f1)
binary_date $(date -r "$BIN" '+%F %T')
META
  echo "$name: $state (exit $rc) after $((end - start)) s"
}

pids=()
# near the time limit: mark the node as timed out and stop every run on it
trap 'touch "$flag"; for p in "${pids[@]}"; do pkill -TERM -P "$p" 2>/dev/null; done' USR2

for name in "${names[@]}"; do
  run_one "$name" &
  pids+=($!)
done

for p in "${pids[@]}"; do
  while :; do wait "$p"; (( $? > 128 )) && kill -0 "$p" 2>/dev/null && continue; break; done
done
rm -f "$flag"
echo "node task $id finished: ${#names[@]} run(s)"
