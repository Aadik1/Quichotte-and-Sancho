#!/bin/bash
# =============================================================================
# status.sh -- one line per run: state, job, walltime, voltages done/expected
#   ./status.sh           all runs
#   ./status.sh FAILED    only runs in that state (NEW QUEUED RUNNING DONE UNCONV FAILED TIMEOUT STALE)
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"
source ./sweep.conf
[[ -f runs.list ]] || { echo "run ./generate.sh first" >&2; exit 1; }

filter=${1:-}
live=" $(squeue -h -u "${USER:-$(whoami)}" -o '%F' 2>/dev/null | sort -u | tr '\n' ' ' || true) "
declare -A count=()

printf '%-34s %-8s %-14s %9s  %s\n' RUN STATE JOB WALL_MIN VOLTAGES
while read -r name; do
  d=$RUNS_DIR/$name
  s=NEW; j=-
  if   [[ -f $d/STATUS ]];    then read -r s j < "$d/STATUS"
  elif [[ -f $d/SUBMITTED ]]; then s=QUEUED; read -r j < "$d/SUBMITTED"; fi
  [[ ($s == QUEUED || $s == RUNNING) && $live != *" ${j%_*} "* ]] && s=STALE
  count[$s]=$(( ${count[$s]:-0} + 1 ))
  [[ -n $filter && $s != "$filter" ]] && continue

  # expected voltages: same rounding as the Fortran (n_volt_steps = nint((Vf-V)/delv))
  expct=$(awk -F= '{k=$1; gsub(/[[:space:]]/,"",k); v=$2; gsub(/[[:space:]]/,"",v); sub(/[dD]/,"e",v); p[k]=v}
                   END{printf "%d", int((p["Vf"]-p["V"])/p["delv"] + 0.5)+1}' "$d/input.dat")
  order=$(awk -F= '{k=$1; gsub(/[[:space:]]/,"",k)} k=="order"{v=$2; gsub(/[[:space:]]/,"",v); print v}' "$d/input.dat")
  done_v=0; f=$d/Volt_Current_$order.dat
  [[ -f $f ]] && done_v=$(awk 'NF>=2{n++; if(NF>=3 && $3==0) u++} END{printf "%d%s", n, (u ? " (" u " unconv)" : "")}' "$f")
  wall=-; [[ -f $d/meta.txt ]] && wall=$(awk '$1=="walltime_s"{printf "%.1f", $2/60}' "$d/meta.txt")
  printf '%-34s %-8s %-14s %9s  %s\n' "$name" "$s" "$j" "$wall" "${done_v%% *}/$expct${done_v#"${done_v%% *}"}"
done < runs.list

echo
for s in NEW QUEUED RUNNING DONE UNCONV FAILED TIMEOUT STALE; do
  [[ -n ${count[$s]:-} ]] && printf '%s=%s  ' "$s" "${count[$s]}"
done; echo
