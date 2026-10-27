#!/bin/bash
# =============================================================================
# collect.sh -- gather every run's I-V curve into one table
#
#   ./collect.sh              -> results_IV.dat  (DONE and UNCONV runs)
#   ./collect.sh --all        also include partial curves (FAILED/TIMEOUT/...)
#
# Column 'conv' is 1 if the SCF converged at that voltage, 0 if not.
#
# One row per (run, voltage). Plain whitespace columns with a '#' header,
# so it loads directly in xmgrace/gnuplot/numpy.loadtxt, and one run can be
# pulled out with awk, e.g. all 12-atom curves at U/t = 2:
#     awk '$1==12 && $4==2' results_IV.dat
# Re-running it rebuilds the table from the run directories.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"
source ./sweep.conf
[[ -f runs.list ]] || { echo "run ./generate.sh first" >&2; exit 1; }

include_all=0; [[ ${1:-} == --all ]] && include_all=1
tmp=$(mktemp)
printf '# %-6s %12s %12s %10s %5s %10s %12s %22s %4s  %-8s %s\n' \
       N U_int t_hop U/t order delta V I conv state run > "$tmp"

n_runs=0; n_rows=0
while read -r name; do
  d=$RUNS_DIR/$name
  s=NEW; [[ -f $d/STATUS ]] && read -r s _ < "$d/STATUS"
  [[ $s != DONE && $s != UNCONV && $include_all == 0 ]] && continue
  order=$(awk -F= '{k=$1; gsub(/[[:space:]]/,"",k)} k=="order"{v=$2; gsub(/[[:space:]]/,"",v); print v}' "$d/input.dat")
  f=$d/Volt_Current_$order.dat
  [[ -s $f ]] || continue
  awk -v run="$name" -v state="$s" '
    FNR == NR { k = $0; sub(/[[:space:]]*=.*/, "", k); gsub(/^[[:space:]]+/, "", k)
                if (index($0, "=")) { v = substr($0, index($0, "=") + 1); gsub(/[[:space:]]/, "", v); sub(/[dD]/, "e", v); p[k] = v }
                next }
    NF >= 2 { t = p["t_hop"] + 0; at = (t < 0 ? -t : t)
              printf "  %-6d %12.6g %12.6g %10.6g %5d %10.6g %12.6f %22.14e %4s  %-8s %s\n",
                     (("Nsites" in p) ? p["Nsites"] : p["Natoms"]), p["U_int"], p["t_hop"], (at > 0 ? p["U_int"] / at : 0), p["order"], p["delta"], $1, $2, (NF >= 3 ? $3 : "-"), state, run }
  ' "$d/input.dat" "$f" >> "$tmp"
  n_runs=$((n_runs + 1)); n_rows=$((n_rows + $(grep -c . "$f")))
done < runs.list

mv "$tmp" "$RESULTS"
echo "$RESULTS: $n_rows rows from $n_runs runs"
