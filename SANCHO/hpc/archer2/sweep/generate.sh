#!/bin/bash
# =============================================================================
# generate.sh -- create one run directory per manifest entry
#
#   ./generate.sh                 use manifest.txt
#   ./generate.sh other.txt       use another manifest
#   FORCE=1 ./generate.sh         overwrite input.dat of runs whose inputs changed
#
# Writes  runs/<name>/input.dat  and the list of all runs to  runs.list
# Existing runs with identical inputs are left alone (results are kept).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./sweep.conf

MANIFEST=${1:-manifest.txt}
[[ -f $MANIFEST ]] || { echo "no manifest: $MANIFEST" >&2; exit 1; }
[[ -f $DEFAULTS ]] || { echo "no defaults: $DEFAULTS" >&2; exit 1; }
mkdir -p "$RUNS_DIR"

# --- 1. expand the manifest into one line per run: "key=val key=val ..."
expanded=$(awk '
  /^[[:space:]]*(#|$)/ { next }
  !hdr { n = NF; for (i = 1; i <= NF; i++) key[i] = $i; hdr = 1; next }
  {
    if (NF != n) { printf "manifest line %d: %d columns, header has %d\n", NR, NF, n > "/dev/stderr"; bad = 1; exit 1 }
    for (i = 1; i <= n; i++) { cnt[i] = split($i, tmp, ","); for (j = 1; j <= cnt[i]; j++) v[i, j] = tmp[j]; idx[i] = 1 }
    while (1) {
      line = ""
      for (i = 1; i <= n; i++) line = line (i > 1 ? " " : "") key[i] "=" v[i, idx[i]]
      print line
      i = n; while (i >= 1 && ++idx[i] > cnt[i]) { idx[i] = 1; i-- }
      if (i < 1) break
    }
  }
  END { if (!hdr && !bad) { print "manifest has no header line" > "/dev/stderr"; exit 1 } }
' "$MANIFEST")

# value of a key in defaults.dat
default_of() { awk -v k="$1" -F= '{ key = $1; gsub(/[[:space:]]/, "", key) } key == k { v = $2; gsub(/[[:space:]]/, "", v); print v; exit }' "$DEFAULTS"; }

# directory name for a run: short, fixed order, zero-padded atom count
run_name() {
  local name="" kv k v
  for kv in $1; do
    k=${kv%%=*}; v=${kv#*=}
    case $k in
      Nsites|Natoms) printf -v v '%02d' "$v"; k=N ;;
      U_int)  k=U ;;
      t_hop)  k=t ;;
      order)  k=o ;;
    esac
    name+="${name:+_}${k}${v}"
  done
  echo "${name//\//-}"
}

: > runs.list.tmp
declare -A seen=()
n_new=0; n_same=0; n_changed=0

while read -r line; do
  [[ -z $line ]] && continue
  name=$(run_name "$line")
  [[ -n ${seen[$name]:-} ]] && { echo "duplicate run in manifest: $name" >&2; exit 1; }
  seen[$name]=1

  # --- 2. resolve Ut into U_int
  overrides=$line
  if [[ " $line " == *" Ut="* ]]; then
    [[ " $line " == *" U_int="* ]] && { echo "$name: give Ut or U_int, not both" >&2; exit 1; }
    ut=$(sed -E 's/.*(^| )Ut=([^ ]*).*/\2/' <<< "$line")
    if [[ " $line " == *" t_hop="* ]]; then t=$(sed -E 's/.*(^| )t_hop=([^ ]*).*/\2/' <<< "$line"); else t=$(default_of t_hop); fi
    u=$(awk -v r="${ut/[dD]/e}" -v t="${t/[dD]/e}" 'BEGIN { if (t < 0) t = -t; printf "%.10g", r * t }')
    overrides="$(sed -E 's/(^| )Ut=[^ ]*//' <<< "$line") U_int=$u"
  fi

  # --- 3. defaults with the overridden keys replaced (new keys appended)
  dir=$RUNS_DIR/$name
  tmp=$(mktemp)
  awk -v ov="$overrides" '
    BEGIN { n = split(ov, a, " "); for (i = 1; i <= n; i++) { p = index(a[i], "="); k = substr(a[i], 1, p - 1); o[k] = substr(a[i], p + 1); ord[i] = k } }
    /=/ && !/^[[:space:]]*#/ { k = $0; sub(/[[:space:]]*=.*/, "", k); gsub(/^[[:space:]]+/, "", k)
                               if (k in o) { print k " = " o[k]; used[k] = 1; next } }
    { print }
    END { for (i = 1; i <= n; i++) if (!(ord[i] in used)) print ord[i] " = " o[ord[i]] }
  ' "$DEFAULTS" > "$tmp"

  if [[ ! -f $dir/input.dat ]]; then
    mkdir -p "$dir"; mv "$tmp" "$dir/input.dat"; n_new=$((n_new + 1))
  elif cmp -s "$tmp" "$dir/input.dat"; then
    rm -f "$tmp"; n_same=$((n_same + 1))
  elif [[ ${FORCE:-0} == 1 ]]; then
    mv "$tmp" "$dir/input.dat"; rm -f "$dir/STATUS" "$dir/SUBMITTED"; n_changed=$((n_changed + 1))
    echo "  inputs changed, reset: $name"
  else
    rm -f "$tmp"
    echo "ERROR: $name already exists with different inputs (FORCE=1 to overwrite and reset it)" >&2
    exit 1
  fi
  echo "$name" >> runs.list.tmp
done <<< "$expanded"

mv runs.list.tmp runs.list
echo "$(wc -l < runs.list) runs in $MANIFEST: $n_new new, $n_same unchanged, $n_changed reset"
