#!/usr/bin/env bash
# =============================================================================
# archive.sh  --  dated source archives of QUICHOTTE and SANCHO
#
#   ./archive.sh              both codes      (same as: make archive)
#   ./archive.sh quichotte    QUICHOTTE only  (everything except SANCHO/)
#   ./archive.sh sancho       SANCHO only     (the SANCHO/ folder)
#
# Writes into archive/ at the repository root:
#   archive/QUICHOTTE_<YYYY-MM-DD>.tar.gz
#   archive/SANCHO_<YYYY-MM-DD>.tar.gz
# A second archive on the same day gets _2, _3, ... rather than overwriting.
# Each unpacks into a folder of the same name, which contains ARCHIVE_INFO.txt
# (date, git version, commit, and whether there were uncommitted changes).
#
# What goes in: the working tree as it is now, INCLUDING uncommitted edits,
# but only files git would track -- build/, bin/, run output and everything
# else in .gitignore is left out. Outside a git repository, build products,
# test scratch and .git are left out instead.
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="$ROOT/archive"
DATE="$(date +%Y-%m-%d)"
WHICH="${1:-all}"

case "$WHICH" in
  all|quichotte|sancho) ;;
  -h|--help) sed -n '3,19p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
  *) echo "usage: $0 [all|quichotte|sancho]" >&2; exit 1 ;;
esac

IN_GIT=0
git -C "$ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1 && IN_GIT=1

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/archive.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$OUTDIR"

# list_files <dir> <name>: newline-separated paths (relative to <dir>) to archive
list_files() {
  local dir="$1" name="$2" f
  if [ "$IN_GIT" = 1 ]; then
    local exclude=()
    [ "$name" = QUICHOTTE ] && exclude=(':!SANCHO' ':!archive')
    git -C "$dir" ls-files -z --cached --others --exclude-standard -- . "${exclude[@]}" |
      while IFS= read -r -d '' f; do
        [ -e "$dir/$f" ] && printf '%s\n' "$f"   # skip files deleted but not yet committed
      done
  else
    local prune=(-path ./.git -o -path ./build -o -path ./bin -o -path ./tests/work -o -path ./archive)
    [ "$name" = QUICHOTTE ] && prune+=(-o -path ./SANCHO)
    (cd "$dir" && find . \( "${prune[@]}" \) -prune -o -type f \
        ! -name '*.o' ! -name '*.mod' ! -name '*.x' -print | sed 's|^\./||')
  fi
}

# make_archive <dir> <name>
make_archive() {
  local dir="$1" name="$2"
  local stem="${name}_${DATE}" n=2
  while [ -e "$OUTDIR/$stem.tar.gz" ]; do stem="${name}_${DATE}_$n"; n=$((n + 1)); done

  local list="$STAGE/$name.list"
  list_files "$dir" "$name" > "$list"
  if [ ! -s "$list" ]; then echo "archive.sh: no files found for $name" >&2; exit 1; fi

  # copy the files into a staging folder named after the archive, so the
  # archive unpacks into QUICHOTTE_<date>/ rather than over the current dir
  mkdir -p "$STAGE/$stem"
  (cd "$dir" && tar -cf - -T "$list") | (cd "$STAGE/$stem" && tar -xf -)

  local version="(not a git repository)" commit="-" state="-"
  if [ "$IN_GIT" = 1 ]; then
    version="$(git -C "$ROOT" describe --always --dirty --tags 2>/dev/null || echo unknown)"
    commit="$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo none)"
    if [ -n "$(git -C "$dir" status --porcelain -- . 2>/dev/null)" ]; then
      state="UNCOMMITTED CHANGES INCLUDED"
      echo "  note: $name has uncommitted changes; they are in the archive (see ARCHIVE_INFO.txt)"
    else
      state="clean (matches the commit)"
    fi
  fi
  cat > "$STAGE/$stem/ARCHIVE_INFO.txt" <<EOF
$name source archive
created   : $(date '+%Y-%m-%d %H:%M:%S %Z')
version   : $version
commit    : $commit
tree      : $state
files     : $(wc -l < "$list" | tr -d ' ')
EOF

  tar -czf "$OUTDIR/$stem.tar.gz" -C "$STAGE" "$stem"
  printf '  %-34s %s files, %s\n' "archive/$stem.tar.gz" \
    "$(wc -l < "$list" | tr -d ' ')" "$(du -h "$OUTDIR/$stem.tar.gz" | cut -f1)"
}

echo "Archiving ($DATE):"
[ "$WHICH" != sancho    ] && make_archive "$ROOT"        QUICHOTTE
[ "$WHICH" != quichotte ] && make_archive "$ROOT/SANCHO" SANCHO
exit 0
