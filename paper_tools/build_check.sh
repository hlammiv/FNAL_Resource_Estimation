#!/bin/bash
# Dual-view build + verification for the tracked-changes workflow.
# Run from doc/:  bash scripts/build_check.sh
# Tracked view gates: 0 errors, 0 undefined refs, verify_table PASS
#   (overfull recorded, not gated -- redline markup perturbs line breaking).
# Clean view gates: 0 errors, 0 overfull, 0 undefined refs, verify_table PASS.
#   The clean-view signature line is what goes in CHANGELOG.md.
set -u
cd "$(dirname "$0")/.."

run_chain () {  # $1 = jobname, $2 = pre-tex definition ('' for tracked)
  local job=$1 pre=$2
  pdflatex -interaction=nonstopmode -halt-on-error -jobname="$job" \
    "${pre}\input{main.tex}" > /dev/null 2>&1 || { echo "$job: pdflatex FAILED (pass 1)"; return 1; }
  biber "$job" > /dev/null 2>&1 || { echo "$job: biber FAILED"; return 1; }
  pdflatex -interaction=nonstopmode -halt-on-error -jobname="$job" \
    "${pre}\input{main.tex}" > /dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error -jobname="$job" \
    "${pre}\input{main.tex}" > /dev/null 2>&1 || { echo "$job: pdflatex FAILED (final pass)"; return 1; }
  local errs undef over pages
  errs=$(grep -c '^!' "$job.log")
  undef=$(grep -ci 'undefined' "$job.log")
  over=$(grep -c 'Overfull' "$job.log")
  pages=$(pdfinfo "$job.pdf" | awk '/^Pages/{print $2}')
  echo "$job: ${pages} pp, ${errs} errors, ${over} overfull, ${undef} undefined refs"
}

echo "== tree containment =="
# Every built source must live inside doc/: the upload package is made from doc/ alone.
# (2026-09-29: main.tex had been building Ch. 6 from ../ch6_rewrite/ for four weeks.)
esc=$(grep -nE '\\(input|include|includegraphics|addbibresource)(\[[^]]*\])?\{[^}]*\.\./' \
        main.tex main-overview.tex commands.tex applications/*.tex 2>/dev/null | grep -v '^[^:]*:[0-9]*:[[:space:]]*%')
if [ -n "$esc" ]; then echo "FAIL: built sources reference paths outside doc/:"; echo "$esc"; exit 1; fi
echo "ok: no built source leaves doc/"
echo "== tracked view (main) =="
run_chain main "" || exit 1
echo "== clean view (main-clean) =="
run_chain main-clean "\def\TrackedView{0}" || exit 1
echo "== verify_table =="
python3 scripts/verify_table.py || exit 1
echo "== estimates (chapter models vs boxes) =="
( cd scripts && python3 -m pytest estimates/tests -q 2>&1 | tail -1 )
echo "== tracked edits =="
python3 scripts/accept_changes.py --list | tail -1
