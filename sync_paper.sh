#!/usr/bin/env bash
# Check the paper's text against the resource models in this repository.
#
# The paper repository holds text, figures and the bibliography only. All code lives here:
# estimates/ (the models) and paper_tools/ (figure scripts and table checks). This script
# assembles a temporary tree in the layout those tools expect, with the paper's sources copied
# in (not linked, so a run is a snapshot of one paper commit), and runs the test suite and the
# table check there. It changes neither repository; figure scripts run from such a tree write to
# the copy's Figures/, not the paper's. The two commit hashes checked are printed first; a dirty
# working tree is reported so the pair in RELEASE.md can be trusted.
#
# Usage: sync_paper.sh [path-to-paper-checkout]
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER="$(cd "${1:-/home/hlamm/Desktop/fnalqc/FNAL_2026_QC_WP}" && pwd)"

TREE="$(mktemp -d)"
trap 'rm -rf "$TREE"' EXIT

echo "paper  $(git -C "$PAPER" rev-parse HEAD) $(git -C "$PAPER" describe --tags --always --dirty 2>/dev/null)"
echo "models $(git -C "$REPO" rev-parse HEAD) $(git -C "$REPO" describe --tags --always --dirty 2>/dev/null)"
rsync -a --exclude .git --exclude scripts --exclude 'main_draft.pdf' --exclude '*.aux' --exclude '*.log' \
  --exclude '*.bcf' --exclude '*.run.xml' --exclude '*.toc' --exclude '*.out' --exclude '*.blg' \
  "$PAPER/" "$TREE/"

mkdir "$TREE/scripts"
rsync -a --exclude __pycache__ --exclude .pytest_cache --exclude test_table.py --exclude test_results_table.py \
  "$REPO/paper_tools/" "$TREE/scripts/"
rsync -a --exclude __pycache__ --exclude .pytest_cache "$REPO/estimates/" "$TREE/scripts/estimates/"
cp "$REPO/paper_tools/test_table.py" "$REPO/paper_tools/test_results_table.py" "$TREE/scripts/estimates/tests/"
cp "$REPO/resources.json" "$TREE/scripts/"
cp "$REPO/docs/CONTRACT.md" "$REPO/docs/CIRCUIT_STATUS.md" "$REPO/docs/PAPER_COSTS.md" \
   "$TREE/scripts/estimates/"
# The Ch. 2 tests read their study files from applications/ beside the package.
ln -s "$REPO/applications" "$TREE/scripts/applications"

cd "$TREE/scripts"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" python3 -m pytest estimates/tests -q -p no:cacheprovider
python3 verify_table.py
