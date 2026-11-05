#!/usr/bin/env bash
# Check the paper's text against the resource models in this repository.
#
# The paper repository holds text, figures and the bibliography only. All code lives here:
# estimates/ (the models) and paper_tools/ (figure scripts and table checks). This script
# assembles a temporary tree in the layout those tools expect, with the paper's sources linked
# in, and runs the test suite and the table check there. It changes neither repository, except
# that figure scripts run from such a tree write to the paper's Figures/.
#
# Usage: sync_paper.sh [path-to-paper-checkout]
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER="$(cd "${1:-/home/hlamm/Desktop/fnalqc/FNAL_2026_QC_WP}" && pwd)"

TREE="$(mktemp -d)"
trap 'rm -rf "$TREE"' EXIT

for f in "$PAPER"/*; do ln -s "$f" "$TREE/"; done
rm -f "$TREE/scripts"

mkdir "$TREE/scripts"
rsync -a --exclude __pycache__ --exclude .pytest_cache --exclude test_table.py \
  "$REPO/paper_tools/" "$TREE/scripts/"
rsync -a --exclude __pycache__ --exclude .pytest_cache "$REPO/estimates/" "$TREE/scripts/estimates/"
cp "$REPO/paper_tools/test_table.py" "$TREE/scripts/estimates/tests/"
cp "$REPO/resources.json" "$TREE/scripts/"
cp "$REPO/docs/CONTRACT.md" "$REPO/docs/CIRCUIT_STATUS.md" "$REPO/docs/PAPER_COSTS.md" \
   "$TREE/scripts/estimates/"
# The Ch. 2 tests read their study files from applications/ beside the package.
ln -s "$REPO/applications" "$TREE/scripts/applications"

cd "$TREE/scripts"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" python3 -m pytest estimates/tests -q -p no:cacheprovider
python3 verify_table.py
