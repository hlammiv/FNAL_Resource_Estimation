#!/usr/bin/env bash
# Push the resource models from this repository into the paper's working copy and run the
# paper's own test suite there.
#
# This repository is the source of truth for estimates/. The paper tree keeps a copy at
# scripts/estimates/ so that its figures (scripts/fig_*.py, resources.py) and its table checks
# (tests/test_table.py, verify_table.py) can import it. Files that live only in the paper tree
# are excluded, so --delete does not remove them:
#   NEEDS_AUTHOR.md, proposals/         the authors' working notes
#   CONTRACT.md, CIRCUIT_STATUS.md,     the paper's copies of docs/ (refreshed below)
#   PAPER_COSTS.md
#   tests/test_table.py                 needs scripts/verify_table.py and the paper's .tex
set -euo pipefail

REPO=/home/hlamm/Desktop/fnalqc/FNAL_Resource_Estimation
PAPER=/home/hlamm/Desktop/fnalqc/overleaf_upload_20261002b

rsync -a --delete \
  --exclude __pycache__ --exclude .pytest_cache \
  --exclude NEEDS_AUTHOR.md --exclude proposals/ \
  --exclude CONTRACT.md --exclude CIRCUIT_STATUS.md --exclude PAPER_COSTS.md \
  --exclude tests/test_table.py \
  "$REPO/estimates/" "$PAPER/scripts/estimates/"

# The three generated/contract documents are kept in docs/ here and beside the package there.
cp "$REPO/docs/CONTRACT.md" "$REPO/docs/CIRCUIT_STATUS.md" "$REPO/docs/PAPER_COSTS.md" \
   "$PAPER/scripts/estimates/"

cd "$PAPER/scripts" && OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" python3 -m pytest estimates/tests -q
