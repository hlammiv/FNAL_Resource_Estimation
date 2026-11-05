#!/usr/bin/env bash
# Regenerate every figure in doc/Figures/ from source.
#
# The Ch. 12 sensing figures (co-located.png, embedded.png) are deliberately
# NOT built here: they predate this pipeline, have no source, and are left
# untouched by request.
set -euo pipefail
cd "$(dirname "$0")"

FIGS=(
  fig_resource_landscape
  fig_qubit_budget
  fig_ch02_nu_nucleus
  fig_ch03_mu2e_0nubb
  fig_ch04_hybrid_lqcd
  fig_ch05_qgp_transport
  fig_ch05_budget_split
  fig_ch06_collider
  fig_ch07_curved_space
  fig_ch08_baryogenesis
  fig_ch09_chiral
  fig_ch10_finite_density
  fig_ch11_qudit
  fig_cover_wilson_hall
)

for f in "${FIGS[@]}"; do
  python3 "$f.py"
done

echo
echo "--- resource dataset ---"
python3 resources.py

echo
echo "--- Table 1.1 consistency ---"
python3 verify_table.py
