#!/usr/bin/env bash
# Single-task job: naive vs equal_area vs equal_rp kappa2_r comparison
# over the full lean_mesh_compare_100shards_2026-08-27 ladder (up to
# n_refine=128, ~10,000 blocks at the top point) -- too slow/heavy for
# the login node (got OOM/resource-killed there), moved to SGE per
# standing convention (prefer detached qsub over interactive runs for
# anything beyond a quick spot check).
set -euo pipefail

ROOT="/projectnb/qfe/misra/IsingS2_precision"
source "$ROOT/.venv_plot/bin/activate"
python3 "$ROOT/scripts/compare_three_mesh_kappa2r.py" --l_max 8
