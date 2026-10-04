#!/usr/bin/env bash
# Run the ising_s2_crit measurement driver across a n_refine ladder at fixed
# q/l_max/coupling_rule/statistics, creating each run's data_dir first (the
# driver does not mkdir -p its output directory).
#
# Usage: scripts/run_ladder.sh <data_dir> <seed> <n_therm> <n_traj> <n_skip> <l_max> <n_refine...>
# Example (this session's push):
#   scripts/run_ladder.sh campaign_runs/push_2026-08-22 4242 2000 20000 2 3 2 4 8 16 32
set -euo pipefail

data_dir=$1; shift
seed=$1; shift
n_therm=$1; shift
n_traj=$1; shift
n_skip=$1; shift
l_max=$1; shift
n_refines=("$@")

cd "$(dirname "$0")/.."

for k in "${n_refines[@]}"; do
  mkdir -p "$data_dir/q5k$k"
  bin/ising_s2_crit --q 5 --n_refine "$k" --l_max "$l_max" \
    --coupling_rule exact_sinh --data_dir "$data_dir" \
    --n_therm "$n_therm" --n_traj "$n_traj" --n_skip "$n_skip" --seed "$seed"
done
