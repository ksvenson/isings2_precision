#!/usr/bin/env bash
# push6: replaces the n_refine={32,64,128} portion of push5 after killing
# it mid-run 2026-08-24, per user direction to prioritize wall-clock over
# push5's originally-planned per-shard statistics. n_refine={8,16} are NOT
# rerun here -- push5 already finished those (32 shards x 1e6 traj/shard
# each, untouched, sitting in campaign_runs/push5_2026-08-24_lmax12/).
#
# Key change from push5: user's actual statistics target is 1,000,000
# TOTAL pooled trajectories per ladder point (not per shard). Since shards
# are independent parallel Markov chains and wall time per shard scales
# with n_traj (not with how many other shards run alongside it), splitting
# a fixed total budget across more/shorter shards is a real wall-clock win
# with no change in total statistics:
# - N_TRAJ=31250/shard (was 1,000,000) -- 32 shards x 31250 = 1,000,000
#   pooled, matching the user's stated target exactly.
# - N_SHARDS=32 kept unchanged (not increased): push5 already proved 32
#   tasks schedule immediately with no queueing on this project's
#   allocation: further increasing shard count doesn't reduce wall time
#   once already fully parallel, it would just risk queueing delay for no
#   benefit.
# - JACK_BLOCK_SIZE=1250 -- native blocks/shard = 31250/1250 = 25; pooled
#   across 32 shards (this campaign's standard reducer convention: pool
#   raw per-block sums in ascending shard order before jackknifing) gives
#   800 total blocks, still well above PLAN.md's Phase 1 gate
#   n_blocks >= 10*(2*l_max+1) = 250 at l_max=12 (3.2x margin) and above
#   the symmetry_test.py --min-blocks-factor floor of 5*(2l+1)=125.
# - Estimated wall time (0.054 s/sweep calibration from push5, scaled by
#   site count for n_refine<128): n_refine=128 ~58min, n_refine=64 ~15min,
#   n_refine=32 ~4min -- vs. push5's ~30h/~28.5h-remaining estimate for
#   the same three points at the old n_traj. h_rt kept generous (4:00:00)
#   for margin, not because the estimate is expected to be exceeded.
# - Tradeoff, stated explicitly: this is 32x fewer pooled trajectories
#   than push5's original per-shard design (1M vs 32M pooled per point).
#   Still formally clears the Phase 1 gate with margin, but error bars on
#   kappa2_r/R_l etc. will be larger than what the killed push5 run would
#   have produced. User confirmed this tradeoff explicitly before this
#   script was written.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/ising_s2_crit"
test -x "$EXE"

LADDER_CSV="32:64:128"
N_SHARDS=32
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_TASKS=$((N_POINTS * N_SHARDS))

L_MAX=12
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=31250
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=1250
SEED_BASE=600000   # disjoint from production_2026-08-23 (90000+96) and push5 (500000+160)
H_RT=${H_RT:-4:00:00}

JOB_TAG="${1:-push6_$(date +%Y-%m-%d)_lmax12_1Mpooled}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_phase1_production",
  "stage": "$JOB_TAG",
  "supersedes": "push5_2026-08-24_lmax12 n_refine={32,64,128} (killed mid-run, tasks 65-160 of job 7301188)",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "n_shards": $N_SHARDS,
  "pooled_total_traj_per_point": $((N_SHARDS * N_TRAJ)),
  "l_max": $L_MAX,
  "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "seed_base": $SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for k in "${LADDER[@]}"; do
  for s in $(seq -w 0 $((N_SHARDS - 1))); do
    mkdir -p "$OUT_ROOT/q5k${k}/shard_${s}/q5k${k}"
  done
done

echo "submitting $N_TASKS-task push6 array ($N_POINTS ladder points x $N_SHARDS shards), l_max=$L_MAX, 1M pooled traj/point, project qfe, unpinned queue, h_rt=$H_RT"
qsub \
  -P qfe \
  -N "s2prec_$JOB_TAG" \
  -j y \
  -o "$LOG_DIR/" \
  -t "1-$N_TASKS" \
  -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/production_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
