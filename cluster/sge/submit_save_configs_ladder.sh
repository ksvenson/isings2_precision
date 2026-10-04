#!/usr/bin/env bash
# First production-scale save_configs ladder for the offline M_l /
# north-pole-two-point-function pipeline (2026-08-26 -- see journal.md's
# "Offline M_l[m,m'] construction..." entries). Ladder/coupling/mesh match
# production_2026-08-23 exactly (n_refine={4,8,16,32,64,128}, exact_sinh,
# naive) so the offline-built M_l/kappa2_r result is directly comparable to
# that run's known l>=3 plateau finding.
#
# N_SHARDS=96 per ladder point (576 tasks total): n_traj=20000/n_skip=2
# gives 10,000 configs/shard (this campaign's standard), so 96 shards =
# 960,000 configs/point -- above the n_blocks>=5*(2*l_max+1)=85 threshold
# scripts/symmetry_test.py's diag_chi2_const needs to actually compute
# chi2/dof (not just R_l/kappa2_r, which have no such requirement), and 4x
# the shard count of every earlier production push in this campaign.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/save_configs"
test -x "$EXE"

LADDER_CSV="4:8:16:32:64:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=96
N_TASKS=$((N_POINTS * N_SHARDS))

COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=20000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
SEED_BASE=900000  # disjoint from every prior push's seed range (800000+ used by real-space ladders)
H_RT=${H_RT:-2:00:00}

JOB_TAG="${1:-save_configs_ladder_$(date +%Y-%m-%d)}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_offline_ml_and_northpole_twopt",
  "stage": "$JOB_TAG",
  "purpose": "save raw MC configs (save_configs) for offline M_l (build_ylm_matrix_from_configs.py) and north-pole two-point function (analyze_north_pole_configs) construction, matched to production_2026-08-23's ladder/coupling/mesh for direct comparison",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "mesh_mode": "naive",
  "coupling_rule": "$COUPLING_RULE",
  "n_shards_per_point": $N_SHARDS,
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "seed_base": $SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for k in "${LADDER[@]}"; do
  for s in $(seq 0 $((N_SHARDS - 1))); do
    mkdir -p "$OUT_ROOT/q5k${k}/shard_$(printf '%03d' "$((10#$s))")"
  done
done

echo "submitting save_configs ladder array: $N_TASKS tasks ($N_POINTS ladder points x $N_SHARDS shards)"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/save_configs_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
