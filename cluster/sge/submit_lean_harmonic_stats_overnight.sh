#!/usr/bin/env bash
# Overnight lean_harmonic_stats push, per user direction 2026-08-26
# ("schedule a much larger stats ladder with larger lattices for an
# overnight run") -- supersedes submit_lean_harmonic_stats.sh's ladder
# (4..32, n_traj=100000): extends the ladder up to n_refine=128 (matching
# this campaign's standard production ladder, e.g.
# save_configs_ladder_2026-08-26) and doubles the per-point statistics.
# Copied from submit_lean_harmonic_stats.sh rather than editing it in
# place, per this campaign's "one submit script per push" convention (see
# CLAUDE.md's cluster conventions section) -- keeps that push's config
# reconstructible from git history.
#
# Runtime estimate (login-node timing: n_refine=32, 10000 meas, 94s ->
# 9.4ms/meas, cost scales ~linearly with n_sites = 10*n_refine^2+2):
# at n_traj=200000/n_skip=2 -> 100,000 meas/shard, n_refine=128
# (163,842 sites, 16x n_refine=32's site count) is the long pole at
# ~150ms/meas * 100,000 ~= 4.2 hours per SHARD (unchanged by N_SHARDS --
# each shard is a full independent run at the same n_traj); all
# comfortably inside the H_RT ceiling below. N_SHARDS independent-seed
# shards per ladder point (per user direction "use more parallelism"),
# mirroring production_shard_task.sh's convention -- more concurrent SGE
# tasks spread across more compute nodes, AND N_SHARDS x more total
# statistics once pooled at analysis time (scripts/plot_lean_ladder.py
# already pools multiple shard files per n_refine automatically).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/lean_harmonic_stats"
test -x "$EXE"

LADDER_CSV="4:6:8:10:12:14:16:20:24:28:32:40:48:56:64:80:96:112:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=16                       # matches submit_production.sh's convention
N_TASKS=$((N_POINTS * N_SHARDS))  # 19 * 16 = 304

L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=200000          # 2x submit_lean_harmonic_stats.sh's 100000, per shard
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
JACK_BLOCK_SIZE=1000    # n_traj/n_skip/jack_block_size = 100 native blocks
SEED_BASE=42             # same base as the earlier lean-stats runs, for continuity
H_RT=${H_RT:-12:00:00}   # generous overnight ceiling; n_refine=128 estimated ~4.2h/shard above

JOB_TAG="${1:-lean_ladder_overnight_$(date +%Y-%m-%d)}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "$JOB_TAG",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "n_shards": $N_SHARDS,
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

echo "submitting $N_TASKS-task lean_harmonic_stats overnight array, project qfe, unpinned queue"
qsub \
  -P qfe \
  -N "s2prec_$JOB_TAG" \
  -j y \
  -o "$LOG_DIR/" \
  -t "1-$N_TASKS" \
  -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",L_MAX="$L_MAX",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",JACK_BLOCK_SIZE="$JACK_BLOCK_SIZE",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
