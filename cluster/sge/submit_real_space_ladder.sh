#!/usr/bin/env bash
# First production-scale resolution ladder for the direct real-space
# Delta_sigma extraction (real_space_2pt_test/fit_real_space_2pt.py),
# per PLAN.md's "next step" after the 2026-08-25 breakthrough: continuum-
# extrapolate Delta(a) from real-space fits, rather than resting on the
# two single-seed login-node spot checks (n_refine=16,32) done so far.
#
# **Config updated 2026-08-26** for real_space_2pt_test.cc's rewrite that
# dropped z=cos(theta) histogram binning (and the --orbit_avg flag it later
# grew) for exact icosahedral-symmetry equivalence classes -- see that
# file's header comment. N_BINS/ORBIT_AVG are gone; THETA_MIN/THETA_MAX
# (still 0.3/2.8, matching the previously-validated fit window) now just
# bound which pairs get built into classes at all, for compute-cost
# control, not which pairs get pooled together (pooling is now always
# symmetry-exact). This script has not yet been re-run against the
# rewritten driver -- the two production ladders described below
# (7311120/7311202) predate the rewrite and used the old binned format.
#
# Design:
# - naive mesh_mode, exact_sinh coupling_rule -- exactly the configuration
#   validated by the two spot checks (journal.md 2026-08-25 "BREAKTHROUGH"
#   entry: Delta=0.128(0) at n_refine=32, Delta=0.1287(3) at n_refine=16).
#   Not a mesh-mode/coupling-rule comparison -- that question was already
#   settled for the l-space method (all combinations gave the same wrong
#   plateau) and isn't re-opened here.
# - Ladder n_refine={8,16,24,32,48,64,96,128} (8 points) -- covers the two
#   already-validated points (16,32) plus enough extra points on both
#   sides for a real continuum (1/n_refine-power) extrapolation of
#   Delta(a). Not extended down to n_refine=2-4 (push7's low end): those
#   are only a few hundred sites and the real-space fit already excludes
#   the smallest angular bins (theta_min=0.3), which would leave very few
#   usable bins at the coarsest resolutions.
# - N_SHARDS=24 independent-seed shards per ladder point. Unlike
#   ising_s2_crit's in-C++ jack_block_size accumulator, real_space_2pt_test
#   has no jackknife machinery at all -- each shard here IS the jackknife
#   unit (a fully independent MC chain from thermalization), pooled via a
#   delete-one-shard jackknife in scripts/fit_real_space_ladder.py. This is
#   arguably a cleaner decorrelation than same-chain jack_block_size blocks
#   since shards never share thermalization history.
# - n_therm/n_traj/n_skip/n_wolff/n_metropolis/n_bins/n_ref all match the
#   validated spot-check configuration exactly (journal.md 2026-08-25):
#   n_therm=2000, n_traj=20000, n_skip=2 (10000 measured configs/shard),
#   n_wolff=5, n_metropolis=4, n_bins=100, n_ref=300.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXE="$ROOT/bin/real_space_2pt_test"
test -x "$EXE"

LADDER_CSV="8:16:24:32:48:64:96:128"
IFS=':' read -r -a LADDER <<<"$LADDER_CSV"
N_POINTS=${#LADDER[@]}
N_SHARDS=24
N_TASKS=$((N_POINTS * N_SHARDS))

COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=20000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
N_REF=300
THETA_MIN=0.3
THETA_MAX=2.8
SEED_BASE=800000  # disjoint from every prior push's seed range (700000+ used by push7)
H_RT=${H_RT:-2:00:00}

JOB_TAG="${1:-real_space_ladder_$(date +%Y-%m-%d)}"
OUT_ROOT="$ROOT/campaign_runs/$JOB_TAG"
mkdir -p "$OUT_ROOT"
LOG_DIR="$ROOT/campaign_runs/${JOB_TAG}_logs"
mkdir -p "$LOG_DIR"

SIM_CHECKSUM=$(sha256sum "$EXE" | awk '{print $1}')
MANIFEST_PATH="$OUT_ROOT/manifest.json"
cat >"$MANIFEST_PATH" <<EOF
{
  "campaign_id": "ising_s2_precision_real_space_delta_extraction",
  "stage": "$JOB_TAG",
  "purpose": "continuum-extrapolate Delta_sigma(a) via direct real-space two-point-function fits across a resolution ladder, jackknifed over independent shards",
  "ladder": [$(IFS=,; echo "${LADDER_CSV//:/,}")],
  "mesh_mode": "naive",
  "coupling_rule": "$COUPLING_RULE",
  "n_shards_per_point": $N_SHARDS,
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "n_ref": $N_REF, "theta_min": $THETA_MIN, "theta_max": $THETA_MAX,
  "seed_base": $SEED_BASE,
  "simulator_checksum": "$SIM_CHECKSUM"
}
EOF
echo "manifest written to $MANIFEST_PATH"

for k in "${LADDER[@]}"; do
  for s in $(seq -w 0 $((N_SHARDS - 1))); do
    mkdir -p "$OUT_ROOT/q5k${k}/shard_${s}"
  done
done

echo "submitting real-space ladder array: $N_TASKS tasks ($N_POINTS ladder points x $N_SHARDS shards)"
qsub \
  -P qfe -N "s2prec_${JOB_TAG}" -j y -o "$LOG_DIR/" \
  -t "1-$N_TASKS" -l h_rt="$H_RT" \
  -v ROOT="$ROOT",LADDER_CSV="$LADDER_CSV",N_SHARDS="$N_SHARDS",OUT_ROOT="$OUT_ROOT",COUPLING_RULE="$COUPLING_RULE",N_THERM="$N_THERM",N_TRAJ="$N_TRAJ",N_SKIP="$N_SKIP",N_WOLFF="$N_WOLFF",N_METROPOLIS="$N_METROPOLIS",N_REF="$N_REF",THETA_MIN="$THETA_MIN",THETA_MAX="$THETA_MAX",SEED_BASE="$SEED_BASE" \
  "$ROOT/cluster/sge/real_space_shard_task.sh"

echo "submitted (detached) -- track with: qstat -u \$USER | grep s2prec_$JOB_TAG"
