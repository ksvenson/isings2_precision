#!/usr/bin/env bash
# naive-vs-equal_area kappa2_r ladder extended one rung past the existing
# n_refine=512 ladder, to n_refine=768. Copied from
# submit_lean_ladder_512.sh per CLAUDE.md's "copy, don't edit in place"
# convention so the config that produced
# campaign_runs/lean_ladder_512_2026-08-27/ stays reconstructible.
#
# Purpose: the objective-1 gate (PLAN.md step 6) still FAILS -- l>=3
# kappa2_r sits on a nonzero, resolution-independent plateau on both
# naive and equal_area out to n_refine=512 (journal.md 2026-09-02). The
# only thing a bigger mesh can settle is whether that plateau *ever*
# starts to shrink. 768 is 2.25x the sites of 512, the largest rung that
# still fits a single SGE task comfortably.
#
# Ladder: {4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768} -- the 512
# ladder's 15 points plus 768. The whole ladder is re-run rather than just
# the new rung, so this push is self-contained: it does not depend on
# having read access to the earlier push's output (which lives in another
# user's project directory). To run ONLY the new rung against existing
# data, set LADDER_SHORT="" LADDER_LONG="768" in the environment.
#
# -------------------------------------------------------------------
# Differences from submit_lean_ladder_512.sh, and why
# -------------------------------------------------------------------
# 1. The production array is SPLIT into a short-rung and a long-rung array
#    instead of one array under a single h_rt. The 512 push used a single
#    h_rt=30:00:00 and lost tasks to walltime kills concentrated entirely
#    at its top rung -- 19 of the eqarea tasks, 3 of which left a shard
#    file with a truncated final jackknife block (journal.md 2026-09-02).
#    Splitting lets the 11 cheap rungs queue under a short h_rt (they are
#    minutes-to-hours) while only the 5 expensive ones pay the scheduling
#    cost of a long one.
# 2. h_rt on the long array is 120:00:00, not 30:00:00. Cost model,
#    calibrated against this campaign's own numbers (~9.4ms/measurement at
#    n_refine=32, scaling linearly in n_sites = 10*n^2+2 -- journal.md
#    2026-08-26, which correctly predicted the 512 rung's ~24h): 768 has
#    5,898,242 sites, ~5.4 s/measurement, x35000 measurements/shard =
#    ~53h. The 512 rung's kills show that estimate runs optimistic, so the
#    margin here is deliberately large rather than tight.
# 3. The long array and the prewarm array request `-pe omp 4`. Both
#    binaries are single-threaded; the slots are purely for the SCC
#    per-core memory allocation, same rationale as
#    submit_gd_stress_1024.sh's `-pe omp 16`. Measured footprint is
#    1.96 kB/site for lean_harmonic_stats and 1.48 kB/site for
#    ising_s2_crit, so n_refine=768 needs ~11.6 GB and ~8.7 GB
#    respectively -- past a default single slot on most SCC nodes, which
#    is why 512 (~5.1 GB) got away without a reservation and 768 will not.
#    If your nodes are memory-lean, `-l mem_per_core=8G` is the
#    alternative knob.
# 4. Seed bases are 920000/930000, disjoint from the 512 push's
#    900000/910000 (+512 tasks each) and from every earlier push.
#
# Everything else -- l_max, coupling rule, n_traj, n_skip, n_wolff,
# n_metropolis, jack_block_size, shard count, mesh-relaxation settings --
# is deliberately IDENTICAL to the 512 push, so the new rung is directly
# comparable to the existing ladder rather than confounded with a
# statistics or parameter change.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Split ladder: cheap rungs vs expensive rungs. Overridable from the
# environment to re-run a subset (e.g. LADDER_SHORT="" LADDER_LONG="768").
LADDER_SHORT="${LADDER_SHORT-4:6:8:12:16:24:32:48:64:96:128}"
LADDER_LONG="${LADDER_LONG-192:256:384:512:768}"
LADDER_ALL="$(echo "${LADDER_SHORT:+$LADDER_SHORT:}${LADDER_LONG}" | sed 's/^://; s/:$//')"

IFS=':' read -r -a LADDER <<<"$LADDER_ALL"
N_POINTS=${#LADDER[@]}

N_SHARDS=32
L_MAX=8
COUPLING_RULE=exact_sinh
N_THERM=2000
N_TRAJ=70000
N_SKIP=2
N_WOLFF=5
N_METROPOLIS=4
# n_traj/n_skip = 35000 measurements/shard; /1000 = 35 native blocks/shard,
# x32 shards = 1120 pooled blocks/point -- comfortably past PLAN.md step
# 5's gate of 10*(2*l_max+1) = 170 at l_max=8. Same value as the 512 push
# (directive.md requires recomputing this per push, not copying it; it is
# recomputed here and lands on the same number because n_traj is unchanged).
JACK_BLOCK_SIZE=1000
EQUAL_AREA_ITERS=20000
EQUAL_AREA_STEP=0.3

H_RT_SHORT=12:00:00
H_RT_LONG=120:00:00
H_RT_PREWARM=24:00:00
PE_LONG=4

# Output directory. A NEW run is stamped with today's date; a SKIP_PREWARM
# recovery run instead adopts the most recent existing run directory.
#
# That asymmetry is deliberate. This ladder's prewarm runs for hours, so a
# recovery is typically attempted the NEXT day -- and a plain date stamp
# would then invent a fresh directory, find an empty mesh_cache, and
# refuse, for no reason other than the clock. Resolving the latest
# existing directory instead means neither path ever needs a filepath
# passed in by hand. OUT_ROOT_BASE still overrides both, for the case of
# several runs on disk and a specific older one to resume.
if [[ -n "${OUT_ROOT_BASE:-}" ]]; then
  : # caller chose explicitly; respect it
elif [[ "${SKIP_PREWARM:-0}" == "1" ]]; then
  # newest first; -maxdepth/-mindepth keep this to the run dirs themselves
  OUT_ROOT_BASE=$(find "$ROOT/campaign_runs" -mindepth 1 -maxdepth 1 -type d \
                    -name 'lean_ladder_768_*' -printf '%T@ %p\n' 2>/dev/null \
                  | sort -rn | head -1 | cut -d' ' -f2-)
  if [[ -z "$OUT_ROOT_BASE" ]]; then
    echo "ERROR: SKIP_PREWARM=1 but no campaign_runs/lean_ladder_768_* dir exists." >&2
    echo "       Nothing to resume -- run without SKIP_PREWARM to start fresh." >&2
    exit 1
  fi
  echo "SKIP_PREWARM: resuming most recent run dir $OUT_ROOT_BASE"
else
  OUT_ROOT_BASE="$ROOT/campaign_runs/lean_ladder_768_$(date +%F)"
fi
OUT_ROOT_NAIVE="$OUT_ROOT_BASE/naive"
OUT_ROOT_EQAREA="$OUT_ROOT_BASE/eqarea"
MESH_CACHE_DIR="$OUT_ROOT_BASE/mesh_cache"
# SGE's -o log directory must exist before qsub runs or every task fails
# immediately (Eqw, "can't open output file ... Is a directory").
mkdir -p "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA" "$MESH_CACHE_DIR" \
         "${OUT_ROOT_BASE}_naive_logs" "${OUT_ROOT_BASE}_eqarea_logs"

SIM_SHA=$(sha256sum "$ROOT/bin/lean_harmonic_stats" | cut -d' ' -f1)
# A recovery re-run must not overwrite the original manifest -- that file is
# this push's provenance record, and the rerun's timestamp/sha would silently
# replace the ones that actually produced the prewarmed meshes.
MANIFEST="$OUT_ROOT_BASE/manifest.json"
if [[ -f "$MANIFEST" ]]; then
  echo "keeping existing $MANIFEST (not overwriting on a re-run)"
  MANIFEST=/dev/null
fi
cat > "$MANIFEST" <<EOF
{
  "campaign_id": "ising_s2_precision_lean_harmonic_stats",
  "stage": "$(basename "$OUT_ROOT_BASE")",
  "ladder": [$(IFS=,; echo "${LADDER[*]}")],
  "ladder_short": "$LADDER_SHORT", "ladder_long": "$LADDER_LONG",
  "n_shards": $N_SHARDS, "l_max": $L_MAX, "coupling_rule": "$COUPLING_RULE",
  "n_therm": $N_THERM, "n_traj": $N_TRAJ, "n_skip": $N_SKIP,
  "n_wolff": $N_WOLFF, "n_metropolis": $N_METROPOLIS,
  "jack_block_size": $JACK_BLOCK_SIZE,
  "mesh_modes_compared": ["naive", "equal_area"],
  "equal_area_iters": $EQUAL_AREA_ITERS, "equal_area_step": $EQUAL_AREA_STEP,
  "seed_base_naive": 920000, "seed_base_eqarea": 930000,
  "h_rt_short": "$H_RT_SHORT", "h_rt_long": "$H_RT_LONG", "pe_long": $PE_LONG,
  "simulator_sha256": "$SIM_SHA",
  "note": "extends lean_ladder_512_2026-08-27 by one rung to n_refine=768; identical stats/params so the new rung is directly comparable"
}
EOF
[[ "$MANIFEST" == /dev/null ]] || echo "wrote $MANIFEST"

test -x "$ROOT/bin/lean_harmonic_stats" || { echo "ERROR: build bin/lean_harmonic_stats first" >&2; exit 1; }
test -x "$ROOT/bin/ising_s2_crit"       || { echo "ERROR: run cluster/build.sh first" >&2; exit 1; }

# ---------------------------------------------------------------
# step 1: prewarm the equal_area mesh cache, serially per point
# ---------------------------------------------------------------
# Without this, the 32 shards of a given (n_refine, equal_area) point all
# relax the same mesh from cold at once and race to write the same cache
# file. `-sync y` is used here deliberately: CLAUDE.md's "prefer detached
# qsub" rule is about *production* arrays, and the whole point of this
# step is that production must not start until it lands. Per that same
# convention, the -sync y return is NOT trusted on its own -- the explicit
# file check below is the real gate.
# RECOVERY PATH -- set SKIP_PREWARM=1 to jump straight to step 2.
#
# `-sync y` blocks this script until the whole prewarm array finishes,
# which at this ladder's top rung is HOURS (the n_refine=768 relaxation is
# 144x the work of n_refine=64; measured 2026-10-04). If the login session
# running this script dies in the meantime, the prewarm array survives --
# it belongs to SGE, not to your shell -- but the production arrays never
# get submitted, because the script that submits them is gone.
#
# Do NOT just re-run this script in that situation: it would submit a
# SECOND prewarm array racing the first one on the same cache files, which
# is the exact race the prewarm exists to prevent. Instead, once the
# original prewarm array has finished, re-run with:
#
#     SKIP_PREWARM=1 bash cluster/sge/submit_lean_ladder_768.sh
#
# which skips the qsub but still runs the cache verification below -- that
# check, not the -sync return, was always the real gate.
if [[ "${SKIP_PREWARM:-0}" == "1" ]]; then
  echo "=== step 1: SKIPPED (SKIP_PREWARM=1) -- verifying existing cache ==="
else
  echo "=== step 1: prewarm equal_area mesh cache ($N_POINTS points) ==="
  echo "    NOTE: this blocks for hours at this ladder's top rung. If this"
  echo "    shell may not survive that, see SKIP_PREWARM in this script."
  PREWARM_JOB=$(qsub -terse \
    -P qfe -N s2prec_lean_ladder_768_prewarm -j y -o "$OUT_ROOT_BASE/" \
    -sync y -t "1-$N_POINTS" -l h_rt=$H_RT_PREWARM -pe omp $PE_LONG \
    -v ROOT="$ROOT",LADDER_CSV="$LADDER_ALL",L_MAX="$L_MAX",EQUAL_AREA_ITERS="$EQUAL_AREA_ITERS",EQUAL_AREA_STEP="$EQUAL_AREA_STEP",MESH_CACHE_DIR="$MESH_CACHE_DIR",MESH_MODE=equal_area \
    "$ROOT/cluster/sge/prewarm_mesh_cache_task.sh")
  echo "prewarm returned: $PREWARM_JOB"
fi

# Verify EXACT file size, not just non-emptiness. S2.h's WritePositions
# writes an int32 n_sites header followed by n_sites * Vec3 (3 doubles),
# so a complete cache for q=5 is exactly 4 + 24*(10*n^2+2) = 240*n^2 + 52
# bytes (verified against real cache files at n_refine=4/8/16).
#
# This matters because a prewarm task killed mid-fwrite -- by a walltime
# limit, or by a qdel -- leaves a TRUNCATED cache file. A `-s` test passes
# it, but ReadPositions correctly rejects it at runtime as a cache miss,
# at which point all 32 shards of that ladder point relax the mesh from
# cold simultaneously and race to rewrite the same file. That is precisely
# the failure the prewarm stage exists to prevent, so it must be caught
# here rather than discovered in production.
STEP_FMT=$(printf '%.3f' "$EQUAL_AREA_STEP")
MISSING=0
for n in "${LADDER[@]}"; do
  f="$MESH_CACHE_DIR/q5k${n}_eqarea_step${STEP_FMT}.dat"
  want=$((240 * n * n + 52))
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing mesh cache for n_refine=$n ($f)" >&2
    MISSING=1
  else
    got=$(stat -c%s "$f" 2>/dev/null || echo 0)
    if [[ "$got" != "$want" ]]; then
      echo "ERROR: TRUNCATED mesh cache for n_refine=$n: $got bytes, expected $want" >&2
      echo "       ($f -- delete it and re-run the prewarm for this point)" >&2
      MISSING=1
    fi
  fi
done
if (( MISSING )); then
  echo "ERROR: prewarm incomplete -- NOT submitting production." >&2
  echo "       Inspect $OUT_ROOT_BASE/s2prec_lean_ladder_768_prewarm.o* ." >&2
  echo "       If the prewarm array is still running, wait for it and then use" >&2
  echo "       SKIP_PREWARM=1 bash \$0  -- do NOT re-run this script as-is," >&2
  echo "       which would submit a second prewarm array racing the first." >&2
  exit 1
fi
echo "all $N_POINTS cache files verified present and byte-complete"

# ---------------------------------------------------------------
# step 2: production arrays (4 total: {short,long} x {naive,eqarea})
# ---------------------------------------------------------------
submit_arm() {
  local name=$1 ladder_csv=$2 out_root=$3 logdir=$4 h_rt=$5 pe=$6 seed_base=$7 mesh_mode=$8
  [[ -z "$ladder_csv" ]] && { echo "skipping $name (empty ladder)"; return; }
  local n_points n_tasks
  IFS=':' read -r -a arr <<<"$ladder_csv"
  n_points=${#arr[@]}
  n_tasks=$((n_points * N_SHARDS))

  local pe_args=()
  [[ -n "$pe" ]] && pe_args=(-pe omp "$pe")

  local vars="ROOT=$ROOT,LADDER_CSV=$ladder_csv,N_SHARDS=$N_SHARDS,OUT_ROOT=$out_root"
  vars+=",L_MAX=$L_MAX,COUPLING_RULE=$COUPLING_RULE,N_THERM=$N_THERM,N_TRAJ=$N_TRAJ"
  vars+=",N_SKIP=$N_SKIP,N_WOLFF=$N_WOLFF,N_METROPOLIS=$N_METROPOLIS"
  vars+=",JACK_BLOCK_SIZE=$JACK_BLOCK_SIZE,SEED_BASE=$seed_base,MESH_MODE=$mesh_mode"
  if [[ "$mesh_mode" != naive ]]; then
    vars+=",EQUAL_AREA_ITERS=$EQUAL_AREA_ITERS,EQUAL_AREA_STEP=$EQUAL_AREA_STEP,MESH_CACHE_DIR=$MESH_CACHE_DIR"
  fi

  echo "submitting $name: $n_points points x $N_SHARDS shards = $n_tasks tasks, h_rt=$h_rt ${pe:+-pe omp $pe}"
  qsub -P qfe -N "$name" -j y -o "$logdir/" \
    -t "1-$n_tasks" -l h_rt="$h_rt" "${pe_args[@]}" \
    -v "$vars" \
    "$ROOT/cluster/sge/lean_harmonic_stats_task.sh"
}

echo "=== step 2: submit production arrays ==="
# Seed bases are offset per arm so that no two arms can collide even
# though each arm restarts its task index at 1. Each arm uses at most
# 11*32 = 352 seeds, so 1000 of headroom per arm is ample.
submit_arm s2prec_ladder768_naive_short  "$LADDER_SHORT" "$OUT_ROOT_NAIVE"  "${OUT_ROOT_BASE}_naive_logs"  "$H_RT_SHORT" ""        920000 naive
submit_arm s2prec_ladder768_naive_long   "$LADDER_LONG"  "$OUT_ROOT_NAIVE"  "${OUT_ROOT_BASE}_naive_logs"  "$H_RT_LONG"  "$PE_LONG" 921000 naive
submit_arm s2prec_ladder768_eqarea_short "$LADDER_SHORT" "$OUT_ROOT_EQAREA" "${OUT_ROOT_BASE}_eqarea_logs" "$H_RT_SHORT" ""        930000 equal_area
submit_arm s2prec_ladder768_eqarea_long  "$LADDER_LONG"  "$OUT_ROOT_EQAREA" "${OUT_ROOT_BASE}_eqarea_logs" "$H_RT_LONG"  "$PE_LONG" 931000 equal_area

cat <<EOF

submitted (detached) -- track with:
  qstat -u \$USER | grep s2prec_ladder768
  qacct -j <jobid> | grep -c 'exit_status  0'

once complete:
  # 1. sanity-check shard counts (expect $N_SHARDS files per rung per mode)
  for d in "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA"; do
    echo "\$d"; ls "\$d" | sed 's/_lean_.*//' | sort | uniq -c
  done
  # 2. the target figure
  .venv_plot/bin/python scripts/plot_mesh_compare_full.py \\
      "$OUT_ROOT_NAIVE" "$OUT_ROOT_EQAREA" --l_max $L_MAX \\
      -o "$OUT_ROOT_BASE/plots"
  # -> $OUT_ROOT_BASE/plots/mesh_compare_kappa2_r_vs_invL.png
  # 3. the per-mode numeric tables behind it
  .venv_plot/bin/python scripts/symmetry_check_lean_kappa2.py "$OUT_ROOT_NAIVE"  --l_max $L_MAX
  .venv_plot/bin/python scripts/symmetry_check_lean_kappa2.py "$OUT_ROOT_EQAREA" --l_max $L_MAX
EOF
