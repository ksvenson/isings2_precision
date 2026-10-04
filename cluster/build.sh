#!/usr/bin/env bash
set -euo pipefail

# Builds this campaign's C++ drivers into bin/.
#
#   cluster/build.sh                      # the production drivers (default)
#   cluster/build.sh --all                # production drivers + diagnostics
#   cluster/build.sh lean_harmonic_stats  # just the named target(s)
#   cluster/build.sh --list               # show what each group contains
#
# Every target is a SINGLE self-contained translation unit: one src/*.cc
# that #includes the headers it needs. There are no object files to link,
# no static/shared library to build, and no Makefile dependency graph --
# each build is one g++ invocation with the same flags. (Never compile
# anything in include/ directly; headers are pulled in by the .cc files.)
#
# Fully self-contained as of 2026-08-22: source, headers (ising.h, S2.h,
# lattice.h, rng.h, grp_o3.h, statistics.h, timer.h, util.h), vendored
# Eigen (header-only), and the grp/ group-theory data directory all live
# under this package -- see journal.md 2026-08-22 for the reversal of the
# earlier "build against ../IsingS2 in place" rule, and its consequence:
# this vendored copy does NOT track ../IsingS2 automatically.
#
# Boost (S2.h needs boost/math/special_functions/spherical_harmonic) is
# the one dependency that is not vendored. On the SCC it comes from the
# module system; off-cluster this script falls back to a system boost
# (e.g. Ubuntu's libboost-math-dev), so the same script works in both
# places -- verified 2026-10-04 on a plain Linux box with all 16 targets.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="$ROOT/bin"

# Drivers that produce production data or are referenced by a cluster task
# script. These are what you need before submitting any job.
PRODUCTION=(
  ising_s2_crit              # measurement driver; also used for mesh-cache prewarm
  lean_harmonic_stats        # the kappa2_r / Delta_s workhorse
  real_space_2pt_test        # real-space Delta_sigma extraction
  save_configs               # raw-config dump for offline M_l reconstruction
  analyze_north_pole_configs # offline single-row <s_i s_j> from those configs
  full_corr_test             # dense n_sites^2 correlator (small n_refine only)
)

# One-off correctness/performance checks. Not production paths; safe to
# build and safe to re-run (all are read-only diagnostics), just rarely
# needed. See each file's header comment for what it checks.
DIAGNOSTICS=(
  test_ico_symmetry test_gd_symmetry test_gd_equal_area test_gd_convergence
  test_gd_stress test_eqrp_symmetry test_orbit_avg fem_scalar_test
  diag_optimize_integrator dump_face_sa_range
)

if [[ "${1:-}" == "--list" ]]; then
  printf 'production:  %s\n' "${PRODUCTION[*]}"
  printf 'diagnostics: %s\n' "${DIAGNOSTICS[*]}"
  exit 0
fi

case "${1:-}" in
  --all) TARGETS=("${PRODUCTION[@]}" "${DIAGNOSTICS[@]}") ;;
  "")    TARGETS=("${PRODUCTION[@]}") ;;
  *)     TARGETS=("$@") ;;
esac

# Boost: prefer the SCC module when the module system is present, else
# assume a system-wide boost on the default include path. `module` is
# usually a shell function rather than a binary, so `type` is the check
# that works for both.
BOOST_INC=""
if type module >/dev/null 2>&1; then
  module load boost/1.83.0
  BOOST_INC="${SCC_BOOST_INCLUDE:-}"
fi

echo "package_root=$ROOT"
echo "compiler=${CXX:-g++}"
"${CXX:-g++}" --version | head -n 1
if [[ -n "$BOOST_INC" ]]; then
  echo "boost_include=$BOOST_INC (module)"
else
  echo "boost_include=<system default>"
fi

mkdir -p "$BIN_DIR"

CXXFLAGS=(-g -O3 -fopenmp -Wall -Wno-deprecated-declarations -Wno-sign-compare)
INCLUDES=(-I "$ROOT/include" -I "$ROOT/include/unsupported")
[[ -n "$BOOST_INC" ]] && INCLUDES+=(-I "$BOOST_INC")

failed=()
for t in "${TARGETS[@]}"; do
  src="$ROOT/src/$t.cc"
  if [[ ! -f "$src" ]]; then
    echo "ERROR: no such target: $t (expected $src)" >&2
    failed+=("$t")
    continue
  fi
  printf '  building %-28s ' "$t"
  if "${CXX:-g++}" "${CXXFLAGS[@]}" "${INCLUDES[@]}" \
       -DGRP_DIR="\"$ROOT/grp\"" "$src" -o "$BIN_DIR/$t" 2>"$BIN_DIR/.$t.buildlog"; then
    echo "ok"
    rm -f "$BIN_DIR/.$t.buildlog"
  else
    echo "FAILED (see $BIN_DIR/.$t.buildlog)"
    failed+=("$t")
  fi
done

if ((${#failed[@]})); then
  echo "build FAILED for: ${failed[*]}" >&2
  exit 1
fi
echo "built ${#TARGETS[@]} target(s) into $BIN_DIR"
