#!/usr/bin/env bash
set -euo pipefail

# Builds src/ising_s2_crit.cc (the K_c-search / Ising-on-S2 driver this
# campaign uses) for the SCC cluster.
#
# Fully self-contained as of 2026-08-22: source, headers (ising.h, S2.h,
# lattice.h, rng.h, grp_o3.h, statistics.h, timer.h, util.h), vendored Eigen
# (header-only), and the grp/ group-theory data directory all live under
# this package (vendored from ../IsingS2 via cp/rsync, not symlinks or a
# build-time dependency) per user request that everything needed to run
# this project live in this directory. This supersedes the earlier
# "build/extend ../IsingS2 in place, do not fork" convention recorded in
# directive.md/CLAUDE.md — see journal.md 2026-08-22 for the reversal and
# its consequence: this vendored copy no longer tracks ../IsingS2
# automatically. If ../IsingS2 gets a fix this campaign needs, re-copy the
# specific file(s) here manually and note it in journal.md.
#
# Boost (S2.h needs boost/math/special_functions for spherical_harmonic) is
# not vendored, so it's pulled from the module system.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="$ROOT/bin"

module load boost/1.83.0

echo "package_root=$ROOT"
echo "compiler=${CXX:-g++}"
"${CXX:-g++}" --version | head -n 1

mkdir -p "$BIN_DIR"

"${CXX:-g++}" -g -O3 -fopenmp -Wall -Wno-deprecated-declarations -Wno-sign-compare \
  -I "$ROOT/include" -I "$ROOT/include/unsupported" -I "$SCC_BOOST_INCLUDE" \
  -DGRP_DIR="\"$ROOT/grp\"" \
  "$ROOT/src/ising_s2_crit.cc" -o "$BIN_DIR/ising_s2_crit"

test -x "$BIN_DIR/ising_s2_crit"
echo "built=$BIN_DIR/ising_s2_crit"
