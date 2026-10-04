// Standalone check (not part of the production driver/cluster/build.sh):
// stress-test QfeLatticeS2::EqualizeFaceAreas at large n_refine (up to
// 2048, ~4.2e7 sites / ~8.4e7 faces) -- construction time, relaxation
// time/iteration, and area-uniformity outcome, to see whether the
// gradient-descent + backtracking-line-search rewrite (verified correct
// and convergent up to n_refine=64, 2026-08-25) still behaves well two
// orders of magnitude further out, or breaks down (timing blowup, memory,
// numerical issues, backtracking never accepting a step, etc).
//
// Compile (mirrors cluster/build.sh minus the driver source file):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/test_gd_stress.cc -o bin/test_gd_stress
// Run: ./bin/test_gd_stress <n_refine> <n_iter> <step>
#include <cstdio>
#include <chrono>
#include <vector>

#include "S2.h"
#include "util.h"

static double AreaStdMean(QfeLatticeS2& lattice) {
  double sum = 0.0, sq_sum = 0.0;
  for (int f = 0; f < lattice.n_faces; f++) {
    double a = lattice.FlatArea(f);
    sum += a;
    sq_sum += a * a;
  }
  double mean = sum / lattice.n_faces;
  double var = sq_sum / lattice.n_faces - mean * mean;
  return sqrt(var) / mean;
}

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = argc > 1 ? atoi(argv[1]) : 2048;
  int n_iter = argc > 2 ? atoi(argv[2]) : 100;
  double step = argc > 3 ? atof(argv[3]) : 0.3;

  fflush(stdout);
  auto t0 = std::chrono::steady_clock::now();
  QfeLatticeS2 lattice(q, n_refine);
  auto t1 = std::chrono::steady_clock::now();
  double construct_sec = std::chrono::duration<double>(t1 - t0).count();
  printf("n_refine=%d n_sites=%d n_faces=%d construct_time=%.3fs\n", n_refine,
         lattice.n_sites, lattice.n_faces, construct_sec);
  fflush(stdout);

  double naive_std_mean = AreaStdMean(lattice);
  printf("naive_std_mean=%.6f\n", naive_std_mean);
  fflush(stdout);

  int iters_used = 0;
  t0 = std::chrono::steady_clock::now();
  lattice.EqualizeFaceAreas(n_iter, step, &iters_used);
  t1 = std::chrono::steady_clock::now();
  double relax_sec = std::chrono::duration<double>(t1 - t0).count();

  double relaxed_std_mean = AreaStdMean(lattice);
  printf("relax_time=%.3fs (%.4fs/iter avg over %d iters actually run)\n",
         relax_sec, relax_sec / std::max(iters_used, 1), iters_used);
  printf(
      "RESULT n_refine=%d n_faces=%d naive_std_mean=%.6f "
      "relaxed_std_mean=%.6f iters_used=%d converged=%d\n",
      n_refine, lattice.n_faces, naive_std_mean, relaxed_std_mean, iters_used,
      iters_used < n_iter ? 1 : 0);
  return 0;
}
