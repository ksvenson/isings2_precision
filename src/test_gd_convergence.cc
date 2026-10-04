// Standalone check (not part of the production driver/cluster/build.sh):
// logs QfeLatticeS2::EqualizeFaceAreas's area std/mean at geometrically
// spaced iteration checkpoints, for one n_refine, to characterize the
// convergence rate (not just the endpoint) -- needed to fit an
// iterations-needed(L) model for the mesh-position cache/library (see
// journal.md 2026-08-25 "mesh position cache" entry). EqualizeFaceAreas is
// called repeatedly with increasing checkpoint deltas rather than in one
// shot, resuming from the current (already-relaxed) site positions each
// time -- note this resets the adaptive backtracking step back to the
// caller-supplied `step` on every call (a minor inefficiency, not a
// correctness issue: the objective/gradient only depend on current
// positions).
//
// Compile (mirrors cluster/build.sh minus the driver source file):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/test_gd_convergence.cc -o bin/test_gd_convergence
// Run: ./bin/test_gd_convergence <n_refine> <step> <max_iter>
#include <cstdio>
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
  int n_refine = argc > 1 ? atoi(argv[1]) : 64;
  double step = argc > 2 ? atof(argv[2]) : 0.3;
  int max_iter = argc > 3 ? atoi(argv[3]) : 4000;

  QfeLatticeS2 lattice(q, n_refine);
  printf("n_refine=%d n_sites=%d n_faces=%d\n", n_refine, lattice.n_sites,
         lattice.n_faces);
  printf("iter std_mean\n");
  printf("%d %.8f\n", 0, AreaStdMean(lattice));
  fflush(stdout);

  int done = 0;
  int checkpoint = 10;
  while (done < max_iter) {
    int step_iters = std::min(checkpoint - done, max_iter - done);
    if (step_iters <= 0) break;
    lattice.EqualizeFaceAreas(step_iters, step);
    done += step_iters;
    printf("%d %.8f\n", done, AreaStdMean(lattice));
    fflush(stdout);
    checkpoint = int(checkpoint * 1.5) + 1;
  }
  return 0;
}
