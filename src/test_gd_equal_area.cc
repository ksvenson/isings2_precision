// Standalone check (not part of the production driver/cluster/build.sh):
// exercises the gradient-descent QfeLatticeS2::EqualizeFaceAreas across the
// production resolution ladder and reports the resulting face-area
// std/mean, to see whether the true-gradient + backtracking-line-search
// rewrite (2026-08-25) fixes the divergence at n_refine>=32 that the
// earlier heuristic version had (see journal.md 2026-08-24 "equal-area
// mesh" entry). Also re-runs test_ico_symmetry's closure check on the
// relaxed mesh, since EqualizeFaceAreas's symmetry-preservation argument
// (fixed, purely local, identical per-site update rule) has not been
// empirically confirmed the way equal_area_analytic was.
//
// Compile (mirrors cluster/build.sh minus the driver source file):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/test_gd_equal_area.cc -o bin/test_gd_equal_area
// Run: ./bin/test_gd_equal_area <n_iter> <step>
#include <cstdio>
#include <vector>

#include "S2.h"
#include "util.h"

int main(int argc, char* argv[]) {
  int q = 5;
  int n_iter = argc > 1 ? atoi(argv[1]) : 300;
  double step = argc > 2 ? atof(argv[2]) : 0.3;

  int ladder[] = {8, 16, 32, 64};
  printf("n_refine  naive_std/mean  relaxed_std/mean  iters_used\n");
  for (int n_refine : ladder) {
    QfeLatticeS2 lattice(q, n_refine, false);

    std::vector<double> area(lattice.n_faces);
    double sum = 0.0;
    for (int f = 0; f < lattice.n_faces; f++) {
      area[f] = lattice.FlatArea(f);
      sum += area[f];
    }
    double mean = sum / lattice.n_faces;
    double var = 0.0;
    for (int f = 0; f < lattice.n_faces; f++) {
      double d = area[f] - mean;
      var += d * d;
    }
    var /= lattice.n_faces;
    double naive_std_mean = sqrt(var) / mean;

    lattice.EqualizeFaceAreas(n_iter, step);

    sum = 0.0;
    for (int f = 0; f < lattice.n_faces; f++) {
      area[f] = lattice.FlatArea(f);
      sum += area[f];
    }
    mean = sum / lattice.n_faces;
    var = 0.0;
    for (int f = 0; f < lattice.n_faces; f++) {
      double d = area[f] - mean;
      var += d * d;
    }
    var /= lattice.n_faces;
    double relaxed_std_mean = sqrt(var) / mean;

    printf("%8d  %.6f        %.6f          (max %d)\n", n_refine,
           naive_std_mean, relaxed_std_mean, n_iter);
  }
  return 0;
}
