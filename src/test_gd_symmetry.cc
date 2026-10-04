// Standalone check (not part of the production driver/cluster/build.sh):
// same brute-force icosahedral-closure check as test_ico_symmetry.cc, but
// applied to the mesh produced by QfeLatticeS2::EqualizeFaceAreas (gradient
// descent + backtracking line search), rather than the naive or
// equal_area_analytic constructions. EqualizeFaceAreas's symmetry argument
// (fixed, identical, purely local per-site update rule -> preserves the
// naive mesh's exact icosahedral symmetry) was only ever argued, not
// empirically confirmed -- do that here before trusting equal_face_area
// for production. O(n_sites^2 * n_group), so only cheap at small n_refine.
//
// Compile (mirrors cluster/build.sh minus the driver source file):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/test_gd_symmetry.cc -o bin/test_gd_symmetry
// Run: ./bin/test_gd_symmetry <n_refine> <n_iter> <step>
#include <cstdio>
#include <vector>

#include "S2.h"
#include "util.h"

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = argc > 1 ? atoi(argv[1]) : 8;
  int n_iter = argc > 2 ? atoi(argv[2]) : 300;
  double step = argc > 3 ? atof(argv[3]) : 0.3;

  QfeLatticeS2 lattice(q, n_refine, false);
  lattice.EqualizeFaceAreas(n_iter, step);
  printf("n_refine=%d n_iter=%d step=%.3f n_sites=%d n_group_elements=%zu\n",
         n_refine, n_iter, step, lattice.n_sites, lattice.G.size());

  double max_mismatch = 0.0;
  int worst_s = -1, worst_g = -1;
  for (int s = 0; s < lattice.n_sites; s++) {
    for (size_t g = 0; g < lattice.G.size(); g++) {
      Vec3 gr = lattice.G[g] * lattice.r[s];
      double best = 1e18;
      for (int s2 = 0; s2 < lattice.n_sites; s2++) {
        double d = (lattice.r[s2] - gr).norm();
        if (d < best) best = d;
      }
      if (best > max_mismatch) {
        max_mismatch = best;
        worst_s = s;
        worst_g = (int)g;
      }
    }
  }
  printf("RESULT n_refine=%d max_mismatch=%.6e (worst s=%d g=%d)\n", n_refine,
         max_mismatch, worst_s, worst_g);
  return 0;
}
