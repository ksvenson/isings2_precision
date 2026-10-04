// Standalone check (not part of cluster/build.sh): verifies that
// QfeLatticeS2::EqualizeCircumPerim preserves icosahedral symmetry, the
// same brute-force closure check test_ico_symmetry.cc uses for
// naive/equal_area_analytic. See IsingS2_precision/journal.md 2026-08-25
// "ROOT CAUSE FOUND" entry -- EqualizeCircumPerim shares
// EqualizeFaceAreas's architecture (fixed, purely local per-vertex update
// rule), so the same symmetry-preservation argument should apply, but
// this is verified directly rather than assumed.
//
// Compile (mirrors test_ico_symmetry.cc):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/test_eqrp_symmetry.cc -o bin/test_eqrp_symmetry
// Run: ./bin/test_eqrp_symmetry <n_refine> <n_iter> <step>
#include <cstdio>
#include <vector>

#include "S2.h"
#include "util.h"

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = argc > 1 ? atoi(argv[1]) : 8;
  int n_iter = argc > 2 ? atoi(argv[2]) : 2000;
  double step = argc > 3 ? atof(argv[3]) : 0.3;

  QfeLatticeS2 lattice(q, n_refine);
  int iters_used = 0;
  lattice.EqualizeCircumPerim(n_iter, step, &iters_used);
  printf("n_refine=%d n_sites=%d n_group_elements=%zu iters_used=%d\n",
         n_refine, lattice.n_sites, lattice.G.size(), iters_used);

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
  printf("RESULT n_refine=%d max_mismatch=%.6e (worst s=%d g=%d)\n",
         n_refine, max_mismatch, worst_s, worst_g);
  return 0;
}
