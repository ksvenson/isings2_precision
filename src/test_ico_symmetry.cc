// Standalone check (not part of the production driver/cluster/build.sh):
// verifies that a QfeLatticeS2 mesh (naive or equal_area_analytic) is
// closed under the icosahedral point group -- i.e. for every symmetry
// element g in G and every site s, G[g]*r[s] lands exactly on some other
// site's position. This is the direct empirical test used in journal.md
// 2026-08-24 to find/fix two bugs in ApplyEqualAreaAnalytic.
//
// Compile (mirrors cluster/build.sh minus the driver source file):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/test_ico_symmetry.cc -o bin/test_ico_symmetry
// Run: ./bin/test_ico_symmetry <n_refine> <0|1 for naive|analytic>
#include <cstdio>
#include <vector>

#include "S2.h"
#include "util.h"

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = argc > 1 ? atoi(argv[1]) : 8;
  bool analytic = argc > 2 ? atoi(argv[2]) != 0 : true;

  QfeLatticeS2 lattice(q, n_refine, analytic);
  printf("n_refine=%d analytic=%d n_sites=%d n_group_elements=%zu\n",
         n_refine, analytic, lattice.n_sites, lattice.G.size());

  // build a lookup: bucket sites by a coarse hash of their coordinates,
  // for fast nearest-match search (n_sites can be tens of thousands).
  // brute-force nearest neighbor is fine at these sizes for a one-off
  // check.
  double max_mismatch = 0.0;
  int worst_s = -1, worst_g = -1;
  for (int s = 0; s < lattice.n_sites; s++) {
    for (size_t g = 0; g < lattice.G.size(); g++) {
      Vec3 gr = lattice.G[g] * lattice.r[s];
      // find nearest actual site to gr
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
    if (s % 20 == 0) {
      printf("  checked site %d/%d, max_mismatch_so_far=%.3e\n", s,
             lattice.n_sites, max_mismatch);
    }
  }
  printf("RESULT n_refine=%d analytic=%d max_mismatch=%.6e (worst s=%d g=%d)\n",
         n_refine, analytic, max_mismatch, worst_s, worst_g);
  return 0;
}
