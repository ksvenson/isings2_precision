// Standalone check (not part of the production driver/cluster/build.sh):
// for each cached equal_rp (--mesh_mode equal_rp, EqualizeCircumPerim)
// mesh in mesh_cache/, computes every face's three flat-triangle interior
// angles (law of cosines on the same chord-length convention FlatArea/
// Circumradius/Perimeter already use), converts to Twist's (s,a)
// angle-deviation invariants (d = theta - 60deg; for a given face
// direction i, s_i = sum of the other two d's = -d_i, a_i = their
// difference), and reports the global min/max of s and a needed to cover
// every triangle actually produced by the equal-area-optimized S2 mesh --
// i.e. the (s,a) domain Twist_xR's/Twist's kappa(s,a) map would need to
// span to be usable as the second-pass empirical coupling rule here (see
// PLAN.md "Local coupling assignment (imported from Twist/)").
//
// Compile (mirrors cluster/build.sh minus the driver source file):
//   module load boost/1.83.0
//   g++ -g -O3 -Wno-deprecated-declarations -Wno-sign-compare \
//     -I include -I include/unsupported -I "$SCC_BOOST_INCLUDE" \
//     -DGRP_DIR="\"$(pwd)/grp\"" \
//     src/dump_face_sa_range.cc -o bin/dump_face_sa_range
// Run: ./bin/dump_face_sa_range <mesh_cache_dir>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "S2.h"

int main(int argc, char** argv) {
  std::string cache_dir = (argc > 1) ? argv[1] : "mesh_cache";
  std::vector<int> ks = {2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128};

  double global_s_min = 1e300, global_s_max = -1e300;
  double global_a_min = 1e300, global_a_max = -1e300;

  for (int k : ks) {
    char path[512];
    snprintf(path, sizeof(path), "%s/q5k%d_eqrp_step0.300.dat", cache_dir.c_str(), k);
    FILE* f = fopen(path, "rb");
    if (f == nullptr) {
      printf("k=%-4d SKIP (no cached eqrp mesh at %s)\n", k, path);
      continue;
    }

    QfeLatticeS2 lattice(5, k);
    bool ok = lattice.ReadPositions(f);
    fclose(f);
    if (!ok) {
      printf("k=%-4d SKIP (cache read failed / n_sites mismatch)\n", k);
      continue;
    }

    // vertex degree, to identify the 12 original icosahedron vertices
    // (degree 5) vs. the generic degree-6 refined-mesh vertices
    std::vector<int> degree(lattice.n_sites, 0);
    for (int l = 0; l < lattice.n_links; l++) {
      degree[lattice.links[l].sites[0]]++;
      degree[lattice.links[l].sites[1]]++;
    }

    double s_min = 1e300, s_max = -1e300;
    double a_min = 1e300, a_max = -1e300;
    double max_abs_d = 0.0;
    int s_max_face = -1, a_max_face = -1, s_min_face = -1;
    double s_max_deg[3] = {0, 0, 0}, a_max_deg[3] = {0, 0, 0}, s_min_deg[3] = {0, 0, 0};
    double s_min_angles_deg[3] = {0, 0, 0};

    for (int face = 0; face < lattice.n_faces; face++) {
      Vec3 v0 = lattice.r[lattice.faces[face].sites[0]];
      Vec3 v1 = lattice.r[lattice.faces[face].sites[1]];
      Vec3 v2 = lattice.r[lattice.faces[face].sites[2]];

      double a_len = (v1 - v2).norm();  // side opposite v0
      double b_len = (v2 - v0).norm();  // side opposite v1
      double c_len = (v0 - v1).norm();  // side opposite v2

      auto clamp_unit = [](double x) { return x < -1.0 ? -1.0 : (x > 1.0 ? 1.0 : x); };
      // interior angle at each vertex via law of cosines (flat/chord
      // triangle, same convention as FlatArea/Circumradius/Perimeter)
      double ang0 = acos(clamp_unit((b_len * b_len + c_len * c_len - a_len * a_len) / (2 * b_len * c_len)));
      double ang1 = acos(clamp_unit((a_len * a_len + c_len * c_len - b_len * b_len) / (2 * a_len * c_len)));
      double ang2 = acos(clamp_unit((a_len * a_len + b_len * b_len - c_len * c_len) / (2 * a_len * b_len)));

      double d0 = ang0 - M_PI / 3.0;
      double d1 = ang1 - M_PI / 3.0;
      double d2 = ang2 - M_PI / 3.0;

      // direction i's (s_i,a_i) uses the OTHER two vertex-angle deviations
      // -- s_i = sum (= -d_i since d0+d1+d2=0 for a Euclidean triangle),
      // a_i = difference -- matching Twist_xR's s_v=d_u+d_w,a_v=d_u-d_w
      // convention (fit_k_sa2_per_coupling.py docstring).
      double triples[3][2] = {{d1, d2}, {d2, d0}, {d0, d1}};
      for (auto& p : triples) {
        double s = p[0] + p[1];
        double a = p[0] - p[1];
        if (s < s_min) {
          s_min = s;
          s_min_face = face;
          s_min_deg[0] = degree[lattice.faces[face].sites[0]];
          s_min_deg[1] = degree[lattice.faces[face].sites[1]];
          s_min_deg[2] = degree[lattice.faces[face].sites[2]];
          s_min_angles_deg[0] = ang0 * 180.0 / M_PI;
          s_min_angles_deg[1] = ang1 * 180.0 / M_PI;
          s_min_angles_deg[2] = ang2 * 180.0 / M_PI;
        }
        if (s > s_max) {
          s_max = s;
          s_max_face = face;
          s_max_deg[0] = degree[lattice.faces[face].sites[0]];
          s_max_deg[1] = degree[lattice.faces[face].sites[1]];
          s_max_deg[2] = degree[lattice.faces[face].sites[2]];
        }
        if (a < a_min) a_min = a;
        if (a > a_max) {
          a_max = a;
          a_max_face = face;
          a_max_deg[0] = degree[lattice.faces[face].sites[0]];
          a_max_deg[1] = degree[lattice.faces[face].sites[1]];
          a_max_deg[2] = degree[lattice.faces[face].sites[2]];
        }
      }
      for (double d : {d0, d1, d2}) {
        if (fabs(d) > max_abs_d) max_abs_d = fabs(d);
      }
    }

    printf("k=%-4d n_faces=%-8d s in [%+.5f, %+.5f]  a in [%+.5f, %+.5f]  max|d|=%.5f\n",
           k, lattice.n_faces, s_min, s_max, a_min, a_max, max_abs_d);
    printf("       s_max face=%-6d vertex_degrees=(%.0f,%.0f,%.0f)  "
           "a_max face=%-6d vertex_degrees=(%.0f,%.0f,%.0f)\n",
           s_max_face, s_max_deg[0], s_max_deg[1], s_max_deg[2],
           a_max_face, a_max_deg[0], a_max_deg[1], a_max_deg[2]);
    printf("       s_min face=%-6d vertex_degrees=(%.0f,%.0f,%.0f)  angles_deg=(%.3f,%.3f,%.3f)\n",
           s_min_face, s_min_deg[0], s_min_deg[1], s_min_deg[2],
           s_min_angles_deg[0], s_min_angles_deg[1], s_min_angles_deg[2]);

    if (s_min < global_s_min) global_s_min = s_min;
    if (s_max > global_s_max) global_s_max = s_max;
    if (a_min < global_a_min) global_a_min = a_min;
    if (a_max > global_a_max) global_a_max = a_max;
  }

  printf("\nGLOBAL (over all cached eqrp meshes): s in [%+.5f, %+.5f]  a in [%+.5f, %+.5f]\n",
         global_s_min, global_s_max, global_a_min, global_a_max);
  return 0;
}
