// analyze_north_pole_configs.cc
// Offline, from bin/save_configs's saved (positions, bit-packed configs)
// shards: the two-point function <s_0 s_j> centered on the north-pole
// site (the site with max z, exactly [0,0,1] for this campaign's meshes),
// averaged over the exact 120-element icosahedral G-orbit of each
// (north_pole, j) pair -- built 2026-08-26 per user direction ("can we
// build both the two point function (centered on the north pole and
// averaged over the 120) and the ylm matrix from these configs").
//
// This deliberately reuses the SAME exact-equivalence-class construction
// as real_space_2pt_test.cc's rewrite and full_corr_test.cc's
// SymmetrizeOverOrbits (Stab(i)-based target dedup + full-G-orbit
// averaging for the MC value) -- but applied to exactly ONE reference row
// (the north pole), computed OFFLINE from already-saved configs rather
// than online during MC. Since it's one row, not the full n_sites x
// n_sites matrix, storage/compute is O(n_sites * |G|) for setup and
// O(n_meas * n_classes) for the measurement pass -- no O(n_sites^2) wall
// at any n_refine, unlike full_corr_test.cc/compute_two_point_stats.py.
//
// Per shard, this driver computes the shard's own MC-time-averaged value
// per class (one pass over that shard's saved configs); shards are then
// combined via Welford's algorithm (mean + variance across independent
// shards), matching compute_two_point_stats.py's pattern and this
// campaign's "accumulated stats, mean/variance of step n updated to step
// n+1" convention.
//
// Positions are read from ONE save_configs positions.dat (text format,
// identical across shards of the same q/n_refine/mesh_mode -- mesh
// relaxation is deterministic, not seed-dependent). The north pole is
// found by argmax(z), not assumed to be site 0, since mesh relaxation
// could in principle renumber/perturb this (it doesn't, in practice, but
// this is cheap to check robustly rather than assume).
//
// Usage: analyze_north_pole_configs --q 5 --n_refine 32 \
//   --positions <run_id>_positions_<seed>.dat \
//   --configs <run_id>_configs_<seed1>.bin --configs <run_id>_configs_<seed2>.bin ... \
//   -o <out_path>
#include <getopt.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "ising.h"
#include "S2.h"
#include "statistics.h"

// Verbatim copy of real_space_2pt_test.cc's SpatialHash (kept in sync by
// hand, not shared code -- same convention as every other standalone
// diagnostic in this campaign).
class SpatialHash {
 public:
  void Build(const QfeLatticeS2& lattice, double cell_size) {
    lattice_ = &lattice;
    cell_size_ = cell_size;
    cells_.clear();
    for (int s = 0; s < lattice.n_sites; s++) {
      cells_[CellKeyForPoint(lattice.r[s].normalized())].push_back(s);
    }
  }

  int FindNearest(const Vec3& p) const {
    Vec3 pn = p.normalized();
    int cx = CellCoord(pn.x());
    int cy = CellCoord(pn.y());
    int cz = CellCoord(pn.z());
    int best = -1;
    double best_dot = -2.0;
    for (int dx = -1; dx <= 1; dx++) {
      for (int dy = -1; dy <= 1; dy++) {
        for (int dz = -1; dz <= 1; dz++) {
          auto it = cells_.find(CellKey(cx + dx, cy + dy, cz + dz));
          if (it == cells_.end()) continue;
          for (int s : it->second) {
            double d = lattice_->r[s].normalized().dot(pn);
            if (d > best_dot) {
              best_dot = d;
              best = s;
            }
          }
        }
      }
    }
    return best;
  }

 private:
  int CellCoord(double x) const { return int(std::floor(x / cell_size_)); }
  int64_t CellKey(int ix, int iy, int iz) const {
    auto enc = [](int v) -> int64_t { return int64_t(v) & 0x1FFFFF; };
    return (enc(ix) << 42) | (enc(iy) << 21) | enc(iz);
  }
  int64_t CellKeyForPoint(const Vec3& p) const {
    return CellKey(CellCoord(p.x()), CellCoord(p.y()), CellCoord(p.z()));
  }

  const QfeLatticeS2* lattice_ = nullptr;
  double cell_size_ = 1.0;
  std::unordered_map<int64_t, std::vector<int>> cells_;
};
#include "util.h"

bool LoadTextPositions(const std::string& path, QfeLatticeS2* lattice) {
  FILE* f = fopen(path.c_str(), "r");
  if (f == nullptr) return false;
  char line[256];
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '#') continue;
    int s;
    double x, y, z, wt;
    if (sscanf(line, "%d %lf %lf %lf %lf", &s, &x, &y, &z, &wt) == 5) {
      if (s < 0 || s >= lattice->n_sites) continue;
      lattice->r[s] = Vec3(x, y, z);
    }
  }
  fclose(f);
  return true;
}

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = 8;
  int l_max_quad = 8;
  std::string positions_path;
  std::vector<std::string> configs_paths;
  std::string out_path;

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"positions", required_argument, 0, 'p'},
      {"configs", required_argument, 0, 'c'},
      {"out", required_argument, 0, 'o'},
      {"l_max_quad", required_argument, 0, 1000},
      {0, 0, 0, 0}};

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, "q:N:p:c:o:", long_options, &o);
    if (c == -1) break;
    switch (c) {
      case 'q': q = atoi(optarg); break;
      case 'N': n_refine = atoi(optarg); break;
      case 'p': positions_path = optarg; break;
      case 'c': configs_paths.push_back(optarg); break;
      case 'o': out_path = optarg; break;
      case 1000: l_max_quad = atoi(optarg); break;
      default: break;
    }
  }
  if (positions_path.empty() || configs_paths.empty() || out_path.empty()) {
    fprintf(stderr, "usage: --q Q --n_refine N --positions FILE --configs FILE [--configs FILE ...] -o OUT\n");
    return 1;
  }

  QfeLatticeS2 lattice(q, n_refine);
  if (!LoadTextPositions(positions_path, &lattice)) {
    fprintf(stderr, "could not read positions file %s\n", positions_path.c_str());
    return 1;
  }
  int n_sites = lattice.n_sites;
  printf("n_refine: %d q: %d total sites: %d\n", n_refine, q, n_sites);

  // Quadrature weights -- added 2026-08-26 per user direction: the
  // per-class Legendre projection was using raw target-site COUNTS
  // (n_targets) as an equal-weight proxy for solid angle, which is only
  // approximately right and most wrong exactly where mesh non-uniformity
  // is largest (coarse n_refine) -- suspected cause of the F_l method's
  // n_refine=4 point going the wrong direction (below 1/8 instead of
  // pushed above it by finite-lattice-spacing effects, unlike every
  // other estimator in this campaign). First fix (UpdateWeights, the
  // cotangent-Laplacian lumped-Voronoi-dual-area scheme every other
  // driver in this codebase uses) only partially closed the gap
  // (0.1196 -> 0.1225 at n_refine=4, still short of 0.125). Switched to
  // OptimizeIntegrator(l_max), which solves for weights that integrate
  // spherical harmonics EXACTLY up to l_max (not just approximately, per
  // CLAUDE.md's note that no other driver in this campaign actually
  // calls it) -- a strictly more accurate quadrature for this
  // Legendre-projection use case. See journal.md.
  lattice.OptimizeIntegrator(l_max_quad);

  // north pole = site with max z (should be exactly [0,0,1]; found
  // robustly rather than assumed to be site 0)
  int np = 0;
  for (int s = 1; s < n_sites; s++) {
    if (lattice.r[s].z() > lattice.r[np].z()) np = s;
  }
  printf("north pole site: %d at (%.6f, %.6f, %.6f)\n", np, lattice.r[np].x(), lattice.r[np].y(), lattice.r[np].z());

  double mean_edge_len = 0.0;
  for (int l = 0; l < lattice.n_links; l++) mean_edge_len += lattice.EdgeLength(l);
  mean_edge_len /= lattice.n_links;
  SpatialHash hash;
  hash.Build(lattice, 2.0 * mean_edge_len);

  // Stab(north pole): elements of G fixing site np exactly
  Vec3 rnp = lattice.r[np].normalized();
  std::vector<int> stab;
  for (size_t g = 0; g < lattice.G.size(); g++) {
    if (hash.FindNearest(lattice.G[g] * rnp) == np) stab.push_back(int(g));
  }
  printf("|Stab(north pole)| = %zu (|G|=%zu)\n", stab.size(), lattice.G.size());

  // Build exact equivalence classes of targets (Stab(np)-dedup) and each
  // class's full 120-element G-orbit of the pair (np, target), same
  // construction as real_space_2pt_test.cc.
  struct OrbitClass {
    std::vector<std::pair<int, int>> pairs;
    int target_rep;
    double theta;
    int n_targets;  // number of distinct target SITES sharing this theta
                    // under Stab(np) -- i.e. this class's solid-angle
                    // multiplicity, needed to Legendre-project <s(theta)>
                    // correctly (added 2026-08-26, per user request to
                    // project the north-pole two-point function onto
                    // Y_lm/F_l -- see scripts/fl_delta_from_north_pole.py).
                    // NOT the same as pairs.size() (the full 120-orbit used
                    // for MC-averaging, which also moves the reference).
    double sum_wt;  // sum of true quadrature weight (UpdateWeights) over
                    // this class's target sites -- the correct solid-angle
                    // weight, replacing the n_targets equal-weight
                    // approximation (added 2026-08-26, same reason).
  };
  std::vector<OrbitClass> classes;
  std::vector<bool> visited(n_sites, false);
  visited[np] = true;
  for (int j = 0; j < n_sites; j++) {
    if (visited[j]) continue;
    Vec3 rj = lattice.r[j].normalized();
    double z = std::min(1.0, std::max(-1.0, rnp.dot(rj)));
    double theta = acos(z);

    std::vector<int> stab_targets;
    for (int g : stab) stab_targets.push_back(hash.FindNearest(lattice.G[g] * rj));
    std::sort(stab_targets.begin(), stab_targets.end());
    stab_targets.erase(std::unique(stab_targets.begin(), stab_targets.end()), stab_targets.end());
    double sum_wt = 0.0;
    for (int t : stab_targets) {
      visited[t] = true;
      sum_wt += lattice.sites[t].wt;
    }

    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(lattice.G.size());
    for (size_t g = 0; g < lattice.G.size(); g++) {
      int sref = hash.FindNearest(lattice.G[g] * rnp);
      int starg = hash.FindNearest(lattice.G[g] * rj);
      pairs.push_back({sref, starg});
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

    classes.push_back({pairs, j, theta, int(stab_targets.size()), sum_wt});
  }
  printf("classes: %zu (vs n_sites-1=%d)\n", classes.size(), n_sites - 1);

  int n_classes = int(classes.size());
  std::vector<double> mean(n_classes, 0.0), M2(n_classes, 0.0);
  int n_shards = 0;
  long total_configs = 0;

  QfeIsing field(&lattice, 1.0);
  int buf_size = (n_sites + 7) / 8;
  std::vector<unsigned char> spin_buf(buf_size);
  std::vector<double> spin(n_sites);

  for (const std::string& cp : configs_paths) {
    FILE* f = fopen(cp.c_str(), "rb");
    if (f == nullptr) {
      fprintf(stderr, "could not open %s\n", cp.c_str());
      return 1;
    }
    std::vector<double> class_sum(n_classes, 0.0);
    long n_meas = 0;
    while (fread(spin_buf.data(), 1, buf_size, f) == size_t(buf_size)) {
      char mask = 1;
      for (int i = 0; i < n_sites; i++) {
        spin[i] = (spin_buf[i / 8] & mask) ? -1.0 : 1.0;
        mask <<= 1;
        if (!mask) mask = 1;
      }
      for (int k = 0; k < n_classes; k++) {
        double acc = 0.0;
        for (const auto& p : classes[k].pairs) acc += spin[p.first] * spin[p.second];
        class_sum[k] += acc / double(classes[k].pairs.size());
      }
      n_meas++;
    }
    fclose(f);
    total_configs += n_meas;

    n_shards++;
    for (int k = 0; k < n_classes; k++) {
      double val = class_sum[k] / double(n_meas);
      double delta = val - mean[k];
      mean[k] += delta / n_shards;
      double delta2 = val - mean[k];
      M2[k] += delta * delta2;
    }
    printf("shard %s: n_meas=%ld\n", cp.c_str(), n_meas);
  }

  if (n_shards < 2) {
    fprintf(stderr, "need >=2 shards for a variance estimate, got %d\n", n_shards);
    return 1;
  }

  FILE* out = fopen(out_path.c_str(), "w");
  if (out == nullptr) {
    fprintf(stderr, "could not open %s for writing\n", out_path.c_str());
    return 1;
  }
  fprintf(out, "# north_pole_site=%d n_sites=%d n_shards=%d total_configs=%ld\n", np, n_sites, n_shards,
          total_configs);
  fprintf(out, "# target_site theta mean sem orbit_size n_targets sum_wt\n");
  for (int k = 0; k < n_classes; k++) {
    double variance = M2[k] / (n_shards - 1);
    double sem = sqrt(variance / n_shards);
    fprintf(out, "%6d %.6f %+.10e %.10e %zu %d %.10e\n", classes[k].target_rep, classes[k].theta, mean[k], sem,
            classes[k].pairs.size(), classes[k].n_targets, classes[k].sum_wt);
  }
  fclose(out);
  printf("wrote %s\n", out_path.c_str());

  return 0;
}
