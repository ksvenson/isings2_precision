// real_space_2pt_test.cc
// Direct real-space measurement of <s_i s_j>, bypassing the
// spherical-harmonic/Legendre-recursion machinery entirely -- built
// 2026-08-25 per explicit user direction after the two-operator F_l(a) fit
// found a stable second dimension that was tested and rejected as epsilon
// (see IsingS2_precision/journal.md 2026-08-25 "REJECTED" entry); user
// wants Delta extracted directly from the correlator's shape vs. angular
// distance, not via F_l ratios/recursion, as an independent check for an
// analysis-pipeline issue.
//
// **Rewritten 2026-08-26, per explicit user direction ("dont do binning
// like this" -> "save the whole two point function averaged over the 120
// orbit" + "accumulated stats way where mean/variance of step n can be
// updated to step n+1"): dropped the z=cos(theta) histogram-binning scheme
// entirely** (both the original pooled-bin path and the bin-representative
// --orbit_avg path it grew into -- see git history for the dropped code).
// Binning by an arbitrary z-window pools pairs that are only approximately
// equivalent (exactly equivalent only in the continuum SO(3) limit, not on
// the discrete icosahedral mesh -- confirmed by the 6-sigma shortest-bin
// disagreement documented in journal.md's "built --orbit_avg" entry). The
// replacement here uses EXACT lattice symmetry instead of an angular
// window to decide which pairs get pooled together, and never mixes pairs
// that are not provably equivalent under QfeLatticeS2::G.
//
// Method: pick n_ref random reference sites. For each reference site i,
// partition every other site j into equivalence classes under Stab(i) (the
// subgroup of the 120-element icosahedral point group G that fixes site i
// exactly, found via nearest-site hash lookup) -- two targets j, j' land in
// the same class iff some g in Stab(i) maps j exactly onto j'. This is the
// full set of targets that MUST have identical <s_i s_j> by exact lattice
// symmetry, given a FIXED reference site i -- no continuum approximation,
// no arbitrary window. For each resulting class, one canonical
// representative pair (i, j) is kept, and its full 120-element G-orbit
// (both i and j moving together, not just Stab(i)) is precomputed once via
// the spatial hash -- this is the "whole two-point function" class: on the
// exact SO(3) sphere, every image of a given (i,j) pair under a rotation
// has the same two-point value, so per-configuration this driver averages
// s[sref]*s[starg] over the ENTIRE orbit (up to 120 pair-images, deduped)
// before feeding it into the accumulator, giving the maximal possible
// symmetry-exact variance reduction (test_orbit_avg.cc's smoke test found
// ~4-5x for a single fixed pair, applied here to every class, not just one
// bin representative per z-window as the old --orbit_avg path did).
//
// Each class's per-configuration orbit-averaged value is fed into a
// QfeMeasReal (include/statistics.h) exactly the way every other measured
// observable in this campaign is -- an online accumulator that only ever
// stores running sum/sum2/n and updates them one configuration at a time
// (Measure() called once per class per config), never buffering raw
// per-configuration samples. This satisfies "accumulated stats, mean/
// variance of step n updated to step n+1" without new machinery.
//
// Optional --theta_min/--theta_max restrict which reference-site targets
// are considered at all (default: full sphere, 0 to pi, minus self-pairs).
// This is NOT a re-introduction of binning -- it does not pool inequivalent
// pairs together, it only skips building classes outside the window
// entirely. It exists purely for compute-cost control: building a class
// costs O(|G|)=120 hash lookups per NEW equivalence class, and for a
// generic (non-high-symmetry) reference site Stab(i) is trivial, so the
// number of classes per reference site is close to n_sites -- at large
// n_refine (e.g. 128, n_sites~1.6e5) with n_ref~300 reference sites and no
// window, this is O(1e10) hash lookups in the setup phase before any MC
// sweep even starts. Narrow the window (e.g. to the fit's usual
// theta in [0.3, 2.8]) for large-mesh production runs; leave it at the
// full-sphere default only for smaller meshes or when the cost is
// acceptable.
//
// Reuses the coupling-assignment loop and mesh_mode/mesh_cache logic
// verbatim from ising_s2_crit.cc (kept in sync by hand, not shared code --
// this is a standalone diagnostic, mirroring fem_scalar_test.cc's
// precedent). Output: one line per (reference site, equivalence class) --
// ref_site, target_site, theta, mean, err, n, orbit_size -- to stdout and
// to <data_dir>/<run_id>_real_space_2pt_<seed_hex>.dat. No z_lo/z_hi/bin
// columns any more -- theta is the pair's own exact angular separation,
// not a bin center.
//
// Usage: real_space_2pt_test --q 5 --n_refine 32 --coupling_rule exact_sinh
//   --mesh_mode equal_rp --n_ref 200 --data_dir <dir> --n_therm 2000
//   --n_traj 20000 --n_skip 2 --seed 1234 [--theta_min 0.3 --theta_max 2.8]
#include <getopt.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "ising.h"
#include "S2.h"
#include "statistics.h"

// Spatial hash grid over unit-sphere site positions, for O(1) nearest-site
// lookup when mapping a rotated (G[g]*r[s]) position back to the lattice
// site it must exactly coincide with (the mesh is constructed to be
// icosahedrally symmetric to machine/near-machine precision -- see
// CLAUDE.md's "Mesh construction modes" section -- so an exact-match
// nearest-neighbor search is the correct, not approximate, operation
// here). Cell size is set from the mesh's own mean edge length so a
// 3x3x3-cell neighborhood always contains the true nearest site.
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
#include "timer.h"
#include "util.h"

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = 8;
  unsigned int seed = 1234u;
  bool cold_start = false;
  int n_therm = 2000;
  int n_traj = 20000;
  int n_skip = 2;
  int n_wolff = 5;
  int n_metropolis = 4;
  int n_ref = 200;
  double theta_min = 0.0;
  double theta_max = M_PI;
  std::string data_dir = "real_space_2pt_test";
  std::string coupling_rule = "exact_sinh";
  std::string mesh_mode = "naive";
  int equal_area_iters = 20000;
  double equal_area_step = 0.3;
  std::string mesh_cache_dir = "";

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"seed", required_argument, 0, 'S'},
      {"cold_start", no_argument, 0, 'C'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"data_dir", required_argument, 0, 'd'},
      {"n_ref", required_argument, 0, 'r'},
      {"coupling_rule", required_argument, 0, 1000},
      {"mesh_mode", required_argument, 0, 1002},
      {"equal_area_iters", required_argument, 0, 1003},
      {"equal_area_step", required_argument, 0, 1004},
      {"mesh_cache_dir", required_argument, 0, 1005},
      {"theta_min", required_argument, 0, 1006},
      {"theta_max", required_argument, 0, 1007},
      {0, 0, 0, 0}};
  const char* short_options = "q:N:S:Ch:t:s:w:e:d:r:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;
    switch (c) {
      case 'q': q = atoi(optarg); break;
      case 'N': n_refine = atoi(optarg); break;
      case 'S': seed = atol(optarg); break;
      case 'C': cold_start = true; break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      case 'd': data_dir = optarg; break;
      case 'r': n_ref = atoi(optarg); break;
      case 1000: coupling_rule = optarg; break;
      case 1002: mesh_mode = optarg; break;
      case 1003: equal_area_iters = atoi(optarg); break;
      case 1004: equal_area_step = std::stod(optarg); break;
      case 1005: mesh_cache_dir = optarg; break;
      case 1006: theta_min = std::stod(optarg); break;
      case 1007: theta_max = std::stod(optarg); break;
      default: break;
    }
  }

  if (mesh_mode != "naive" && mesh_mode != "equal_area" && mesh_mode != "equal_rp") {
    fprintf(stderr, "unknown --mesh_mode '%s'\n", mesh_mode.c_str());
    return 1;
  }
  if (coupling_rule != "duality" && coupling_rule != "exact_sinh" && coupling_rule != "owen_dual") {
    fprintf(stderr, "unknown --coupling_rule '%s'\n", coupling_rule.c_str());
    return 1;
  }

  std::string mesh_mode_suffix = (mesh_mode == "equal_area") ? "_eqarea"
                                  : (mesh_mode == "equal_rp") ? "_eqrp" : "";
  std::string run_id = string_format("q%dk%d%s", q, n_refine, mesh_mode_suffix.c_str());
  printf("run_id: %s coupling_rule: %s mesh_mode: %s n_ref: %d theta: [%.4f,%.4f]\n",
         run_id.c_str(), coupling_rule.c_str(), mesh_mode.c_str(), n_ref, theta_min, theta_max);

  QfeLatticeS2 lattice(q, n_refine);

  if (mesh_mode == "equal_area" || mesh_mode == "equal_rp") {
    std::string cache_path;
    bool cache_hit = false;
    if (!mesh_cache_dir.empty()) {
      cache_path = string_format("%s/q%dk%d%s_step%.3f.dat", mesh_cache_dir.c_str(), q,
                                  n_refine, mesh_mode_suffix.c_str(), equal_area_step);
      FILE* cache_in = fopen(cache_path.c_str(), "rb");
      if (cache_in != nullptr) {
        cache_hit = lattice.ReadPositions(cache_in);
        fclose(cache_in);
        printf("mesh_cache: %s %s\n", cache_hit ? "hit" : "stale/unreadable, recomputing", cache_path.c_str());
      } else {
        printf("mesh_cache: miss %s\n", cache_path.c_str());
      }
    }
    if (!cache_hit) {
      int iters_used = 0;
      Timer relax_timer;
      if (mesh_mode == "equal_area") {
        lattice.EqualizeFaceAreas(equal_area_iters, equal_area_step, &iters_used);
      } else {
        lattice.EqualizeCircumPerim(equal_area_iters, equal_area_step, &iters_used);
      }
      relax_timer.Stop();
      printf("%s: iters_used=%d time=%.3fs\n", mesh_mode.c_str(), iters_used, relax_timer.Duration());
      if (!mesh_cache_dir.empty()) {
        FILE* cache_out = fopen(cache_path.c_str(), "wb");
        if (cache_out != nullptr) {
          lattice.WritePositions(cache_out);
          fclose(cache_out);
        }
      }
    }
  }

  lattice.SeedRng(seed);
  lattice.UpdateWeights();
  printf("n_refine: %d q: %d total sites: %d\n", n_refine, q, lattice.n_sites);

  QfeIsing field(&lattice, 1.0);
  if (cold_start) {
    field.ColdStart();
  } else {
    field.HotStart();
  }

  // coupling assignment (verbatim copy of ising_s2_crit.cc's link loop --
  // see that file for the derivation comments)
  for (int l0 = 0; l0 < lattice.n_links; l0++) {
    int f0 = lattice.links[l0].faces[0];
    int f1 = lattice.links[l0].faces[1];
    int s0 = lattice.links[l0].sites[0];
    int s1 = lattice.links[l0].sites[1];
    int s2;
    for (int e = 0; e < 3; e++) {
      s2 = lattice.faces[f0].sites[e];
      if (s2 != s0 && s2 != s1) break;
    }
    int s3;
    for (int e = 0; e < 3; e++) {
      s3 = lattice.faces[f1].sites[e];
      if (s3 != s0 && s3 != s1) break;
    }
    int l1 = lattice.FindLink(s0, s2);
    int l2 = lattice.FindLink(s1, s2);
    int l3 = lattice.FindLink(s0, s3);
    int l4 = lattice.FindLink(s1, s3);

    Vec3 v_f0 = lattice.FaceCircumcenter(f0);
    Vec3 v_f1 = lattice.FaceCircumcenter(f1);
    Vec3 v_l0 = 0.5 * (lattice.r[s0] + lattice.r[s1]);
    Vec3 v_l1 = 0.5 * (lattice.r[s0] + lattice.r[s2]);
    Vec3 v_l2 = 0.5 * (lattice.r[s1] + lattice.r[s2]);
    Vec3 v_l3 = 0.5 * (lattice.r[s0] + lattice.r[s3]);
    Vec3 v_l4 = 0.5 * (lattice.r[s1] + lattice.r[s3]);
    Vec3 v_l0_f0 = (v_f0 - v_l0).normalized();
    Vec3 v_l0_f1 = (v_f1 - v_l0).normalized();
    Vec3 v_f0_l1 = (v_l1 - v_f0).normalized();
    Vec3 v_f0_l2 = (v_l2 - v_f0).normalized();
    Vec3 v_f1_l3 = (v_l3 - v_f1).normalized();
    Vec3 v_f1_l4 = (v_l4 - v_f1).normalized();
    double cos1 = sqrt(0.5 * (1.0 + v_l0_f0.dot(v_f0_l1)));
    double cos2 = sqrt(0.5 * (1.0 + v_l0_f0.dot(v_f0_l2)));
    double cos3 = sqrt(0.5 * (1.0 + v_l0_f1.dot(v_f1_l3)));
    double cos4 = sqrt(0.5 * (1.0 + v_l0_f1.dot(v_f1_l4)));
    double cos12 = sqrt(0.5 * (1.0 - v_f0_l1.dot(v_f0_l2)));
    double cos34 = sqrt(0.5 * (1.0 - v_f1_l3.dot(v_f1_l4)));
    double cos_prod_num = cos1 * cos2 * cos3 * cos4;
    double cos_prod_den = cos12 * cos34;

    double len_l0_sq = lattice.EdgeSquared(l0);
    double len_l0 = sqrt(len_l0_sq);
    double len_l1 = lattice.EdgeLength(l1);
    double len_l2 = lattice.EdgeLength(l2);
    double len_l3 = lattice.EdgeLength(l3);
    double len_l4 = lattice.EdgeLength(l4);
    double len_f0 = len_l0 + len_l1 + len_l2;
    double len_f1 = len_l0 + len_l3 + len_l4;

    double tanh_sq_L_num = 4.0 * len_l0_sq * cos_prod_num;
    double tanh_sq_L_den = len_f0 * len_f1 * cos_prod_den;
    double L = atanh(sqrt(tanh_sq_L_num / tanh_sq_L_den));
    double K_duality = 0.5 * asinh(1.0 / sinh(2.0 * L));

    double cos_theta_f0 = (len_l1 * len_l1 + len_l2 * len_l2 - len_l0_sq) / (2.0 * len_l1 * len_l2);
    double cos_theta_f1 = (len_l3 * len_l3 + len_l4 * len_l4 - len_l0_sq) / (2.0 * len_l3 * len_l4);
    double theta_f0 = acos(cos_theta_f0);
    double theta_f1 = acos(cos_theta_f1);
    double K_exact_sinh_f0 = 0.5 * asinh(1.0 / tan(theta_f0));
    double K_exact_sinh_f1 = 0.5 * asinh(1.0 / tan(theta_f1));
    double K_exact_sinh = 0.5 * (K_exact_sinh_f0 + K_exact_sinh_f1);

    double len_dual = (v_f0 - v_f1).norm();
    double K_owen_dual = 0.5 * asinh(len_dual / len_l0);

    double K = (coupling_rule == "exact_sinh") ? K_exact_sinh
               : (coupling_rule == "owen_dual") ? K_owen_dual : K_duality;
    lattice.links[l0].wt = K;
  }

  printf("initial action: %.12f\n", field.Action());

  // pick n_ref reference sites
  n_ref = std::min(n_ref, lattice.n_sites);
  std::mt19937 ref_rng(seed + 2);
  std::vector<int> all_sites(lattice.n_sites);
  for (int s = 0; s < lattice.n_sites; s++) all_sites[s] = s;
  std::shuffle(all_sites.begin(), all_sites.end(), ref_rng);
  std::vector<int> ref_sites(all_sites.begin(), all_sites.begin() + n_ref);

  double mean_edge_len = 0.0;
  for (int l = 0; l < lattice.n_links; l++) mean_edge_len += lattice.EdgeLength(l);
  mean_edge_len /= lattice.n_links;
  SpatialHash hash;
  hash.Build(lattice, 2.0 * mean_edge_len);

  // For each reference site i, build the exact-symmetry equivalence
  // classes described at the top of this file: partition targets j under
  // Stab(i) (cheap dedup, no continuum approximation), then for each class
  // precompute its canonical representative's FULL 120-element G-orbit
  // (both i and j moving together) for per-configuration averaging.
  struct OrbitClass {
    std::vector<std::pair<int, int>> pairs;  // (sref, starg) images, deduped
    int target_rep;
    double theta;
  };
  std::vector<std::vector<OrbitClass>> ref_classes(n_ref);

  Timer setup_timer;
  long total_orbit_size = 0, total_classes = 0;
  for (int r = 0; r < n_ref; r++) {
    int i = ref_sites[r];
    Vec3 ri = lattice.r[i].normalized();

    // Stab(i): elements of G that fix site i exactly (identity always
    // included since G[0] is the identity element).
    std::vector<int> stab;
    for (size_t g = 0; g < lattice.G.size(); g++) {
      if (hash.FindNearest(lattice.G[g] * ri) == i) stab.push_back(int(g));
    }

    std::vector<bool> visited(lattice.n_sites, false);
    visited[i] = true;
    for (int j = 0; j < lattice.n_sites; j++) {
      if (visited[j]) continue;

      double z = ri.dot(lattice.r[j].normalized());
      z = std::min(1.0, std::max(-1.0, z));
      double theta = acos(z);
      if (theta < theta_min || theta > theta_max) {
        visited[j] = true;  // out of window: skip, don't bother deduping
        continue;
      }

      // Dedup this reference site's remaining unvisited targets that are
      // Stab(i)-equivalent to j (all have identical <s_i s_j> by exact
      // lattice symmetry given fixed i).
      Vec3 rj = lattice.r[j].normalized();
      for (int g : stab) visited[hash.FindNearest(lattice.G[g] * rj)] = true;

      // Full 120-element G-orbit of the (i,j) pair, deduped, for maximal
      // per-configuration variance reduction.
      std::vector<std::pair<int, int>> pairs;
      pairs.reserve(lattice.G.size());
      for (size_t g = 0; g < lattice.G.size(); g++) {
        int sref = hash.FindNearest(lattice.G[g] * ri);
        int starg = hash.FindNearest(lattice.G[g] * rj);
        pairs.push_back({sref, starg});
      }
      std::sort(pairs.begin(), pairs.end());
      pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

      ref_classes[r].push_back({pairs, j, theta});
      total_orbit_size += pairs.size();
      total_classes++;
    }
  }
  setup_timer.Stop();
  printf("exact-orbit setup: %ld classes, mean orbit size %.2f (|G|=%zu), time=%.3fs\n",
         total_classes, double(total_orbit_size) / std::max(total_classes, 1L), lattice.G.size(),
         setup_timer.Duration());

  std::vector<std::vector<QfeMeasReal>> ref_meas(n_ref);
  for (int r = 0; r < n_ref; r++) ref_meas[r].resize(ref_classes[r].size());

  Timer timer;
  for (int n = 0; n < (n_therm + n_traj); n++) {
    for (int j = 0; j < n_metropolis; j++) field.Metropolis();
    for (int j = 0; j < n_wolff; j++) field.WolffUpdate();

    if (n % n_skip || n < n_therm) continue;

    for (int r = 0; r < n_ref; r++) {
      const std::vector<OrbitClass>& classes = ref_classes[r];
      for (size_t c = 0; c < classes.size(); c++) {
        double acc = 0.0;
        for (const auto& p : classes[c].pairs) acc += field.spin[p.first] * field.spin[p.second];
        acc /= double(classes[c].pairs.size());
        // Online mean/variance update (QfeMeasReal::Measure), one call per
        // class per measured configuration -- no raw per-config buffering.
        ref_meas[r][c].Measure(acc);
      }
    }
  }
  timer.Stop();
  printf("measurement loop time: %.3fs\n", timer.Duration());

  std::string out_path = string_format("%s/%s_real_space_2pt_%08X.dat", data_dir.c_str(), run_id.c_str(), seed);
  FILE* out = fopen(out_path.c_str(), "w");
  if (out == nullptr) {
    fprintf(stderr, "could not open %s for writing\n", out_path.c_str());
    out = stdout;
  }
  fprintf(out, "# ref_site target_site theta mean err n orbit_size\n");
  for (int r = 0; r < n_ref; r++) {
    int i = ref_sites[r];
    const std::vector<OrbitClass>& classes = ref_classes[r];
    for (size_t c = 0; c < classes.size(); c++) {
      fprintf(out, "%6d %6d %.6f %+.10e %.10e %d %zu\n", i, classes[c].target_rep, classes[c].theta,
              ref_meas[r][c].Mean(), ref_meas[r][c].Error(), ref_meas[r][c].n, classes[c].pairs.size());
    }
  }
  if (out != stdout) fclose(out);

  return 0;
}
