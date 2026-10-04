// full_corr_test.cc
// Accumulates the MC-averaged full all-pairs two-point function <s_i s_j>
// (every site pair, not a sampled subset) -- built 2026-08-26 per explicit
// user direction: "we will build the [spherical-harmonic decomposition]
// matrices from the MC averaged two point functions. Building them in the
// simulator was a mistake since the two point functions were not smoothed
// by high stats yet." I.e. the M_l[m,m'] matrices this campaign's
// spherical-symmetry gate needs (previously accumulated directly inside
// ising_s2_crit.cc's --l_max path, contracting each raw per-configuration
// spin snapshot against Y_lm before ever forming a converged real-space
// correlator) should instead be built OFFLINE from an already-converged,
// full-statistics <s_i s_j> matrix. This driver produces exactly that
// matrix; scripts/build_ylm_matrix_from_corr.py does the offline Y_lm
// contraction into M_l.
//
// Per user follow-up ("the 120 averaging can always be done for free
// statistics"): once every (i,j) pair has actually been measured (unlike
// real_space_2pt_test.cc's n_ref-sampled reference sites, here EVERY site
// is its own reference), averaging each entry over its exact 120-element
// icosahedral G-orbit combines multiple genuinely-independent measured
// entries (not redundant re-use of a single measurement) -- a real, free
// variance reduction, applied once as a post-processing pass after MC
// accumulation (SymmetrizeOverOrbits() below), not per-configuration.
//
// Method: for each measured configuration, form the spin vector s in
// {-1,+1}^n_sites and accumulate the outer product s*s^T into the upper
// triangle of a dense Eigen::MatrixXd via a symmetric rank-1 update
// (Eigen's selfadjointView<Upper>().rankUpdate, the same operation as
// BLAS's dsyr). This is O(n_sites^2) per measured configuration -- there
// is no way around that cost if every pair is genuinely measured, unlike
// real_space_2pt_test.cc's O(n_ref*n_sites) sampling. Dense storage means
// this driver is memory-bound (~8*n_sites^2 bytes) and only practical up
// to modest n_refine (roughly n_refine<=32-48 on typical SCC memory
// limits; DO NOT run this at n_refine=128 without checking available
// memory first -- 163842^2*8 bytes ~ 215 GB).
//
// Output: <data_dir>/<run_id>_full_corr_<seed_hex>.dat -- a small text
// header (n_sites, n_meas) followed by the site positions (reuses
// QfeLatticeS2::WritePositions' binary format), the quadrature weights
// (sites[s].wt, binary doubles, needed downstream for the Y_lm
// contraction), and the symmetrized upper-triangle (i<=j) of the
// correlator matrix (binary doubles, row-major over i).
//
// Usage: full_corr_test --q 5 --n_refine 16 --coupling_rule exact_sinh
//   --mesh_mode naive --data_dir <dir> --n_therm 2000 --n_traj 20000
//   --n_skip 2 --seed 1234
#include <getopt.h>

#include <Eigen/Dense>
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

// Spatial hash grid over unit-sphere site positions, for O(1) nearest-site
// lookup when mapping a rotated (G[g]*r[s]) position back to the lattice
// site it must exactly coincide with. Verbatim copy of
// real_space_2pt_test.cc's SpatialHash (kept in sync by hand, not shared
// code -- see that file's header comment for the same convention).
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

// Post-processing pass: replace each upper-triangle entry corr(i,j) with
// the mean of corr over the exact 120-element G-orbit of the pair (i,j).
// Every orbit is visited exactly once (marked in `visited`, indexed by
// lo*n_sites+hi so lower/upper images of the same pair collide); since an
// orbit's members are always ALL visited or ALL unvisited together (orbit
// membership is symmetric), reading every member's original value before
// writing any of them back is safe -- no member is ever read after being
// overwritten.
void SymmetrizeOverOrbits(const QfeLatticeS2& lattice, const SpatialHash& hash,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>* corr) {
  int n = lattice.n_sites;
  std::vector<uint8_t> visited(size_t(n) * n, 0);
  auto lo_hi = [](int a, int b) { return a <= b ? std::make_pair(a, b) : std::make_pair(b, a); };

  for (int i = 0; i < n; i++) {
    Vec3 ri = lattice.r[i].normalized();
    for (int j = i; j < n; j++) {
      size_t vidx = size_t(i) * n + j;
      if (visited[vidx]) continue;

      Vec3 rj = lattice.r[j].normalized();
      std::vector<std::pair<int, int>> pairs;
      pairs.reserve(lattice.G.size());
      for (size_t g = 0; g < lattice.G.size(); g++) {
        int sref = hash.FindNearest(lattice.G[g] * ri);
        int starg = hash.FindNearest(lattice.G[g] * rj);
        pairs.push_back(lo_hi(sref, starg));
      }
      std::sort(pairs.begin(), pairs.end());
      pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

      double acc = 0.0;
      for (const auto& p : pairs) acc += (*corr)(p.first, p.second);
      acc /= double(pairs.size());

      for (const auto& p : pairs) {
        (*corr)(p.first, p.second) = acc;
        visited[size_t(p.first) * n + p.second] = 1;
      }
    }
  }
}

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
  std::string data_dir = "full_corr_test";
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
      {"coupling_rule", required_argument, 0, 1000},
      {"mesh_mode", required_argument, 0, 1002},
      {"equal_area_iters", required_argument, 0, 1003},
      {"equal_area_step", required_argument, 0, 1004},
      {"mesh_cache_dir", required_argument, 0, 1005},
      {0, 0, 0, 0}};
  const char* short_options = "q:N:S:Ch:t:s:w:e:d:";

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
      case 1000: coupling_rule = optarg; break;
      case 1002: mesh_mode = optarg; break;
      case 1003: equal_area_iters = atoi(optarg); break;
      case 1004: equal_area_step = std::stod(optarg); break;
      case 1005: mesh_cache_dir = optarg; break;
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
  printf("run_id: %s coupling_rule: %s mesh_mode: %s\n", run_id.c_str(), coupling_rule.c_str(),
         mesh_mode.c_str());

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

  int n = lattice.n_sites;
  size_t bytes = size_t(n) * n * sizeof(double);
  printf("allocating %d x %d correlator matrix (%.2f GB)\n", n, n, bytes / 1e9);
  // Row-major so a row's upper-triangle tail (corr.row(i).data()+i,
  // n-i doubles) is contiguous for the fwrite below -- Eigen's default
  // MatrixXd is column-major, which would make that pointer arithmetic
  // read across columns, not along the row.
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> corr(n, n);
  corr.setZero();
  Eigen::VectorXd s(n);
  long n_meas = 0;

  Timer timer;
  for (int step = 0; step < (n_therm + n_traj); step++) {
    for (int j = 0; j < n_metropolis; j++) field.Metropolis();
    for (int j = 0; j < n_wolff; j++) field.WolffUpdate();

    if (step % n_skip || step < n_therm) continue;

    for (int i = 0; i < n; i++) s(i) = field.spin[i];
    corr.selfadjointView<Eigen::Upper>().rankUpdate(s, 1.0);
    n_meas++;
  }
  corr /= double(n_meas);
  timer.Stop();
  printf("measurement loop: n_meas=%ld time=%.3fs\n", n_meas, timer.Duration());

  double mean_edge_len = 0.0;
  for (int l = 0; l < lattice.n_links; l++) mean_edge_len += lattice.EdgeLength(l);
  mean_edge_len /= lattice.n_links;
  SpatialHash hash;
  hash.Build(lattice, 2.0 * mean_edge_len);

  Timer sym_timer;
  SymmetrizeOverOrbits(lattice, hash, &corr);
  sym_timer.Stop();
  printf("orbit symmetrization: time=%.3fs\n", sym_timer.Duration());

  // sanity check: diagonal must be exactly 1.0 (s_i^2 = 1 every config)
  double max_diag_err = 0.0;
  for (int i = 0; i < n; i++) max_diag_err = std::max(max_diag_err, std::abs(corr(i, i) - 1.0));
  printf("sanity check: max |corr(i,i) - 1| = %.3e\n", max_diag_err);

  std::string out_path = string_format("%s/%s_full_corr_%08X.dat", data_dir.c_str(), run_id.c_str(), seed);
  FILE* out = fopen(out_path.c_str(), "wb");
  if (out == nullptr) {
    fprintf(stderr, "could not open %s for writing\n", out_path.c_str());
    return 1;
  }
  int32_t n32 = n;
  int64_t n_meas64 = n_meas;
  fwrite(&n32, sizeof(n32), 1, out);
  fwrite(&n_meas64, sizeof(n_meas64), 1, out);
  lattice.WritePositions(out);
  std::vector<double> wt(n);
  for (int i = 0; i < n; i++) wt[i] = lattice.sites[i].wt;
  fwrite(wt.data(), sizeof(double), n, out);
  for (int i = 0; i < n; i++) {
    fwrite(corr.row(i).data() + i, sizeof(double), n - i, out);
  }
  fclose(out);
  printf("wrote %s\n", out_path.c_str());

  return 0;
}
