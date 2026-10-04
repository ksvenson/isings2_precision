// lean_harmonic_stats.cc
// Data-volume-lean alternative to save_configs.cc + offline M_l
// reconstruction: instead of dumping raw per-configuration spins
// (O(n_meas*n_sites) bits) or the full (l_max+1)^2 x (l_max+1)^2 M_l[m,m']
// matrix, accumulate ONLINE during the MC run, per configuration:
//   - |S_lm|^2 for every (l,m), l=0..l_max, m=-l..l -- (l_max+1)^2
//     diagonal harmonic-operator expectation values (Trace(M_l) building
//     blocks; <S_lm> itself is exactly zero by the spin-flip symmetry, so
//     the second moment is the physical quantity, not a raw expectation)
//   - powers of the magnetization M = (1/n_sites) sum_i s_i: M, M^2, M^3,
//     M^4 (the existing bulk observable set _bulk_*.dat already reports,
//     duplicated here so this driver's output is self-contained)
// Never materializes a full n_sites x n_sites correlator or an M_l[m,m']
// matrix, and never writes a raw configuration to disk -- output is a
// handful of numbers per jackknife block. Built 2026-08-26 per user
// direction ("do a small simulator run with these operators accumulating"
// / "expectation values of (l+1)^2 operators plus powers of the
// magnetization to keep the data volume lean").
//
// Reuses save_configs.cc's mesh-construction/coupling-assignment code
// verbatim (kept in sync by hand, same convention as every other
// standalone driver in this campaign).
//
// Quadrature weight: uses QfeLatticeS2::UpdateWeights() (the same scheme
// save_configs.cc/build_ylm_matrix_from_configs.py already use), NOT
// OptimizeIntegrator -- OptimizeIntegrator is a more accurate quadrature
// for this Legendre-projection use case (see journal.md's 2026-08-26 "F_l
// normalization bug" entry) but has an unresolved bug where its
// least-squares solve returns a negative weight (tripping its own
// assert(real(wt(id)) > 0.0)) at n_refine=32 and above -- not fixed here,
// flagged as a separate open item. Delta_s(a) from this driver should be
// expected to show the same small residual bias UpdateWeights showed
// elsewhere in this campaign (partially, not fully, corrected).
//
// Output: <data_dir>/<run_id>_lean_<seed_hex>.dat, text, one header line
// plus one line per jackknife block:
//   block_index n_meas M_sum M2_sum M3_sum M4_sum <abs2_sum for each
//   (l,m) in row-major l=0..l_max,m=-l..l order>
#include <getopt.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "ising.h"
#include "S2.h"
#include "statistics.h"
#include "timer.h"
#include "util.h"

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = 8;
  int l_max = 8;
  unsigned int seed = 1234u;
  bool cold_start = false;
  int n_therm = 2000;
  int n_traj = 20000;
  int n_skip = 2;
  int n_wolff = 5;
  int n_metropolis = 4;
  int jack_block_size = 100;
  std::string data_dir = "lean_harmonic_stats";
  std::string coupling_rule = "exact_sinh";
  // --mesh_mode: ported from ising_s2_crit.cc 2026-08-27, per user
  // correction ("we should be comparing the naive to the equal areas
  // ALWAYS") after every result in this driver's high-stats overnight
  // ladder turned out to have used the naive mesh only, with no way to
  // compare against equal_area/equal_rp -- see ising_s2_crit.cc's
  // "mesh_mode" comment block and journal.md/CLAUDE.md 2026-08-25/27 for
  // the full motivation (naive/equal-area meshes are NOT exactly critical
  // per arXiv:2407.00459's coupling-rule derivation; equal_rp is the
  // paper's fix). Same three modes, same flags, same semantics as
  // ising_s2_crit.cc -- kept in sync by hand like the rest of this file's
  // mesh/coupling code.
  std::string mesh_mode = "naive";
  int equal_area_iters = 20000;
  double equal_area_step = 0.3;
  std::string mesh_cache_dir = "";

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"l_max", required_argument, 0, 'L'},
      {"seed", required_argument, 0, 'S'},
      {"cold_start", no_argument, 0, 'C'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"jack_block_size", required_argument, 0, 'j'},
      {"data_dir", required_argument, 0, 'd'},
      {"coupling_rule", required_argument, 0, 1000},
      {"mesh_mode", required_argument, 0, 1002},
      {"equal_area_iters", required_argument, 0, 1003},
      {"equal_area_step", required_argument, 0, 1004},
      {"mesh_cache_dir", required_argument, 0, 1005},
      {0, 0, 0, 0}};
  const char* short_options = "q:N:L:S:Ch:t:s:w:e:j:d:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;
    switch (c) {
      case 'q': q = atoi(optarg); break;
      case 'N': n_refine = atoi(optarg); break;
      case 'L': l_max = atoi(optarg); break;
      case 'S': seed = atol(optarg); break;
      case 'C': cold_start = true; break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      case 'j': jack_block_size = atoi(optarg); break;
      case 'd': data_dir = optarg; break;
      case 1000: coupling_rule = optarg; break;
      case 1002: mesh_mode = optarg; break;
      case 1003: equal_area_iters = atoi(optarg); break;
      case 1004: equal_area_step = std::stod(optarg); break;
      case 1005: mesh_cache_dir = optarg; break;
      default: break;
    }
  }

  if (coupling_rule != "duality" && coupling_rule != "exact_sinh" && coupling_rule != "owen_dual") {
    fprintf(stderr, "unknown --coupling_rule '%s'\n", coupling_rule.c_str());
    return 1;
  }

  if (mesh_mode != "naive" && mesh_mode != "equal_area" && mesh_mode != "equal_rp") {
    fprintf(stderr,
            "unknown --mesh_mode '%s' (expected 'naive', 'equal_area', or 'equal_rp')\n",
            mesh_mode.c_str());
    return 1;
  }

  std::string mesh_mode_suffix = (mesh_mode == "equal_area")
                                      ? "_eqarea"
                                      : (mesh_mode == "equal_rp") ? "_eqrp" : "";
  std::string run_id = string_format("q%dk%d%s", q, n_refine, mesh_mode_suffix.c_str());
  printf("run_id: %s coupling_rule: %s mesh_mode: %s l_max: %d\n", run_id.c_str(),
         coupling_rule.c_str(), mesh_mode.c_str(), l_max);

  QfeLatticeS2 lattice(q, n_refine);

  if (mesh_mode == "equal_area" || mesh_mode == "equal_rp") {
    std::string cache_path;
    bool cache_hit = false;
    if (!mesh_cache_dir.empty()) {
      cache_path = string_format("%s/q%dk%d%s_step%.3f.dat",
                                  mesh_cache_dir.c_str(), q, n_refine,
                                  mesh_mode_suffix.c_str(), equal_area_step);
      FILE* cache_in = fopen(cache_path.c_str(), "rb");
      if (cache_in != nullptr) {
        cache_hit = lattice.ReadPositions(cache_in);
        fclose(cache_in);
        printf("mesh_cache: %s %s\n", cache_hit ? "hit" : "stale/unreadable, recomputing",
               cache_path.c_str());
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
      printf("%s: iters_used=%d (cap=%d) time=%.3fs\n", mesh_mode.c_str(), iters_used,
             equal_area_iters, relax_timer.Duration());
      if (!mesh_cache_dir.empty()) {
        FILE* cache_out = fopen(cache_path.c_str(), "wb");
        if (cache_out != nullptr) {
          lattice.WritePositions(cache_out);
          fclose(cache_out);
          printf("mesh_cache: wrote %s\n", cache_path.c_str());
        } else {
          fprintf(stderr,
                  "mesh_cache: could not open %s for writing (does --mesh_cache_dir exist?), "
                  "continuing without caching this result\n",
                  cache_path.c_str());
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

  // coupling assignment (verbatim copy of ising_s2_crit.cc/save_configs.cc's
  // link loop -- see those files for the derivation comments)
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

  // Precompute (l,m) mode list and Y_lm(site) once -- fixed for the whole
  // run since mesh positions don't change after construction.
  std::vector<int> l_list, m_list;
  for (int l = 0; l <= l_max; l++) {
    for (int m = -l; m <= l; m++) {
      l_list.push_back(l);
      m_list.push_back(m);
    }
  }
  int n_lm = int(l_list.size());
  printf("n_lm = (l_max+1)^2 = %d\n", n_lm);

  // Off-diagonal pairs (k,kp) with l_list[k]==l_list[kp], m_list[k] <
  // m_list[kp] -- added 2026-08-27 per user direction ("make that
  // change") so this driver can compute the REAL kappa2_r/R_l (needs the
  // full M_l[m,m'] matrix, not just the diagonal |S_lm|^2 this driver
  // saved before). Diagonal (m==mp) entries are NOT duplicated here --
  // they're exactly the existing abs2_sum columns, real-valued, already
  // written below. Same (l,m,mp) convention as
  // scripts/symmetry_test.py's full_pairs_for_l / ising_s2_crit.cc's
  // full_pairs (M_l[m,mp] = <S_lm * conj(S_lmp)>, m<=mp, Hermitian so
  // this determines the rest). Only 444 pairs at l_max=8 (vs 81
  // diagonal) -- see CLAUDE.md's "how many operators" table -- cheap to
  // add since it's O(n_pairs) per measurement, not O(n_sites) (S_lm
  // itself is already computed for the diagonal).
  struct OffDiagPair { int k, kp, l, m, mp; };
  std::vector<OffDiagPair> offdiag_pairs;
  for (int k = 0; k < n_lm; k++) {
    for (int kp = k + 1; kp < n_lm; kp++) {
      if (l_list[k] != l_list[kp]) continue;
      if (m_list[kp] <= m_list[k]) continue;
      offdiag_pairs.push_back({k, kp, l_list[k], m_list[k], m_list[kp]});
    }
  }
  int n_offdiag = int(offdiag_pairs.size());
  printf("n_offdiag = %d\n", n_offdiag);

  int n_sites = lattice.n_sites;
  std::vector<std::vector<std::complex<double>>> wY(n_lm, std::vector<std::complex<double>>(n_sites));
  for (int k = 0; k < n_lm; k++) {
    for (int s = 0; s < n_sites; s++) {
      double theta = acos(lattice.r[s].z());
      double phi = atan2(lattice.r[s].y(), lattice.r[s].x());
      // conj(Y_lm) folded into the precomputed weight*Ylm table, matching
      // S_lm = sum_i w_i * conj(Y_lm(x_i)) * s_i (analyze_saved_configs.py's
      // convention)
      std::complex<double> ylm = boost::math::spherical_harmonic(l_list[k], m_list[k], theta, phi);
      wY[k][s] = lattice.sites[s].wt * std::conj(ylm);
    }
  }

  std::string out_path = string_format("%s/%s_lean_%08X.dat", data_dir.c_str(), run_id.c_str(), seed);
  FILE* out = fopen(out_path.c_str(), "w");
  if (out == nullptr) {
    fprintf(stderr, "could not open %s for writing\n", out_path.c_str());
    return 1;
  }
  fprintf(out, "# q=%d n_refine=%d n_sites=%d l_max=%d n_lm=%d n_offdiag=%d coupling_rule=%s jack_block_size=%d\n",
          q, n_refine, n_sites, l_max, n_lm, n_offdiag, coupling_rule.c_str(), jack_block_size);
  fprintf(out, "# lm_order:");
  for (int k = 0; k < n_lm; k++) fprintf(out, " (%d,%d)", l_list[k], m_list[k]);
  fprintf(out, "\n");
  fprintf(out, "# offdiag_order:");
  for (int i = 0; i < n_offdiag; i++) {
    fprintf(out, " (%d,%d,%d)", offdiag_pairs[i].l, offdiag_pairs[i].m, offdiag_pairs[i].mp);
  }
  fprintf(out, "\n");
  fprintf(out, "# block_index n_meas M_sum M2_sum M3_sum M4_sum abs2_sum[0..n_lm-1] "
               "offdiag_re[0..n_offdiag-1] offdiag_im[0..n_offdiag-1]\n");

  int block_n = 0;
  double block_M = 0, block_M2 = 0, block_M3 = 0, block_M4 = 0;
  std::vector<double> block_abs2(n_lm, 0.0);
  std::vector<double> block_offdiag_re(n_offdiag, 0.0);
  std::vector<double> block_offdiag_im(n_offdiag, 0.0);
  int block_index = 0;

  auto write_block = [&]() {
    fprintf(out, "%d %d %.16e %.16e %.16e %.16e", block_index, block_n, block_M, block_M2, block_M3, block_M4);
    for (int k = 0; k < n_lm; k++) fprintf(out, " %.16e", block_abs2[k]);
    for (int i = 0; i < n_offdiag; i++) fprintf(out, " %.16e", block_offdiag_re[i]);
    for (int i = 0; i < n_offdiag; i++) fprintf(out, " %.16e", block_offdiag_im[i]);
    fprintf(out, "\n");
  };

  Timer timer;
  int n_meas = 0;
  for (int n = 0; n < (n_therm + n_traj); n++) {
    for (int j = 0; j < n_metropolis; j++) field.Metropolis();
    for (int j = 0; j < n_wolff; j++) field.WolffUpdate();

    if (n % n_skip || n < n_therm) continue;

    double M = 0.0;
    for (int s = 0; s < n_sites; s++) M += field.spin[s];
    M /= n_sites;
    double M2 = M * M;
    block_M += M;
    block_M2 += M2;
    block_M3 += M2 * M;
    block_M4 += M2 * M2;

    std::vector<std::complex<double>> S(n_lm);
    for (int k = 0; k < n_lm; k++) {
      std::complex<double> S_lm = 0.0;
      for (int s = 0; s < n_sites; s++) S_lm += wY[k][s] * field.spin[s];
      S[k] = S_lm;
      block_abs2[k] += std::norm(S_lm);  // |S_lm|^2
    }
    for (int i = 0; i < n_offdiag; i++) {
      std::complex<double> val = S[offdiag_pairs[i].k] * std::conj(S[offdiag_pairs[i].kp]);
      block_offdiag_re[i] += val.real();
      block_offdiag_im[i] += val.imag();
    }

    block_n++;
    n_meas++;
    if (block_n == jack_block_size) {
      write_block();
      block_index++;
      block_n = 0;
      block_M = block_M2 = block_M3 = block_M4 = 0.0;
      std::fill(block_abs2.begin(), block_abs2.end(), 0.0);
      std::fill(block_offdiag_re.begin(), block_offdiag_re.end(), 0.0);
      std::fill(block_offdiag_im.begin(), block_offdiag_im.end(), 0.0);
    }
  }
  // final partial block (only if it has any measurements -- jackknife code
  // downstream should treat blocks of unequal size correctly since n_meas
  // is written explicitly per block)
  if (block_n > 0) {
    write_block();
    block_index++;
  }
  fclose(out);
  timer.Stop();
  printf("wrote %d blocks (%d measurements) to %s (%.3fs)\n", block_index, n_meas, out_path.c_str(),
         timer.Duration());

  return 0;
}
