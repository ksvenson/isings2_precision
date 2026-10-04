// save_configs.cc
// Pure data-collection driver: builds the mesh, assigns couplings, runs
// Metropolis+Wolff, and dumps (1) site positions + quadrature weights once,
// (2) a sequence of packed raw spin configurations -- NO analysis of any
// kind baked in. Built 2026-08-25 per user direction ("let us just save
// the real-space correlators and then do the spherical symmetry tests by
// projecting afterwards") after real_space_2pt_test.cc's real-space fit
// found a clean Delta~0.128 while every harmonic (F_l/Delta_l_pair)
// estimator gave a resolution-independent wrong answer -- decouples
// expensive MC generation from analysis so any post-processing scheme
// (real-space fit with any window, exact M_l/F_l harmonic projection,
// UV-pair-excluded M_l, anything else) can be tried offline against the
// same saved ensemble without rerunning the simulation. See
// IsingS2_precision/journal.md 2026-08-25 "save_configs" entry.
//
// Reuses ising_s2_crit.cc's mesh-construction/coupling-assignment code
// verbatim (kept in sync by hand, mirrors real_space_2pt_test.cc's
// precedent).
//
// Output files (both under --data_dir, named by run_id/seed):
//   <run_id>_positions_<seed_hex>.dat -- text, one line per site:
//     "s x y z wt" (wt = QfeLatticeS2::UpdateWeights's quadrature weight)
//   <run_id>_configs_<seed_hex>.bin -- binary, n_saved back-to-back packed
//     spin configs, each ceil(n_sites/8) bytes, same bit-packing as
//     QfeIsing::WriteField/ReadField (bit i set => spin[i]=-1), no header
//     -- n_sites (needed to know the per-config byte stride) comes from
//     the positions file / command-line args, not stored in this file.
#include <getopt.h>

#include <Eigen/Dense>
#include <cmath>
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
  unsigned int seed = 1234u;
  bool cold_start = false;
  int n_therm = 2000;
  int n_traj = 20000;
  int n_skip = 2;
  int n_wolff = 5;
  int n_metropolis = 4;
  std::string data_dir = "save_configs";
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

  // dump site positions + quadrature weights once
  std::string pos_path = string_format("%s/%s_positions_%08X.dat", data_dir.c_str(), run_id.c_str(), seed);
  FILE* pos_file = fopen(pos_path.c_str(), "w");
  if (pos_file == nullptr) {
    fprintf(stderr, "could not open %s for writing\n", pos_path.c_str());
    return 1;
  }
  fprintf(pos_file, "# s x y z wt\n");
  for (int s = 0; s < lattice.n_sites; s++) {
    fprintf(pos_file, "%d %.16e %.16e %.16e %.16e\n", s, lattice.r[s].x(), lattice.r[s].y(),
            lattice.r[s].z(), lattice.sites[s].wt);
  }
  fclose(pos_file);
  printf("wrote positions: %s\n", pos_path.c_str());

  std::string configs_path = string_format("%s/%s_configs_%08X.bin", data_dir.c_str(), run_id.c_str(), seed);
  FILE* configs_file = fopen(configs_path.c_str(), "wb");
  if (configs_file == nullptr) {
    fprintf(stderr, "could not open %s for writing\n", configs_path.c_str());
    return 1;
  }

  Timer timer;
  int n_saved = 0;
  for (int n = 0; n < (n_therm + n_traj); n++) {
    for (int j = 0; j < n_metropolis; j++) field.Metropolis();
    for (int j = 0; j < n_wolff; j++) field.WolffUpdate();

    if (n % n_skip || n < n_therm) continue;

    field.WriteField(configs_file);
    n_saved++;
  }
  fclose(configs_file);
  timer.Stop();
  printf("wrote %d configs to %s (%.3fs)\n", n_saved, configs_path.c_str(), timer.Duration());

  return 0;
}
