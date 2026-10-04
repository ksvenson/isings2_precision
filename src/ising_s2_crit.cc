// ising_s2_crit.cc
// By Evan Owen 
#include <getopt.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
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
  // default parameters
  int q = 5;
  int n_refine = 1;
  unsigned int seed = 1234u;
  bool cold_start = false;
  int l_max = 6;
  int n_therm = 2000;
  int n_traj = 20000;
  int n_skip = 10;
  int n_wolff = 5;
  int n_metropolis = 4;
  double wall_time = 0.0;
  std::string data_dir = "ising_s2_crit";
  // Block size (in measured samples, i.e. post n_skip) for the M_l[m,m']
  // block-jackknife accumulator used by the spherical-symmetry test. Kept
  // configurable rather than hardcoded: Wolff-cluster updates make tau_int
  // for these observables plausibly O(1) sample, but that must be checked
  // empirically per run (e.g. by comparing jackknife error at a few block
  // sizes), not assumed. See IsingS2_precision/journal.md 2026-08-22.
  int jack_block_size = 1;
  // "duality": the pre-existing per-link formula below (dual-link-length
  // based, K=0.5*asinh(1/sinh(2L))). "exact_sinh": the flat-triangle
  // star-triangle rule K(theta)=0.5*asinh(cot(theta)) imported from Twist/,
  // applied per adjacent face and averaged across the link's two faces.
  // "owen_dual": sinh(2K_ij) = l*_ij / l_ij (Owen/Brower, arXiv:2503.05621
  // Eq 15) -- l_ij is the direct chord edge length, l*_ij is the straight-
  // line distance between the two adjacent faces' circumcenters (the dual
  // edge). Added 2026-08-25 after IsingS2_precision found exact_sinh's
  // Delta_s(a) plateaus at ~0.33 (not the exact 1/8) and does not shrink
  // with refinement -- neither duality nor exact_sinh implements this
  // paper's actual formula (duality uses an unrelated "kink" cross-ratio
  // construction; exact_sinh is a flat-lattice angle formula ported to the
  // curved mesh as an unverified judgment call). See
  // IsingS2_precision/journal.md 2026-08-25.
  // See IsingS2_precision/journal.md 2026-08-22 for why duality/exact_sinh
  // are two distinct, not-yet-reconciled formulas.
  std::string coupling_rule = "duality";
  // "naive": the constructor's plain flat-subdivide-then-radially-project
  // mesh, unmodified. "equal_area": run QfeLatticeS2::EqualizeFaceAreas
  // (gradient descent + backtracking line search) as a post-process before
  // weights/couplings are computed from site positions -- verified
  // 2026-08-25 to preserve icosahedral symmetry to machine precision and to
  // give the best area-uniformity of the mesh-construction options tried
  // (see journal.md 2026-08-25). Simplified 2026-08-25 (per user direction)
  // to just these two: the driver no longer exposes "equal_area_analytic"
  // (per-point Snyder-style placement, boundary-layer-only correction,
  // worse area-uniformity than equal_area at matched n_refine) or the
  // grp/orbit-group-element machinery it and --orbit_path relied on
  // (QfeLatticeS2::ReadSymmetryData/UpdateOrbits/ReadOrbits/G/site_g,
  // needlessly complicated for now per user direction) -- that code is
  // untouched in S2.h for reference, just never invoked from here (the
  // constructor's equal_area_analytic argument always defaults to false).
  // See journal.md 2026-08-24/2026-08-25 "equal-area mesh" entries --
  // motivated by the l>=3 spherical-symmetry-breaking plateau found not to
  // shrink with refinement on the naive mesh.
  // "equal_rp": added 2026-08-25 after finding neither `naive` nor
  // `equal_area` fix the Delta_s(a)/delta_l(a) scaling-dimension plateau
  // found via the Owen_section_D diagnostic (IsingS2_precision/
  // scripts/cft_symmetry_test.py) -- arXiv:2407.00459 Sec 3.1/3.3/Appendix
  // B derives that this file's `duality`/`owen_dual` coupling rules are
  // only exactly critical when every triangle has equal circumradius AND
  // equal perimeter (not equal area; the paper explicitly tested and
  // rejected an equal-area-only smoother for the same reason this
  // campaign's `equal_area` mode was found not to fix the plateau). Runs
  // QfeLatticeS2::EqualizeCircumPerim, minimizing the joint objective
  // E_R+E_P+E_A (area kept in the mix per user direction -- it was
  // already shown to help SO(3)-restoration speed, see journal.md
  // 2026-08-25 push7 entry; the paper's finding is that area uniformity
  // alone is insufficient, not that it hurts jointly). See
  // journal.md 2026-08-25 "ROOT CAUSE FOUND" entry.
  std::string mesh_mode = "naive";
  // Cap on GD iterations -- since EqualizeFaceAreas now has a practical
  // relative-improvement early stop (S2.h, added 2026-08-25; see
  // journal.md "iters-needed model" entry), this is a generous safety
  // ceiling, not the iteration count actually used at most resolutions.
  int equal_area_iters = 20000;
  double equal_area_step = 0.3;
  // If set, --mesh_mode equal_area first checks
  // <mesh_cache_dir>/q<q>_k<n_refine>_step<step>.dat for a previously
  // relaxed mesh at this exact (q, n_refine, step) and loads it instead of
  // re-running the (self-terminating but still nontrivial-cost, e.g.
  // O(minutes) at n_refine=128) gradient descent; on a cache miss, it
  // relaxes as usual and writes the result there for future runs. Empty
  // (default) disables caching entirely -- unchanged behavior. Not
  // data_dir/run-scoped: this is meant to be a single shared library
  // reused across every campaign run at a given resolution, not per-shard
  // output. See CLAUDE.md "Equal-area mesh position cache".
  std::string mesh_cache_dir = "";

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"seed", required_argument, 0, 'S'},
      {"cold_start", no_argument, 0, 'C'},
      {"l_max", required_argument, 0, 'l'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"data_dir", required_argument, 0, 'd'},
      {"wall_time", required_argument, 0, 'W'},
      {"coupling_rule", required_argument, 0, 1000},
      {"jack_block_size", required_argument, 0, 1001},
      {"mesh_mode", required_argument, 0, 1002},
      {"equal_area_iters", required_argument, 0, 1003},
      {"equal_area_step", required_argument, 0, 1004},
      {"mesh_cache_dir", required_argument, 0, 1005},
      {0, 0, 0, 0}};

  const char* short_options = "q:N:S:Cl:h:t:s:w:e:d:W:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'q':
        q = atoi(optarg);
        break;
      case 'N':
        n_refine = atoi(optarg);
        break;
      case 'S':
        seed = atol(optarg);
        break;
      case 'C':
        cold_start = true;
        break;
      case 'l':
        l_max = atoi(optarg);
        break;
      case 'h':
        n_therm = atoi(optarg);
        break;
      case 't':
        n_traj = atoi(optarg);
        break;
      case 's':
        n_skip = atoi(optarg);
        break;
      case 'w':
        n_wolff = atoi(optarg);
        break;
      case 'e':
        n_metropolis = atoi(optarg);
        break;
      case 'd':
        data_dir = optarg;
        break;
      case 'W':
        wall_time = std::stod(optarg);
        break;
      case 1000:
        coupling_rule = optarg;
        break;
      case 1001:
        jack_block_size = atoi(optarg);
        break;
      case 1002:
        mesh_mode = optarg;
        break;
      case 1003:
        equal_area_iters = atoi(optarg);
        break;
      case 1004:
        equal_area_step = std::stod(optarg);
        break;
      case 1005:
        mesh_cache_dir = optarg;
        break;
      default:
        break;
    }
  }

  if (mesh_mode != "naive" && mesh_mode != "equal_area" &&
      mesh_mode != "equal_rp") {
    fprintf(stderr,
            "unknown --mesh_mode '%s' (expected 'naive', 'equal_area', or "
            "'equal_rp')\n",
            mesh_mode.c_str());
    return 1;
  }
  printf("mesh_mode: %s\n", mesh_mode.c_str());

  if (coupling_rule != "duality" && coupling_rule != "exact_sinh" &&
      coupling_rule != "owen_dual") {
    fprintf(stderr,
            "unknown --coupling_rule '%s' (expected 'duality', 'exact_sinh', "
            "or 'owen_dual')\n",
            coupling_rule.c_str());
    return 1;
  }
  printf("coupling_rule: %s\n", coupling_rule.c_str());

  if (jack_block_size < 1) {
    fprintf(stderr, "--jack_block_size must be >= 1 (got %d)\n", jack_block_size);
    return 1;
  }
  printf("jack_block_size: %d\n", jack_block_size);

  std::string mesh_mode_suffix = (mesh_mode == "equal_area")
                                      ? "_eqarea"
                                      : (mesh_mode == "equal_rp") ? "_eqrp" : "";
  std::string run_id =
      string_format("q%dk%d%s", q, n_refine, mesh_mode_suffix.c_str());
  printf("run_id: %s\n", run_id.c_str());

  printf("n_therm: %d\n", n_therm);
  printf("n_traj: %d\n", n_traj);
  printf("n_skip: %d\n", n_skip);
  printf("n_wolff: %d\n", n_wolff);
  printf("n_metropolis: %d\n", n_metropolis);

  // number of spherical harmonics to measure
  int n_ylm = ((l_max + 1) * (l_max + 2)) / 2;
  printf("l_max: %d\n", l_max);
  printf("n_ylm: %d\n", n_ylm);

  // create a refined triangular lattice (equal_area_analytic and the
  // grp/orbit-group-element machinery are unused here -- see the mesh_mode
  // comment above)
  QfeLatticeS2 lattice(q, n_refine);

  // report face area/circumradius/perimeter non-uniformity, before and
  // after any relaxation, for comparison -- circumradius/perimeter added
  // 2026-08-25 (see mesh_mode "equal_rp" comment above) since those, not
  // area, are the quantities arXiv:2407.00459's coupling-rule derivation
  // actually requires to be uniform.
  auto print_mesh_uniformity = [&](const char* label) {
    double a_sum = 0.0, a_sq = 0.0, r_sum = 0.0, r_sq = 0.0, p_sum = 0.0, p_sq = 0.0;
    for (int f = 0; f < lattice.n_faces; f++) {
      double a = lattice.FlatArea(f);
      double rad = lattice.Circumradius(f);
      double per = lattice.Perimeter(f);
      a_sum += a; a_sq += a * a;
      r_sum += rad; r_sq += rad * rad;
      p_sum += per; p_sq += per * per;
    }
    double n = double(lattice.n_faces);
    double a_mean = a_sum / n, r_mean = r_sum / n, p_mean = p_sum / n;
    double a_std = sqrt(a_sq / n - a_mean * a_mean);
    double r_std = sqrt(r_sq / n - r_mean * r_mean);
    double p_std = sqrt(p_sq / n - p_mean * p_mean);
    printf("%s face area: mean=%.6e std/mean=%.6e | circumradius: "
           "mean=%.6e std/mean=%.6e | perimeter: mean=%.6e std/mean=%.6e\n",
           label, a_mean, a_std / a_mean, r_mean, r_std / r_mean, p_mean,
           p_std / p_mean);
  };
  print_mesh_uniformity("pre-relax");

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
        lattice.EqualizeFaceAreas(equal_area_iters, equal_area_step,
                                   &iters_used);
      } else {
        lattice.EqualizeCircumPerim(equal_area_iters, equal_area_step,
                                     &iters_used);
      }
      relax_timer.Stop();
      printf("%s: iters_used=%d (cap=%d) time=%.3fs\n", mesh_mode.c_str(),
             iters_used, equal_area_iters, relax_timer.Duration());
      if (!mesh_cache_dir.empty()) {
        FILE* cache_out = fopen(cache_path.c_str(), "wb");
        if (cache_out != nullptr) {
          lattice.WritePositions(cache_out);
          fclose(cache_out);
          printf("mesh_cache: wrote %s\n", cache_path.c_str());
        } else {
          fprintf(stderr,
                  "mesh_cache: could not open %s for writing (does "
                  "--mesh_cache_dir exist?), continuing without caching this "
                  "result\n",
                  cache_path.c_str());
        }
      }
    }
  }

  if (mesh_mode != "naive") {
    print_mesh_uniformity("post-relax");
  }

  lattice.SeedRng(seed);
  lattice.UpdateWeights();

  printf("n_refine: %d\n", n_refine);
  printf("q: %d\n", q);
  printf("total sites: %d\n", lattice.n_sites);

  double vol = lattice.vol;
  double vol_sq = vol * vol;

  QfeIsing field(&lattice, 1.0);

  // check if there is an rng file
  std::string rng_path =
      string_format("%s/%s/%s_rng_%08X.dat", data_dir.c_str(), run_id.c_str(),
                    run_id.c_str(), seed);
  printf("opening file: %s\n", rng_path.c_str());
  FILE* rng_file = fopen(rng_path.c_str(), "r");
  if (rng_file != nullptr) {
    printf("loading rng state from file\n");
    //lattice.rng.ReadRng(rng_file);
    fclose(rng_file);
  }

  // check if there is a field checkpoint file
  std::string field_path =
      string_format("%s/%s/%s_field_%08X.dat", data_dir.c_str(), run_id.c_str(),
                    run_id.c_str(), seed);
  printf("opening file: %s\n", field_path.c_str());
  FILE* field_file = fopen(field_path.c_str(), "rb");
  if (field_file != nullptr) {
    printf("loading field from checkpoint file\n");
    field.ReadField(field_file);
    fclose(field_file);
  } else if (cold_start) {
    printf("cold start\n");
    field.ColdStart();
  } else {
    printf("hot start\n");
    field.HotStart();
  }

  printf("initial action: %.12f\n", field.Action());

  // numbering for adjacent sites, links, and faces
  //
  //                      s2
  //                    /    \
  //        f1a       /        \       f2a
  //               l1            l2
  //              /       f0       \
  //            /                    \
  //          s0 - - - -  l0  - - - - s1
  //            \                    /
  //              \       f1       /
  //               l3            l4
  //        f3b       \        /       f4b
  //                    \    /
  //                      s3
  //
  //

  // 3. this way uses dual links on the trianglular lattice with a kink
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

    // calculate unit vectors from links to/from face circumcenters
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

    // Exact flat-triangle star-triangle rule (Twist/'s
    // K(theta)=0.5*asinh(cot(theta))), applied per adjacent face using that
    // face's angle opposite this link (law of cosines from the face's own
    // three chord lengths, treating each face as a flat Euclidean triangle),
    // then averaged over the link's two faces since each face independently
    // predicts a K for the shared edge and there is no single-triangle
    // convention on a curved mesh where every edge borders two distinct
    // faces (unlike Twist's flat lattice, where translational symmetry
    // makes the two faces equivalent). See journal.md 2026-08-22.
    double cos_theta_f0 = (len_l1 * len_l1 + len_l2 * len_l2 - len_l0_sq) /
                           (2.0 * len_l1 * len_l2);
    double cos_theta_f1 = (len_l3 * len_l3 + len_l4 * len_l4 - len_l0_sq) /
                           (2.0 * len_l3 * len_l4);
    double theta_f0 = acos(cos_theta_f0);
    double theta_f1 = acos(cos_theta_f1);
    double K_exact_sinh_f0 = 0.5 * asinh(1.0 / tan(theta_f0));
    double K_exact_sinh_f1 = 0.5 * asinh(1.0 / tan(theta_f1));
    double K_exact_sinh = 0.5 * (K_exact_sinh_f0 + K_exact_sinh_f1);

    // Owen/Brower arXiv:2503.05621 Eq 15: sinh(2K_ij) = l*_ij/l_ij, with
    // l*_ij the straight-line (chord) distance between the two adjacent
    // faces' circumcenters (v_f0, v_f1 above) -- the actual dual edge, not
    // the "kink" cross-ratio construction K_duality uses.
    double len_dual = (v_f0 - v_f1).norm();
    double K_owen_dual = 0.5 * asinh(len_dual / len_l0);

    double K = (coupling_rule == "exact_sinh")
                   ? K_exact_sinh
                   : (coupling_rule == "owen_dual") ? K_owen_dual : K_duality;

    // printf("%.12f %.12f %.12f\n", lattice.links[l0].wt, K_duality, K_exact_sinh);
    lattice.links[l0].wt = K;
  }

  // calculate spherical harmonics at each site. In all of these loops, (y_l,
  // y_m) are the integer eigenvalues of the spherical harmonics and y_i is a
  // sequential index over all of eigenvalues up to l_max (excluding negative
  // y_m).
  Eigen::MatrixXcd ylm(lattice.n_sites, n_ylm);
  for (int s = 0; s < lattice.n_sites; s++) {
    for (int y_i = 0, y_l = 0, y_m = 0; y_i < n_ylm; y_i++) {
      ylm(s, y_i) = lattice.CalcYlm(s, y_l, y_m);

      y_m++;
      if (y_m > y_l) {
        y_l++;
        y_m = 0;
      }
    }
  }

  // measurements
  std::vector<QfeMeasReal> legendre_2pt(l_max + 1);
  std::vector<QfeMeasReal> ylm_2pt(n_ylm);

  // Full (2l+1)x(2l+1) M_l[m,m'] = <S_lm S*_lm'> matrix (upper triangle
  // m<=m' only; M_l is Hermitian, so this determines the rest) for the
  // spherical-symmetry test. ylm_2pt above only ever measured m=m'>=0, so
  // it can't detect off-diagonal SO(3) breaking. S_l,-m for negative m is
  // reconstructed from S_l,m via the same (-1)^m*conj() relation GetYlm
  // uses for Y_l,-m (S2.h), not separately computed.
  struct FullPair { int l, m, mp; };
  std::vector<FullPair> full_pairs;
  for (int l = 0; l <= l_max; l++) {
    for (int m = -l; m <= l; m++) {
      for (int mp = m; mp <= l; mp++) {
        full_pairs.push_back({l, m, mp});
      }
    }
  }
  std::vector<QfeMeasComplex> ylm_2pt_full(full_pairs.size());

  // Block-jackknife accumulator for the symmetry-test statistics computed
  // offline (R_l off-diagonal power ratio, chi-square of the diagonal
  // across m) — see IsingS2_precision/journal.md 2026-08-22. Stores only
  // completed blocks of jack_block_size measured samples each (plus a
  // final partial block flushed at the end); each block's per-pair sums
  // are synchronized across (l,m,mp), which is what lets the offline
  // jackknife capture cross-entry correlations that per-entry errors miss.
  struct JackBlock {
    int n;
    std::vector<Complex> sum;  // one entry per full_pairs[i]
  };
  std::vector<JackBlock> jack_blocks;
  std::vector<Complex> jack_block_cur_sum(full_pairs.size(), 0.0);
  int jack_block_cur_n = 0;

  QfeMeasReal mag;     // magnetization
  QfeMeasReal mag_2;   // magnetization^2
  QfeMeasReal mag_4;   // magnetization^4
  QfeMeasReal mag_6;   // magnetization^6
  QfeMeasReal mag_8;   // magnetization^8
  QfeMeasReal mag_10;  // magnetization^10
  QfeMeasReal mag_12;  // magnetization^12
  QfeMeasReal action;
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;


  //Save Lattice Configuration //Commented out by JYL, not needed to save this for now


  // char config_path[200];
  
  // sprintf(config_path, "%s/%s/%s_%08X_config.dat", \
  //     data_dir.c_str(), run_id.c_str(), run_id.c_str(), seed);
  // FILE* configfile = fopen(config_path, "w");
  // printf("opening file: %s\n", config_path);
  // assert(configfile!= nullptr);

  // for (int l0 = 0; l0<lattice.n_links; l0++){
  //   int s0 = lattice.links[l0].sites[0];
  //   int s1 = lattice.links[l0].sites[1];

  //   fprintf(configfile, "%d %.12f %.12f %.12f %.12f %.12f %.12f\n", l0, lattice.r[s0].x(), \
  //   lattice.r[s0].y(), lattice.r[s0].z(), lattice.r[s1].x(), lattice.r[s1].y(), lattice.r[s1].z());
  // }
  


  // read bulk measurements
  std::string bulk_path =
      string_format("%s/%s/%s_bulk_%08X.dat", data_dir.c_str(), run_id.c_str(),
                    run_id.c_str(), seed);
  printf("opening file: %s\n", bulk_path.c_str());
  FILE* bulk_file = fopen(bulk_path.c_str(), "r");
  if (bulk_file != nullptr) {
    while (!feof(bulk_file)) {
      char meas_name[40];
      fscanf(bulk_file, "%s ", meas_name);
      if (strcmp(meas_name, "action") == 0) {
        action.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag") == 0) {
        mag.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag") == 0) {
        mag.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^2") == 0) {
        mag_2.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^4") == 0) {
        mag_4.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^6") == 0) {
        mag_6.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^8") == 0) {
        mag_8.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^10") == 0) {
        mag_10.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^12") == 0) {
        mag_12.ReadMeasurement(bulk_file);
      } else {
        printf("unknown measurement: %s\n", meas_name);
      }
    }
    fclose(bulk_file);
  }

  // read 2-point function legendre coefficients
  std::string legendre_2pt_path =
      string_format("%s/%s/%s_legendre_2pt_%08X.dat", data_dir.c_str(),
                    run_id.c_str(), run_id.c_str(), seed);
  printf("opening file: %s\n", legendre_2pt_path.c_str());
  FILE* legendre_2pt_file = fopen(legendre_2pt_path.c_str(), "r");
  if (legendre_2pt_file != nullptr) {
    int l;
    while (!feof(legendre_2pt_file)) {
      fscanf(legendre_2pt_file, "%d ", &l);
      if (l > l_max) {
        fscanf(legendre_2pt_file, "%*[^\n]");
      } else {
        legendre_2pt[l].ReadMeasurement(legendre_2pt_file);
      }
      fscanf(legendre_2pt_file, "\n");
    }
    fclose(legendre_2pt_file);
  }

  // read 2-point function ylm coefficients
  std::string ylm_2pt_path =
      string_format("%s/%s/%s_ylm_2pt_%08X.dat", data_dir.c_str(),
                    run_id.c_str(), run_id.c_str(), seed);
  printf("opening file: %s\n", ylm_2pt_path.c_str());
  FILE* ylm_2pt_file = fopen(ylm_2pt_path.c_str(), "r");
  if (ylm_2pt_file != nullptr) {
    int i_ylm, l, m;
    while (!feof(ylm_2pt_file)) {
      fscanf(ylm_2pt_file, "%d %d %d ", &i_ylm, &l, &m);
      if (i_ylm >= n_ylm) {
        fscanf(ylm_2pt_file, "%*[^\n]");
      } else {
        ylm_2pt[i_ylm].ReadMeasurement(ylm_2pt_file);
      }
      fscanf(ylm_2pt_file, "\n");
    }
    fclose(ylm_2pt_file);
  }

  // read full off-diagonal M_l[m,m'] checkpoint, if present
  std::string ylm_2pt_full_path =
      string_format("%s/%s/%s_ylm_2pt_full_%08X.dat", data_dir.c_str(),
                    run_id.c_str(), run_id.c_str(), seed);
  printf("opening file: %s\n", ylm_2pt_full_path.c_str());
  FILE* ylm_2pt_full_file = fopen(ylm_2pt_full_path.c_str(), "r");
  if (ylm_2pt_full_file != nullptr) {
    for (size_t i = 0; i < full_pairs.size(); i++) {
      int l, m, mp;
      double mean_re, mean_im, err_re, err_im;
      int meas_n;
      int n_read = fscanf(ylm_2pt_full_file, "%d %d %d %lf %lf %lf %lf %d\n",
                           &l, &m, &mp, &mean_re, &mean_im, &err_re, &err_im,
                           &meas_n);
      if (n_read != 8) break;
      double count = double(meas_n);
      ylm_2pt_full[i].real_part.n = meas_n;
      ylm_2pt_full[i].real_part.sum = mean_re * count;
      ylm_2pt_full[i].real_part.sum2 = (err_re * err_re * count + mean_re * mean_re) * count;
      ylm_2pt_full[i].imag_part.n = meas_n;
      ylm_2pt_full[i].imag_part.sum = mean_im * count;
      ylm_2pt_full[i].imag_part.sum2 = (err_im * err_im * count + mean_im * mean_im) * count;
    }
    fclose(ylm_2pt_full_file);
  }

  // read completed block-jackknife blocks, if present. Resuming a run
  // assumes the same --jack_block_size and --l_max as the run that wrote
  // this file; changing either between resumes silently produces a mix of
  // block sizes / pair sets, which the offline jackknife does not detect —
  // keep --jack_block_size/--l_max fixed for a given data_dir/run_id/seed.
  std::string jack_blocks_path =
      string_format("%s/%s/%s_ylm_2pt_full_jackblocks_%08X.dat",
                    data_dir.c_str(), run_id.c_str(), run_id.c_str(), seed);
  printf("opening file: %s\n", jack_blocks_path.c_str());
  FILE* jack_blocks_file = fopen(jack_blocks_path.c_str(), "r");
  if (jack_blocks_file != nullptr) {
    int n_pairs_prev = 0, block_size_prev = 0, n_blocks_prev = 0;
    int n_read = fscanf(jack_blocks_file, "# n_pairs=%d jack_block_size=%d n_blocks=%d\n",
                         &n_pairs_prev, &block_size_prev, &n_blocks_prev);
    if (n_read == 3 && n_pairs_prev == (int)full_pairs.size()) {
      for (int k = 0; k < n_blocks_prev; k++) {
        JackBlock blk;
        blk.sum.assign(full_pairs.size(), 0.0);
        for (size_t i = 0; i < full_pairs.size(); i++) {
          int k_read, l, m, mp, n_k;
          double sum_re, sum_im;
          int n_fields = fscanf(jack_blocks_file, "%d %d %d %d %lf %lf %d\n",
                                 &k_read, &l, &m, &mp, &sum_re, &sum_im, &n_k);
          if (n_fields != 7) break;
          blk.sum[i] = Complex(sum_re, sum_im);
          blk.n = n_k;
        }
        jack_blocks.push_back(blk);
      }
    } else {
      fprintf(stderr,
              "warning: %s pair count/header mismatch, ignoring stale jackblocks file\n",
              jack_blocks_path.c_str());
    }
    fclose(jack_blocks_file);
  }

  Timer timer;

  for (int n = 0; n < (n_traj + n_therm); n++) {
    if (wall_time > 0.0 && timer.Duration() > wall_time) break;

    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) {
      metropolis_sum += field.Metropolis();
    }
    accept_metropolis.Measure(metropolis_sum);

    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) {
      cluster_size_sum += field.WolffUpdate();
    }
    cluster_size.Measure(double(cluster_size_sum) / vol);

    if (n % n_skip || n < n_therm) continue;

    // measure correlators
    std::vector<Complex> ylm_2pt_sum(n_ylm, 0.0);
    std::vector<Complex> ylm_4pt_sum(n_ylm, 0.0);
    double spin_sum = 0.0;

    for (int s = 0; s < lattice.n_sites; s++) {
      double wt_2pt = field.spin[s] * lattice.sites[s].wt;

      spin_sum += wt_2pt;

      for (int ylm_i = 0; ylm_i < n_ylm; ylm_i++) {
        Complex y = ylm(s, ylm_i);
        ylm_2pt_sum[ylm_i] += y * wt_2pt;
      }
    }

    // S_l,m for m<0, reconstructed from the m>=0 values above (same
    // relation as GetYlm's Y_l,-m = (-1)^m * conj(Y_l,m)).
    auto Slm = [&](int l, int m) -> Complex {
      int am = abs(m);
      Complex s = ylm_2pt_sum[l * (l + 1) / 2 + am];
      if (m < 0) {
        s = std::conj(s);
        if (am & 1) s *= -1.0;
      }
      return s;
    };
    for (size_t i = 0; i < full_pairs.size(); i++) {
      const FullPair& p = full_pairs[i];
      Complex val = Slm(p.l, p.m) * std::conj(Slm(p.l, p.mp)) / vol_sq;
      ylm_2pt_full[i].Measure(val);
      jack_block_cur_sum[i] += val;
    }
    jack_block_cur_n++;
    if (jack_block_cur_n == jack_block_size) {
      jack_blocks.push_back({jack_block_cur_n, jack_block_cur_sum});
      std::fill(jack_block_cur_sum.begin(), jack_block_cur_sum.end(), Complex(0.0, 0.0));
      jack_block_cur_n = 0;
    }

    double legendre_2pt_sum = 0.0;
    for (int ylm_i = 0, l = 0, m = 0; ylm_i < n_ylm; ylm_i++) {
      ylm_2pt[ylm_i].Measure(std::norm(ylm_2pt_sum[ylm_i]) / vol_sq);

      legendre_2pt_sum += ylm_2pt[ylm_i].last * (m == 0 ? 1.0 : 2.0);

      m++;
      if (m > l) {
        double coeff = 4.0 * M_PI / double(2 * l + 1);
        legendre_2pt[l].Measure(legendre_2pt_sum * coeff);
        legendre_2pt_sum = 0.0;
        l++;
        m = 0;
      }
    }

    double m = spin_sum / vol;
    double m_sq = m * m;
    mag.Measure(fabs(m));
    mag_2.Measure(m_sq);
    mag_4.Measure(m_sq * m_sq);
    mag_6.Measure(mag_4.last * m_sq);
    mag_8.Measure(mag_6.last * m_sq);
    mag_10.Measure(mag_8.last * m_sq);
    mag_12.Measure(mag_10.last * m_sq);
    action.Measure(field.Action());
    // printf("%06d %.12f %.4f %.4f\n", \
    //     n, action.last, \
    //     accept_metropolis.last, \
    //     cluster_size.last);
  }

  // flush a trailing partial block (n < jack_block_size) rather than
  // dropping those samples; the offline jackknife weights blocks by their
  // own n, so a smaller final block is handled correctly, not specially.
  if (jack_block_cur_n > 0) {
    jack_blocks.push_back({jack_block_cur_n, jack_block_cur_sum});
  }

  timer.Stop();
  printf("duration: %.6f\n", timer.Duration());

  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());

  // write rng state to file //Commented out by JYL, not needed to save this for now
  // printf("writing rng state to file: %s\n", rng_path.c_str());
  // rng_file = fopen(rng_path.c_str(), "w");
  // assert(rng_file != nullptr);
  // lattice.rng.WriteRng(rng_file);
  // fclose(rng_file);

  // // write field configuration to file
  // printf("writing field to file: %s\n", field_path.c_str());
  // field_file = fopen(field_path.c_str(), "wb");
  // assert(field_file != nullptr);
  // field.WriteField(field_file);
  // fclose(field_file);

  double m_mean = mag.Mean();
  double m_err = mag.Error();
  double m2_mean = mag_2.Mean();
  double m2_err = mag_2.Error();
  double m4_mean = mag_4.Mean();
  double m4_err = mag_4.Error();

  // open an output file
  printf("opening file: %s\n", bulk_path.c_str());
  bulk_file = fopen(bulk_path.c_str(), "w");
  assert(bulk_file != nullptr);

  printf("action: %+.12e %.12e %.4f %.4f\n", action.Mean(), action.Error(),
         action.AutocorrFront(), action.AutocorrBack());
  fprintf(bulk_file, "action ");
  action.WriteMeasurement(bulk_file);

  printf("mag: %.12e %.12e %.4f %.4f\n", m_mean, m_err, mag.AutocorrFront(),
         mag.AutocorrBack());
  fprintf(bulk_file, "mag ");
  mag.WriteMeasurement(bulk_file);

  printf("m^2: %.12e %.12e %.4f %.4f\n", m2_mean, m2_err, mag_2.AutocorrFront(),
         mag_2.AutocorrBack());
  fprintf(bulk_file, "mag^2 ");
  mag_2.WriteMeasurement(bulk_file);

  printf("m^4: %.12e %.12e %.4f %.4f\n", m4_mean, m4_err, mag_4.AutocorrFront(),
         mag_4.AutocorrBack());
  fprintf(bulk_file, "mag^4 ");
  mag_4.WriteMeasurement(bulk_file);

  printf("m^6: %.12e %.12e %.4f %.4f\n", mag_6.Mean(), mag_6.Error(),
         mag_6.AutocorrFront(), mag_6.AutocorrBack());
  fprintf(bulk_file, "mag^6 ");
  mag_6.WriteMeasurement(bulk_file);

  printf("m^8: %.12e %.12e %.4f %.4f\n", mag_8.Mean(), mag_8.Error(),
         mag_8.AutocorrFront(), mag_8.AutocorrBack());
  fprintf(bulk_file, "mag^8 ");
  mag_8.WriteMeasurement(bulk_file);

  printf("m^10: %.12e %.12e %.4f %.4f\n", mag_10.Mean(), mag_10.Error(),
         mag_10.AutocorrFront(), mag_10.AutocorrBack());
  fprintf(bulk_file, "mag^10 ");
  mag_10.WriteMeasurement(bulk_file);

  printf("m^12: %.12e %.12e %.4f %.4f\n", mag_12.Mean(), mag_12.Error(),
         mag_12.AutocorrFront(), mag_12.AutocorrBack());
  fprintf(bulk_file, "mag^12 ");
  mag_12.WriteMeasurement(bulk_file);

  double U4_mean = 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));
  double U4_err =
      0.5 * U4_mean *
      sqrt(pow(m4_err / m4_mean, 2.0) + pow(2.0 * m2_err / m2_mean, 2.0));
  printf("U4: %.12e %.12e\n", U4_mean, U4_err);

  double m_susc_mean = (m2_mean - m_mean * m_mean) * vol;
  double m_susc_err =
      sqrt(pow(m2_err, 2.0) + pow(2.0 * m_mean * m_err, 2.0)) * vol;
  printf("m_susc: %.12e %.12e\n", m_susc_mean, m_susc_err);

  // print 2-point function legendre coefficients
  printf("opening file: %s\n", legendre_2pt_path.c_str());
  legendre_2pt_file = fopen(legendre_2pt_path.c_str(), "w");
  assert(legendre_2pt_file != nullptr);
  for (int l = 0; l <= l_max; l++) {
    fprintf(legendre_2pt_file, "%02d ", l);
    legendre_2pt[l].WriteMeasurement(legendre_2pt_file);
  }
  fclose(legendre_2pt_file);

  // print 2-point function spherical harmonic coefficients
  printf("opening file: %s\n", ylm_2pt_path.c_str());
  ylm_2pt_file = fopen(ylm_2pt_path.c_str(), "w");
  assert(ylm_2pt_file != nullptr);
  for (int ylm_i = 0, l = 0, m = 0; ylm_i < n_ylm; ylm_i++) {
    fprintf(ylm_2pt_file, "%04d %02d %02d ", ylm_i, l, m);
    ylm_2pt[ylm_i].WriteMeasurement(ylm_2pt_file);
    m++;
    if (m > l) {
      l++;
      m = 0;
    }
  }
  fclose(ylm_2pt_file);

  // print full off-diagonal M_l[m,m'] matrix (symmetry-test observable)
  printf("opening file: %s\n", ylm_2pt_full_path.c_str());
  ylm_2pt_full_file = fopen(ylm_2pt_full_path.c_str(), "w");
  assert(ylm_2pt_full_file != nullptr);
  for (size_t i = 0; i < full_pairs.size(); i++) {
    const FullPair& p = full_pairs[i];
    Complex mean = ylm_2pt_full[i].Mean();
    Complex err = ylm_2pt_full[i].Error();
    fprintf(ylm_2pt_full_file, "%02d %+03d %+03d %.16e %.16e %.16e %.16e %d\n",
            p.l, p.m, p.mp, real(mean), imag(mean), real(err), imag(err),
            ylm_2pt_full[i].real_part.n);
  }
  fclose(ylm_2pt_full_file);

  // print block-jackknife blocks for the offline symmetry-test analysis
  // (scripts/symmetry_test.py in IsingS2_precision/ — see that campaign's
  // journal.md 2026-08-22 for the design)
  printf("opening file: %s\n", jack_blocks_path.c_str());
  jack_blocks_file = fopen(jack_blocks_path.c_str(), "w");
  assert(jack_blocks_file != nullptr);
  fprintf(jack_blocks_file, "# n_pairs=%zu jack_block_size=%d n_blocks=%zu\n",
          full_pairs.size(), jack_block_size, jack_blocks.size());
  for (size_t k = 0; k < jack_blocks.size(); k++) {
    const JackBlock& blk = jack_blocks[k];
    for (size_t i = 0; i < full_pairs.size(); i++) {
      const FullPair& p = full_pairs[i];
      fprintf(jack_blocks_file, "%zu %02d %+03d %+03d %.16e %.16e %d\n", k, p.l,
              p.m, p.mp, real(blk.sum[i]), imag(blk.sum[i]), blk.n);
    }
  }
  fclose(jack_blocks_file);

  return 0;
}