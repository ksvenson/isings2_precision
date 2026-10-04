// test_orbit_avg.cc
// Small-scale test of icosahedral-orbit averaging for the real-space
// two-point function, per user suggestion 2026-08-25: instead of binning
// pairs purely by z=cos(theta) (real_space_2pt_test.cc's method, which
// implicitly treats all pairs at the same z as statistically equivalent --
// only exactly true in the continuum SO(3) limit, not on the discrete
// icosahedral mesh), average <s_ref s_targ> over the exact |G|=120-element
// icosahedral point-group orbit of the (ref, targ) pair. This is an exact
// zero-extra-simulation-cost symmetrization: for the SAME measured
// configuration, apply every g in G to both site positions, look up the
// resulting site indices, and average the correlator over all group
// images. Two questions this answers: (1) does orbit-averaging noticeably
// reduce the per-config variance vs. a single pair or vs. n_ref random
// reference sites (variance-reduction value), and (2) does the orbit-
// averaged mean differ from the naive z-binned mean (a direct measurement
// of icosahedral-symmetry contamination in the existing real-space
// method)?
//
// Deliberately small-scale (per user request "do a smaller run to test")
// -- brute-force O(n_sites) position lookup per group element is fine at
// n_refine<=16 but would not scale to production n_refine (128) without a
// KD-tree; this is a feasibility/value test, not a production driver.
//
// Usage: test_orbit_avg --q 5 --n_refine 8 --coupling_rule exact_sinh
//   --n_therm 2000 --n_traj 5000 --n_skip 2 --seed 1234 --n_ref_sites 5
#include <getopt.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "ising.h"
#include "S2.h"
#include "statistics.h"
#include "timer.h"
#include "util.h"

// Find the site whose position most closely matches p (brute force).
int FindSiteByPosition(const QfeLatticeS2& lattice, const Vec3& p) {
  int best = -1;
  double best_dot = -2.0;
  for (int s = 0; s < lattice.n_sites; s++) {
    double d = lattice.r[s].normalized().dot(p.normalized());
    if (d > best_dot) {
      best_dot = d;
      best = s;
    }
  }
  return best;
}

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = 8;
  unsigned int seed = 1234u;
  int n_therm = 2000;
  int n_traj = 5000;
  int n_skip = 2;
  int n_wolff = 5;
  int n_metropolis = 4;
  int n_ref_sites = 5;
  std::string coupling_rule = "exact_sinh";

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"seed", required_argument, 0, 'S'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"n_ref_sites", required_argument, 0, 'r'},
      {"coupling_rule", required_argument, 0, 1000},
      {0, 0, 0, 0}};
  const char* short_options = "q:N:S:h:t:s:w:e:r:";
  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;
    switch (c) {
      case 'q': q = atoi(optarg); break;
      case 'N': n_refine = atoi(optarg); break;
      case 'S': seed = atol(optarg); break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      case 'r': n_ref_sites = atoi(optarg); break;
      case 1000: coupling_rule = optarg; break;
      default: break;
    }
  }

  QfeLatticeS2 lattice(q, n_refine);
  lattice.SeedRng(seed);
  lattice.UpdateWeights();
  printf("n_refine: %d q: %d total sites: %d |G|: %zu\n", n_refine, q,
         lattice.n_sites, lattice.G.size());

  QfeIsing field(&lattice, 1.0);
  field.HotStart();

  // coupling assignment (verbatim copy of real_space_2pt_test.cc /
  // ising_s2_crit.cc's link loop)
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
    double len_l0_sq = lattice.EdgeSquared(l0);
    int l1 = lattice.FindLink(s0, s2);
    int l2 = lattice.FindLink(s1, s2);
    int l3 = lattice.FindLink(s0, s3);
    int l4 = lattice.FindLink(s1, s3);
    double len_l1 = lattice.EdgeLength(l1);
    double len_l2 = lattice.EdgeLength(l2);
    double len_l3 = lattice.EdgeLength(l3);
    double len_l4 = lattice.EdgeLength(l4);
    double cos_theta_f0 = (len_l1 * len_l1 + len_l2 * len_l2 - len_l0_sq) / (2.0 * len_l1 * len_l2);
    double cos_theta_f1 = (len_l3 * len_l3 + len_l4 * len_l4 - len_l0_sq) / (2.0 * len_l3 * len_l4);
    double theta_f0 = acos(cos_theta_f0);
    double theta_f1 = acos(cos_theta_f1);
    double K_exact_sinh_f0 = 0.5 * asinh(1.0 / tan(theta_f0));
    double K_exact_sinh_f1 = 0.5 * asinh(1.0 / tan(theta_f1));
    lattice.links[l0].wt = 0.5 * (K_exact_sinh_f0 + K_exact_sinh_f1);
  }
  (void)coupling_rule;  // only exact_sinh wired for this quick test

  // pick n_ref_sites reference sites, and for each, one target site at a
  // moderate angular separation (roughly a third of the way to the
  // antipode) -- enough to be away from both short-distance contamination
  // and the antipodal special point.
  std::mt19937 pick_rng(seed + 7);
  std::vector<int> all_sites(lattice.n_sites);
  for (int s = 0; s < lattice.n_sites; s++) all_sites[s] = s;
  std::shuffle(all_sites.begin(), all_sites.end(), pick_rng);
  n_ref_sites = std::min(n_ref_sites, lattice.n_sites);

  struct PairTest {
    int ref, targ;
    std::vector<std::pair<int, int>> orbit_pairs;  // deduped (ref', targ') images
  };
  std::vector<PairTest> tests;
  for (int i = 0; i < n_ref_sites; i++) {
    int ref = all_sites[i];
    // target: the site whose position is closest to a fixed target
    // z = ref . target ~ 0.0 (roughly 90 degrees away)
    Vec3 rref = lattice.r[ref].normalized();
    Vec3 arbitrary(1.0, 0.0, 0.0);
    if (std::abs(rref.dot(arbitrary)) > 0.9) arbitrary = Vec3(0.0, 1.0, 0.0);
    Vec3 perp = (arbitrary - rref * rref.dot(arbitrary)).normalized();
    int targ = FindSiteByPosition(lattice, perp);

    PairTest t;
    t.ref = ref;
    t.targ = targ;
    std::vector<std::pair<int, int>> raw_pairs;
    for (size_t g = 0; g < lattice.G.size(); g++) {
      Vec3 pref = lattice.G[g] * rref;
      Vec3 ptarg = lattice.G[g] * lattice.r[targ].normalized();
      int sref = FindSiteByPosition(lattice, pref);
      int starg = FindSiteByPosition(lattice, ptarg);
      raw_pairs.push_back({sref, starg});
    }
    std::sort(raw_pairs.begin(), raw_pairs.end());
    raw_pairs.erase(std::unique(raw_pairs.begin(), raw_pairs.end()), raw_pairs.end());
    t.orbit_pairs = raw_pairs;
    tests.push_back(t);
    printf("ref=%d targ=%d z=%.6f orbit_size=%zu (|G|=%zu)\n", ref, targ,
           rref.dot(lattice.r[targ].normalized()), t.orbit_pairs.size(), lattice.G.size());
  }

  QfeMeasReal single_pair_meas;   // mean of the n_ref_sites raw pairs, per config
  QfeMeasReal orbit_avg_meas;     // mean of the n_ref_sites orbit-averaged pairs, per config

  Timer timer;
  for (int n = 0; n < (n_therm + n_traj); n++) {
    for (int j = 0; j < n_metropolis; j++) field.Metropolis();
    for (int j = 0; j < n_wolff; j++) field.WolffUpdate();

    if (n % n_skip || n < n_therm) continue;

    double single_sum = 0.0;
    double orbit_sum = 0.0;
    for (const PairTest& t : tests) {
      single_sum += field.spin[t.ref] * field.spin[t.targ];
      double acc = 0.0;
      for (const auto& p : t.orbit_pairs) acc += field.spin[p.first] * field.spin[p.second];
      orbit_sum += acc / double(t.orbit_pairs.size());
    }
    single_pair_meas.Measure(single_sum / tests.size());
    orbit_avg_meas.Measure(orbit_sum / tests.size());
  }
  timer.Stop();

  printf("\nmeasurement loop time: %.3fs (n_configs=%d)\n", timer.Duration(), single_pair_meas.n);
  printf("single-pair estimator (avg over %d random ref/targ pairs, no orbit averaging):\n", n_ref_sites);
  printf("  mean=%.6e err=%.6e\n", single_pair_meas.Mean(), single_pair_meas.Error());
  printf("orbit-averaged estimator (same configs, same pairs, exact icosahedral-orbit average):\n");
  printf("  mean=%.6e err=%.6e\n", orbit_avg_meas.Mean(), orbit_avg_meas.Error());
  printf("error ratio (single/orbit): %.3f  (sqrt(typical orbit size) ~ %.3f)\n",
         single_pair_meas.Error() / orbit_avg_meas.Error(),
         sqrt(double(tests[0].orbit_pairs.size())));
  double diff = orbit_avg_meas.Mean() - single_pair_meas.Mean();
  double diff_err = sqrt(pow(single_pair_meas.Error(), 2) + pow(orbit_avg_meas.Error(), 2));
  printf("mean difference (orbit - single): %.6e +/- %.6e (%.2f sigma)\n", diff, diff_err,
         diff / diff_err);

  return 0;
}
