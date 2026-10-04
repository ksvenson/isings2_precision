// fem_scalar_test.cc
//
// Deterministic (no Monte Carlo) diagnostic for IsingS2_precision's
// Delta_s(a) plateau (journal.md 2026-08-25 "ROOT CAUSE FOUND" entry and
// the 2026-08-25 follow-up question of whether the plateau is a bug in the
// shared mesh/quadrature/harmonic-projection code vs. a genuine physics
// effect specific to the Ising coupling assignment).
//
// Method: build the free massless scalar field's lattice propagator on
// the *same* QfeLatticeS2 mesh (same UpdateWeights cotangent-Laplacian
// weights, same UpdateYlm/GetYlm harmonics) this campaign's Ising driver
// uses, and project it onto spherical harmonics directly via one sparse
// linear solve per (l,m) -- no field configurations, no thermalization,
// no statistical error at all. This isolates whether the mesh/weight/
// harmonic-projection infrastructure itself breaks SO(3) symmetry or the
// known 1/(l(l+1)) scaling of the free-field 2-point function, independent
// of any question about whether the Ising coupling rule is exactly
// critical.
//
// Stiffness matrix: QfeLatticeS2::UpdateWeights() already computes, per
// link, links[l].wt = sum over the link's two adjacent faces of
// cot(opposite angle)/4 (see S2.h UpdateWeights, half_wt formula) -- this
// is exactly the standard cotangent-Laplacian / Regge FEM off-diagonal
// stiffness weight, and sites[s].wt is the matching lumped mass matrix
// (Voronoi dual area). No new geometry code needed: K_ii = sum_j w_ij,
// K_ij = -w_ij for neighbors j, M = diag(sites[s].wt).
//
// Regularization: K alone is singular (constant function is the zero
// mode). Solve (K + reg_eps*M) x = v instead; reg_eps is chosen far below
// the smallest physical eigenvalue l(l+1)=2 at l=1, so the induced
// eigenvalue shift is negligible and does not need to be corrected for.
//
// Continuum comparison: for a canonically normalized massless free scalar
// on the unit S^2, <phi(n)phi(n')> = sum_{l>=1} 1/(l(l+1)) * (2l+1)/(4pi)
// P_l(cos gamma) (l=0 mode removed/divergent, as usual). In the same
// M_l[m,m'] convention this campaign uses elsewhere (S_lm = integral phi
// Y_lm* dOmega), this means C_l = <S_lm S_lm*> should satisfy
// C_l * l(l+1) = const (l-independent) for l>=1, exactly analogous to how
// the Ising CFT test checks Delta rather than an absolute normalization.
#include <getopt.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "S2.h"
#include "util.h"

typedef std::complex<double> Complex;

int main(int argc, char* argv[]) {
  int q = 5;
  int n_refine = 4;
  int l_max = 8;
  double reg_eps = 1e-8;
  std::string mesh_mode = "naive";
  int equal_area_iters = 20000;
  double equal_area_step = 0.3;
  std::string mesh_cache_dir = "";

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"l_max", required_argument, 0, 'l'},
      {"reg_eps", required_argument, 0, 'r'},
      {"mesh_mode", required_argument, 0, 1002},
      {"equal_area_iters", required_argument, 0, 1003},
      {"equal_area_step", required_argument, 0, 1004},
      {"mesh_cache_dir", required_argument, 0, 1005},
      {0, 0, 0, 0},
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "q:N:l:r:", long_options, nullptr)) !=
         -1) {
    switch (opt) {
      case 'q':
        q = atoi(optarg);
        break;
      case 'N':
        n_refine = atoi(optarg);
        break;
      case 'l':
        l_max = atoi(optarg);
        break;
      case 'r':
        reg_eps = std::stod(optarg);
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
  printf("q: %d\n", q);
  printf("n_refine: %d\n", n_refine);
  printf("l_max: %d\n", l_max);
  printf("reg_eps: %.3e\n", reg_eps);

  QfeLatticeS2 lattice(q, n_refine);

  std::string mesh_mode_suffix = (mesh_mode == "equal_area")
                                      ? "_eqarea"
                                      : (mesh_mode == "equal_rp") ? "_eqrp" : "";

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
        printf("mesh_cache: %s %s\n",
               cache_hit ? "hit" : "stale/unreadable, recomputing",
               cache_path.c_str());
      } else {
        printf("mesh_cache: miss %s\n", cache_path.c_str());
      }
    }
    if (!cache_hit) {
      int iters_used = 0;
      if (mesh_mode == "equal_area") {
        lattice.EqualizeFaceAreas(equal_area_iters, equal_area_step,
                                   &iters_used);
      } else {
        lattice.EqualizeCircumPerim(equal_area_iters, equal_area_step,
                                     &iters_used);
      }
      printf("relax: iters_used=%d (cap=%d)\n", iters_used,
             equal_area_iters);
      if (!mesh_cache_dir.empty()) {
        FILE* cache_out = fopen(cache_path.c_str(), "wb");
        if (cache_out != nullptr) {
          lattice.WritePositions(cache_out);
          fclose(cache_out);
          printf("mesh_cache: wrote %s\n", cache_path.c_str());
        }
      }
    }
  }

  lattice.UpdateWeights();
  lattice.UpdateYlm(l_max);

  const int n = lattice.n_sites;
  printf("total sites: %d\n", n);

  // Assemble sparse stiffness K (cotangent Laplacian) and lumped mass M.
  std::vector<double> mass(n);
  for (int s = 0; s < n; s++) mass[s] = lattice.sites[s].wt;

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(size_t(n) * 8);
  std::vector<double> diag(n, 0.0);
  for (int l = 0; l < lattice.n_links; l++) {
    double w = lattice.links[l].wt;
    int a = lattice.links[l].sites[0];
    int b = lattice.links[l].sites[1];
    diag[a] += w;
    diag[b] += w;
    triplets.emplace_back(a, b, -w);
    triplets.emplace_back(b, a, -w);
  }
  for (int s = 0; s < n; s++) {
    triplets.emplace_back(s, s, diag[s] + reg_eps * mass[s]);
  }

  Eigen::SparseMatrix<double> A(n, n);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    fprintf(stderr, "sparse Cholesky factorization failed\n");
    return 1;
  }

  printf("%3s %4s %16s %16s %16s\n", "l", "dim", "C_l (avg diag)", "R_l (off-diag)",
         "C_l*l(l+1)");

  for (int l = 0; l <= l_max; l++) {
    int dim = 2 * l + 1;

    // v[m] = mass-weighted Y_lm at every site; x[m] = A^{-1} v[m].
    std::vector<Eigen::VectorXd> v_re(dim), v_im(dim), x_re(dim), x_im(dim);
    for (int mi = 0; mi < dim; mi++) {
      int m = mi - l;
      Eigen::VectorXd vr(n), vi(n);
      for (int s = 0; s < n; s++) {
        Complex y = lattice.GetYlm(s, l, m);
        vr(s) = mass[s] * y.real();
        vi(s) = mass[s] * y.imag();
      }
      v_re[mi] = vr;
      v_im[mi] = vi;
      x_re[mi] = solver.solve(vr);
      x_im[mi] = solver.solve(vi);
    }

    // M_l[m,m'] = v_m^dagger . x_m' (v^dagger since S_lm uses Y_lm*, and
    // x_m' = A^{-1} v_m' already carries the mass weight on the RHS side).
    std::vector<std::vector<Complex>> M(dim, std::vector<Complex>(dim));
    for (int mi = 0; mi < dim; mi++) {
      for (int mj = 0; mj < dim; mj++) {
        double re = v_re[mi].dot(x_re[mj]) + v_im[mi].dot(x_im[mj]);
        double im = v_re[mi].dot(x_im[mj]) - v_im[mi].dot(x_re[mj]);
        M[mi][mj] = Complex(re, im);
      }
    }

    double diag_sum = 0.0, diag_sq_sum = 0.0, off_sq_sum = 0.0;
    for (int mi = 0; mi < dim; mi++) {
      double d = M[mi][mi].real();
      diag_sum += d;
      diag_sq_sum += d * d;
      for (int mj = 0; mj < dim; mj++) {
        if (mi == mj) continue;
        off_sq_sum += std::norm(M[mi][mj]);
      }
    }
    double c_l = diag_sum / dim;
    double r_l = diag_sq_sum > 0.0 ? off_sq_sum / diag_sq_sum : 0.0;

    printf("%3d %4d %16.8e %16.8e %16.8e\n", l, dim, c_l, r_l,
           l > 0 ? c_l * l * (l + 1) : 0.0);
  }

  return 0;
}
