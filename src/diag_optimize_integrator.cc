// Ad hoc diagnostic (not part of any campaign pipeline): checks whether
// QfeLatticeS2::OptimizeIntegrator's least-squares solve actually moves
// away from its uniform wt=1.0 initial guess, or converges trivially.
// Build: see cluster/build.sh for the flags, swap in this source file.
#include <cstdio>
#include <string>
#include "S2.h"

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

int main(int argc, char** argv) {
  int q = 5;
  int n_refine = argc > 1 ? atoi(argv[1]) : 4;
  int l_max = argc > 2 ? atoi(argv[2]) : 8;
  const char* positions_path = argc > 3 ? argv[3] : nullptr;

  QfeLatticeS2 lattice(q, n_refine);
  if (positions_path) {
    if (!LoadTextPositions(positions_path, &lattice)) {
      printf("could not load %s\n", positions_path);
      return 1;
    }
    printf("loaded positions from %s\n", positions_path);
  }
  printf("n_sites=%d n_distinct=%d\n", lattice.n_sites, lattice.n_distinct);

  // replicate OptimizeIntegrator's setup here to print diagnostics
  std::vector<int> l_relevant, m_relevant;
  int m_spacing = (q == 5) ? 5 : 4;
  for (int l = 0; l <= l_max; l++) {
    if ((l % 2) && (q != 3)) continue;
    int n_l = (l / 2) + (l / 3) + (l / q) - l + 1;
    for (int i = 0; i < n_l; i++) {
      l_relevant.push_back(l);
      int m = i * m_spacing;
      if (l % 2) m = ((2 * i + 1) * m_spacing) / 2;
      m_relevant.push_back(m);
    }
  }
  int n_relevant = l_relevant.size();
  printf("n_relevant=%d: ", n_relevant);
  for (int i = 0; i < n_relevant; i++) printf("(l=%d,m=%d) ", l_relevant[i], m_relevant[i]);
  printf("\n");

  ComplexMat S = ComplexMat::Zero(n_relevant, lattice.n_distinct);
  for (int s = 0; s < lattice.n_sites; s++) {
    int id = lattice.sites[s].id;
    double theta = acos(lattice.r[s].z());
    double phi = atan2(lattice.r[s].y(), lattice.r[s].x());
    for (int i = 0; i < n_relevant; i++) {
      S(i, id) += boost::math::spherical_harmonic(l_relevant[i], m_relevant[i], theta, phi);
    }
  }

  ComplexVec x0(lattice.n_distinct);
  for (int id = 0; id < lattice.n_distinct; id++) {
    int s = lattice.distinct_first[id];
    x0(id) = lattice.sites[s].wt;
    if (id < 5) printf("x0(%d)=%.6f (from site %d)\n", id, real(x0(id)), s);
  }

  ComplexVec b = ComplexVec::Zero(n_relevant);
  b(0) = 0.28209479177387814347 * lattice.vol;
  printf("vol=%.6f n_sites=%d b(0)=%.6f\n", lattice.vol, lattice.n_sites, real(b(0)));

  ComplexVec resid0 = S * x0 - b;
  printf("||S*x0 - b|| = %.6e (per-row: ", resid0.norm());
  for (int i = 0; i < n_relevant; i++) printf("%.6e ", std::abs(resid0(i)));
  printf(")\n");

  Eigen::LeastSquaresConjugateGradient<ComplexMat> cg;
  cg.compute(S);
  printf("cg.info() after compute() = %d (0=Success)\n", (int)cg.info());
  ComplexVec wt = cg.solveWithGuess(b, x0);
  printf("cg.info() after solve = %d, iterations=%d, error=%.6e\n",
         (int)cg.info(), (int)cg.iterations(), cg.error());

  ComplexVec resid1 = S * wt - b;
  printf("||S*wt - b|| = %.6e\n", resid1.norm());

  double max_dev = 0;
  for (int id = 0; id < lattice.n_distinct; id++) {
    double dev = std::abs(real(wt(id)) - real(x0(id)));
    if (dev > max_dev) max_dev = dev;
  }
  printf("max|wt - x0| over all distinct ids = %.6e\n", max_dev);
  for (int id = 0; id < std::min(5, lattice.n_distinct); id++) {
    printf("  id=%d: x0=%.6f wt=%.6f\n", id, real(x0(id)), real(wt(id)));
  }

  if (argc > 4) {
    FILE* out = fopen(argv[4], "w");
    fprintf(out, "# site opt_wt (from OptimizeIntegrator(l_max=%d))\n", l_max);
    for (int s = 0; s < lattice.n_sites; s++) {
      int id = lattice.sites[s].id;
      fprintf(out, "%d %.15e\n", s, real(wt(id)));
    }
    fclose(out);
    printf("wrote %s\n", argv[4]);
  }
  return 0;
}
