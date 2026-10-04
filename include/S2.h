// s2.h
// Mostly by Evan Owen, but maybe minor additions by JYL
#pragma once

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "grp_o3.h"
#include "lattice.h"
#include "util.h"

// symmetry group data directory must be defined
#ifndef GRP_DIR
#error Error: GRP_DIR is not defined (it should be in the Makefile)
#endif

typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> ComplexMat;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> ComplexVec;

/// @brief Simplicial lattice discretization of a 2-sphere
class QfeLatticeS2 : public QfeLattice {
 public:
  QfeLatticeS2(int q = 5, int k = 1, bool equal_area_analytic = false);
  Vec3 EqualAreaGnomonicPoint(Vec3 face_r[3], int x, int y, int k);
  void ReadBaseLattice(int q);
  void WriteSite(FILE* file, int s);
  void ReadSite(FILE* file, int s);
  int CreateOrbit(double xi1, double xi2);
  Vec3 CalcOrbitPos(int o);
  void ReadOrbits(FILE* file);
  void WriteOrbits(FILE* file);
  void UpdateOrbits();
  void ReadSymmetryData(int q, int k);
  void ApplyEqualAreaAnalytic(int k);
  void ResizeSites(int n_sites);
  void Inflate();
  void UpdateAntipodes();
  Vec3 FaceCircumcenter(int f);
  double EdgeSquared(int l);
  double EdgeLength(int l);
  double FlatArea(int f);
  double Perimeter(int f);
  void EqualizeFaceAreas(int n_iter, double step, int* iters_used = nullptr,
                          double rel_tol = 1e-5, int patience = 20);
  void EqualizeDualAreas(int n_iter, double step);
  double Circumradius(int f);
  void EqualizeCircumPerim(int n_iter, double step, int* iters_used = nullptr,
                            double rel_tol = 1e-5, int patience = 20);
  void WritePositions(FILE* file);
  bool ReadPositions(FILE* file);


  //double DeficitAngle(int s);
  //double DualArea(int s);
  void UpdateWeights();
  double CalcLatticeSpacing();
  void OptimizeIntegrator(int l_max);
  void UpdateYlm(int l_max);
  Complex GetYlm(int s, int l, int m);
  Complex CalcYlm(int s, int l, int m);
  double CosTheta(int s1, int s2);
  void PrintCoordinates();

  int q;                        // base polyhedron parameter
  std::vector<Vec3> r;          // vertex coordinates
  std::vector<int> antipode;    // antipode of each site (0 by default)
  std::vector<int> site_orbit;  // orbit id for each site
  std::vector<int> face_orbit;  // orbit id for each face
  std::vector<Vec3> orbit_xi;   // barycentric coordinates for each orbit
  Vec3 first_face_r[4];         // coordinates of first face vertices
  std::vector<GrpElemO3> G;     // symmetry group elements
  std::vector<int> site_g;      // group element for each site
  std::unordered_map<std::string, int> orbit_map; //Added by JYL 

  std::vector<std::vector<Complex>> ylm;  // spherical harmonics
};

/// @brief Spherical excess (= area on the unit sphere, steradians) of the
/// spherical triangle with unit-vector vertices a, b, c, via L'Huilier's
/// theorem. Used by EqualAreaGnomonicPoint to locate points by exact area
/// fraction rather than by chord-length or azimuth fraction (neither of
/// which is area-preserving on the sphere). See IsingS2_precision/journal.md
/// 2026-08-24 "analytic equal-area mesh" entry for the derivation.
static inline double ClampUnit(double x) {
  return x < -1.0 ? -1.0 : (x > 1.0 ? 1.0 : x);
}

static double SphericalTriangleArea(const Vec3& a, const Vec3& b,
                                     const Vec3& c) {
  double side_a = acos(ClampUnit(b.dot(c)));
  double side_b = acos(ClampUnit(a.dot(c)));
  double side_c = acos(ClampUnit(a.dot(b)));
  double s = 0.5 * (side_a + side_b + side_c);
  double p = std::max(0.0, tan(0.5 * s) * tan(0.5 * (s - side_a)) *
                                tan(0.5 * (s - side_b)) *
                                tan(0.5 * (s - side_c)));
  return 4.0 * atan(sqrt(p));
}

/// @brief Exact equal-area placement of grid point (x,y) (x,y >= 0,
/// x+y <= k) within the spherical triangle face_r[0],face_r[1],face_r[2],
/// generalizing Snyder's (1992) equal-area polyhedral projection to a
/// direct per-point construction (no global mesh relaxation, unlike
/// EqualizeFaceAreas/EqualizeDualAreas below).
///
/// Method: target barycentric weights are beta0=(k-x-y)/k on face_r[0],
/// beta1=x/k on face_r[1], beta2=y/k on face_r[2] (matching the naive
/// flat-subdivide scheme's convention exactly, so face/mesh combinatorics
/// are unchanged). Write h = beta1+beta2 = 1-beta0. The point v0 = face_r[0]
/// is treated as an apex; p is the point on the opposite edge
/// face_r[1]-face_r[2] (a geodesic) with
/// SphericalTriangleArea(v0,face_r[1],p) / SphericalTriangleArea(v0,face_r[1],face_r[2])
/// = beta2/h exactly -- found by bisection along the edge (monotonic area,
/// closed-form area evaluation per step, not a mesh-wide relaxation).
/// The target point v then sits on the geodesic ray v0->p at the exact
/// equal-area radial fraction h, i.e. (1-cos(angle(v0,v))) =
/// h*(1-cos(angle(v0,p))) -- the standard Lambert-azimuthal-equal-area
/// radial law, valid here because p was already chosen so azimuthal
/// "wedges" toward the edge carry exactly proportional area. This
/// composition reproduces exact spherical-area partitioning in both
/// coordinates (edge-position and radial fraction), matching Snyder's
/// method.
///
/// IMPORTANT: this construction treats face_r[0] as a distinguished apex
/// and is NOT symmetric under permuting the three vertices -- the
/// opposite edge (x+y=k) is placed by area-ratio relative to face_r[0],
/// while the two near edges (x=0, y=0) are placed by radial equal-area
/// fraction from face_r[0], a different rule. Two faces sharing an edge
/// generally disagree on which vertex is "face_r[0]" for that edge, so
/// applying this formula AT the boundary (x=0, y=0, or x+y=k) makes
/// adjacent faces compute different points for the same shared site --
/// confirmed empirically 2026-08-24 (mesh construction assertion failure,
/// `s_next <= n_sites`, at n_refine=4). The caller therefore only invokes
/// this for strictly interior points (0<x, 0<y, x+y<k); boundary points
/// keep the naive linear-interpolation-then-normalize placement, which
/// depends only on the two shared edge endpoints and t (not on the third
/// vertex), so it is inherently face-order-independent and safe. This
/// means only interior-to-face triangles get the equal-area correction;
/// the boundary layer keeps naive-mesh area distortion, a smaller and
/// smaller fraction of the mesh (O(1/n_refine)) as resolution increases.
/// @param face_r The face's 3 vertices (unit vectors), same convention as
/// the flat-subdivide loop that calls this (face_r[0]+x*n_x+y*n_y).
/// @param x,y Grid indices, 0 < x, 0 < y, x+y < k (strictly interior only
/// -- see note above).
/// @param k Refinement level of this face.
Vec3 QfeLatticeS2::EqualAreaGnomonicPoint(Vec3 face_r[3], int x, int y,
                                           int k) {
  const Vec3& v0 = face_r[0];
  const Vec3& v1 = face_r[1];
  const Vec3& v2 = face_r[2];

  double h = double(x + y) / double(k);
  double r_target = (double(y) / double(k)) / h;  // beta2 / h, in [0,1]

  double total_area = SphericalTriangleArea(v0, v1, v2);
  double target_area = r_target * total_area;

  // bisect for t such that p(t) = normalize((1-t)*v1 + t*v2) satisfies
  // SphericalTriangleArea(v0, v1, p(t)) == target_area. Area(v0,v1,p(t))
  // is continuous and monotonically increasing from 0 (t=0, p=v1) to
  // total_area (t=1, p=v2), so bisection converges unconditionally.
  double t_lo = 0.0, t_hi = 1.0;
  Vec3 p = v2;
  for (int iter = 0; iter < 60; iter++) {
    double t_mid = 0.5 * (t_lo + t_hi);
    p = ((1.0 - t_mid) * v1 + t_mid * v2).normalized();
    double area_mid = SphericalTriangleArea(v0, v1, p);
    if (area_mid < target_area) {
      t_lo = t_mid;
    } else {
      t_hi = t_mid;
    }
  }

  double cos_theta_p = ClampUnit(v0.dot(p));
  double theta_p = acos(cos_theta_p);
  if (theta_p < 1e-14) return v0;  // degenerate (k too small to matter)

  // NOTE: uses h*h, not h -- h=(x+y)/k behaves like a radius-type
  // coordinate on the flat (x,y) grid (linear grid spacing), not a literal
  // area fraction. A uniform (x,y) grid needs area to grow like h^2 along
  // a fixed-azimuth ray for the map to be equal-area (exactly analogous to
  // ordinary polar coordinates needing an "r dr" measure, not "dr", to
  // preserve area for a uniformly-spaced radial index) -- see journal.md
  // 2026-08-24 for the Jacobian derivation. Using plain h here was the
  // bug behind the first (reverted) version of this function, which made
  // area non-uniformity get *worse*, not better, with mesh refinement.
  double cos_theta_v = ClampUnit(1.0 - h * h * (1.0 - cos_theta_p));
  double theta_v = acos(cos_theta_v);

  Vec3 e_dir = (p - cos_theta_p * v0).normalized();
  return (cos(theta_v) * v0 + sin(theta_v) * e_dir).normalized();
}

/// @brief Create a simplicial discretization of a 2-sphere.
/// @param q Number of links meeting at each site. Valid values for @p q are 3,
/// 4, and 5 for a tetrahedron, octahedron, and icosahedron, respectively.
/// @param k Refinement level.
QfeLatticeS2::QfeLatticeS2(int q, int k, bool equal_area_analytic) {
  // refinement level must be positive
  assert(k >= 1);

  // tetrahedron, octahedron, and icosahedran are the only valid base
  // polyhedrons
  assert(q >= 3 && q <= 5);

  this->q = q;

  // return base lattice if unrefined
  if (k == 1) {
    ReadBaseLattice(q);
    UpdateDistinct();
    ReadSymmetryData(q, k);
    return;
  }

  if (q == 3) {
    // k-refined tetrahedron has V = 2 k^2 + 2 vertices
    ResizeSites(2 * k * k + 2);
  } else if (q == 4) {
    // k-refined octahedron has V = 4 k^2 + 2 vertices
    ResizeSites(4 * k * k + 2);
  } else if (q == 5) {
    // k-refined icosahedron has V = 10 k^2 + 2 vertices
    ResizeSites(10 * k * k + 2);
  }

  // set the lattice volume
  vol = double(n_sites);

  // create an unrefined base lattice
  QfeLatticeS2 base_lattice(q);

  // maps to identify vertices, orbits, and faces
  std::unordered_map<std::string, int> coord_map;
  //std::unordered_map<std::string, int> orbit_map;
  std::unordered_map<std::string, int> face_map;

  // index of next site and face to create
  int s_next = 0;
  int f_next = 0;

  // xy coordinates of each site
  std::vector<int> site_x(n_sites);
  std::vector<int> site_y(n_sites);

  // loop over faces of base polyhedron
  for (int f = 0; f < base_lattice.n_faces; f++) {
    // get the vertices of the base polyhedron face
    Vec3 face_r[3];
    for (int i = 0; i < 3; i++) {
      int b = base_lattice.faces[f].sites[i];
      face_r[i] = base_lattice.r[b];
      if (f == 0) {
        // set the coordinates of the first face's vertices
        first_face_r[i] = face_r[i];
      }
    }

    // unit vectors in the face in xy basis
    Vec3 n_x = (face_r[1] - face_r[0]) / double(k);
    Vec3 n_y = (face_r[2] - face_r[0]) / double(k);

    // list of sites in this face labeled by xy positions on the face
    int k1 = k + 1;
    int xy_max = k1 * k1;
    int xy_list[k1][k1];

    // loop over xy to find all sites and set their positions
    for (int xy = 0; xy <= xy_max; xy++) {
      int x = xy % k1;
      int y = xy / k1;

      // skip xy values outside the triangle
      if ((x + y) > k) continue;

      // calculate the coordinates of this vertex. Always the naive
      // placement here, regardless of equal_area_analytic -- see
      // QfeLatticeS2::ApplyEqualAreaAnalytic (called after
      // ReadSymmetryData, once site_g/orbit bookkeeping exists) for why:
      // computing the analytic placement independently per-face broke
      // icosahedral symmetry (EqualAreaGnomonicPoint is not symmetric
      // under permuting which vertex is the apex), so instead a single
      // analytic representative is computed per orbit on face 0 and
      // propagated to every other site as an exact group-element image,
      // which is symmetric by construction regardless of the per-point
      // formula's own symmetry. See journal.md 2026-08-24 "icosahedral
      // symmetry check" entries for the full story.
      Vec3 v = (face_r[0] + x * n_x + y * n_y).normalized();

      // deal with negative zero
      std::string vec_name = Vec3ToString(v);

      // check if the site already exists
      if (coord_map.find(vec_name) == coord_map.end()) {
        // sorted barycentric coordinates define the site orbit
        int xi[3];
        xi[0] = x;
        xi[1] = y;
        xi[2] = k - x - y;
        std::sort(xi, xi + 3, std::greater<int>());
        std::string orbit_name = string_format("%d_%d_%d", xi[0], xi[1], xi[2]);

        // check if the orbit already exists
        if (orbit_map.find(orbit_name) == orbit_map.end()) {
          // create a new orbit
          double xi1 = double(xi[0]) / double(k);
          double xi2 = double(xi[1]) / double(k);
          orbit_map[orbit_name] = CreateOrbit(xi1, xi2);
        }

        // create a new site
        int orbit_id = orbit_map[orbit_name];
        coord_map[vec_name] = s_next;
        r[s_next] = v;
        sites[s_next].nn = 0;
        sites[s_next].wt = 1.0;
        sites[s_next].id = orbit_id;
        site_orbit[s_next] = orbit_id;
        s_next++;
        assert(s_next <= n_sites);
      }

      // get the site index
      int s = coord_map[vec_name];

      // save this site in the xy list
      xy_list[x][y] = s;
      site_x[s] = x;
      site_y[s] = y;
    }

    // add faces and links
    for (int xy = 0; xy <= xy_max; xy++) {
      int x = xy % k1;
      int y = xy / k1;
      if ((x + y) > k) continue;

      if ((x + y) != k) {
        // "forward" triangle
        int s1 = xy_list[x][y];
        int s2 = xy_list[x][y + 1];
        int s3 = xy_list[x + 1][y];
        AddFace(s1, s2, s3);
      }

      if ((x != 0) && (y != 0)) {
        // "backward" triangle
        int s1 = xy_list[x][y];
        int s2 = xy_list[x][y - 1];
        int s3 = xy_list[x - 1][y];
        AddFace(s1, s2, s3);
      }
    }

    // set the orbit id for each face
    int n_distinct_faces = 0;
    face_orbit.resize(n_faces);
    while (f_next != n_faces) {
      int xi[3] = {0, 0, 0};

      // compute the barycentric coordinates of this face
      for (int i = 0; i < 3; i++) {
        int s = faces[f_next].sites[i];
        int x = site_x[s];
        int y = site_y[s];

        xi[0] += x;
        xi[1] += y;
        xi[2] += k - x - y;
      }

      // sorted barycentric coordinates define the face orbit
      std::sort(xi, xi + 3, std::greater<int>());
      std::string face_name = string_format("%d_%d_%d", xi[0], xi[1], xi[2]);

      // check if the face already exists
      if (face_map.find(face_name) == face_map.end()) {
        face_map[face_name] = n_distinct_faces++;
      }
      face_orbit[f_next] = face_map[face_name];
      f_next++;
    }
  }

  // check that all of the sites and faces have been created
  assert(s_next == n_sites);
  assert(f_next == n_faces);

  // read the symmetry group data (always against the naive positions just
  // built above -- see ApplyEqualAreaAnalytic for why the analytic
  // placement is applied as a separate pass afterward rather than during
  // construction).
  UpdateDistinct();
  ReadSymmetryData(q, k);

  if (equal_area_analytic) {
    ApplyEqualAreaAnalytic(k);
  }
}

/// @brief Read base polyhedron data from grp directory
/// @param q polyhedron parameter
void QfeLatticeS2::ReadBaseLattice(int q) {
  // tetrahedron, octahedron, and icosahedron are the only valid base
  // polyhedrons
  assert(q >= 3 && q <= 5);
  this->q = q;

  // read the base lattice file in the symmetry group directory
  std::string lattice_path = string_format("%s/lattice/o3q%d.dat", GRP_DIR, q);
  FILE* file = fopen(lattice_path.c_str(), "r");
  assert(file != nullptr);
  ReadLattice(file);
  fclose(file);

  // set the lattice volume
  vol = double(n_sites);

  // create a single orbit
  CreateOrbit(0.0, 0.0);

  // initialize the first face vertex coordinates
  for (int i = 0; i < 3; i++) {
    int s = faces[0].sites[i];
    first_face_r[i] = r[s];
  }

  face_orbit.resize(n_faces);
  for (int f = 0; f < n_faces; f++) face_orbit[f] = 0;
}

/// @brief Write a site to a lattice file
/// @param file Lattice file
/// @param s Site index
void QfeLatticeS2::WriteSite(FILE* file, int s) {
  QfeLattice::WriteSite(file, s);
  double theta = acos(r[s].z());
  double phi = atan2(r[s].y(), r[s].x());
  fprintf(file, " %+.20f %+.20f", theta, phi);
}

/// @brief Read a site from a lattice file
/// @param file Lattice file
/// @param s Site index
void QfeLatticeS2::ReadSite(FILE* file, int s) {
  QfeLattice::ReadSite(file, s);
  double theta, phi;
  fscanf(file, " %lf %lf", &theta, &phi);

  r[s][0] = sin(theta) * cos(phi);
  r[s][1] = sin(theta) * sin(phi);
  r[s][2] = cos(theta);
  r[s].normalize();

  Vec3 north_pole(0.0, 0.0, 1.0);
  Vec3 south_pole(0.0, 0.0, -1.0);
  if (AlmostEq(r[s], north_pole)) r[s] = north_pole;
  if (AlmostEq(r[s], south_pole)) r[s] = south_pole;
}

/// @brief Create an orbit
/// @param xi1 1st barycentric coordinate
/// @param xi2 2nd barycentric coordinate
int QfeLatticeS2::CreateOrbit(double xi1, double xi2) {
  // barycentric coordinates, sorted to account for degeneracies
  double xi[3] = {xi1, xi2, 1.0 - xi1 - xi2};
  std::sort(xi, xi + 3, std::greater<double>());
  int o = orbit_xi.size();
  orbit_xi.push_back(Vec3(xi));
  return o;
}

/// @brief Calculate the coordinates of the first site in an orbit
/// @param o Orbit index
/// @return Normalized orbit coordinates
Vec3 QfeLatticeS2::CalcOrbitPos(int o) {
  Vec3 r = Vec3::Zero();
  Vec3 xi = orbit_xi[o];
  for (int i = 0; i < 3; i++) {
    r += xi(i) * first_face_r[i];
  }
  return r.normalized();
}

/// @brief Read an orbit file and convert to site coordinates
/// @param file Orbit file
void QfeLatticeS2::ReadOrbits(FILE* file) {
  // make sure we're at the beginning of the file
  fseek(file, 0L, SEEK_SET);

  // read orbit data
  orbit_xi.resize(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    // read barycentric coordinates
    int o_check;
    fscanf(file, "%d", &o_check);
    assert(o_check == o);

    // read dof values
    double xi_sum = 0.0;
    for (int i = 0; i < 2; i++) {
      double temp;
      fscanf(file, "%lf", &temp);
      xi_sum += temp;
      orbit_xi[o][i] = temp;
    }
    orbit_xi[o][2] = 1.0 - xi_sum;
    fscanf(file, "\n");
  }
  assert(feof(file));

  UpdateOrbits();
}

/// @brief Write orbit barycentric coordinates to a file that can be read
/// via ReadOrbits
/// @param file Orbit file
void QfeLatticeS2::WriteOrbits(FILE* file) {
  // read orbit data
  for (int o = 0; o < n_distinct; o++) {
    // read barycentric coordinates
    fprintf(file, "%d", o);

    // read dof values
    for (int i = 0; i < 2; i++) {
      fprintf(file, " %.16f", orbit_xi[o][i]);
    }
    fprintf(file, "\n");
  }
}

/// @brief Re-place every site using EqualAreaGnomonicPoint for interior
/// grid points, while guaranteeing exact icosahedral symmetry by
/// construction: only ONE representative site per orbit (on face 0) is
/// ever evaluated by the (apex-asymmetric, not itself point-group-
/// symmetric) analytic formula; every other site in that orbit is set to
/// an exact group-element image of that one representative
/// (`G[site_g[s]] * orbit_r[site_orbit[s]]`, the same pattern
/// UpdateOrbits() uses). This is why the constructor calls this as a
/// separate pass *after* ReadSymmetryData -- site_g/site_orbit must
/// already be valid, and they're only computed correctly against the
/// naive (exactly-symmetric) placement (see the constructor and
/// ReadSymmetryData). Must be called with @p k equal to the constructor's
/// refine level (not stored as a member).
/// @param k Refinement level of this lattice (same as constructor's k).
void QfeLatticeS2::ApplyEqualAreaAnalytic(int k) {
  // Orbit multiplicity (how many sites actually carry a given orbit index)
  // reveals orbits with a nontrivial icosahedral stabilizer -- points
  // sitting exactly on a symmetry axis/mirror plane, generic multiplicity
  // is G.size(), special orbits have a proper divisor of it. Propagating
  // via "one analytic representative + existing site_g coset
  // representatives" is only valid if the representative respects that
  // same stabilizer (else the propagated set doesn't close under G --
  // confirmed empirically 2026-08-24, an interior orbit with
  // multiplicity 60 < n_group=120 gave a real, non-round-off mismatch).
  // The naive formula is always safe (it's what site_g/the stabilizer
  // structure was computed against in the first place), so special orbits
  // fall back to it, same as boundary orbits.
  std::vector<int> orbit_mult(n_distinct, 0);
  for (int s = 0; s < n_sites; s++) orbit_mult[site_orbit[s]]++;

  std::vector<Vec3> orbit_r(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    // orbit_xi[o] holds SORTED-DESCENDING barycentric fractions, and
    // CalcOrbitPos assigns them to face_r[0],[1],[2] in that order (largest
    // weight on face_r[0]). The naive grid formula's weight on face_r[0] is
    // (k-x-y)/k, on face_r[1] is x/k, on face_r[2] is y/k -- so matching
    // CalcOrbitPos's convention requires x from orbit_xi[o](1) (the
    // *second* entry) and y from orbit_xi[o](2) (the third), NOT (0) and
    // (1). Got this backwards on the first attempt (confirmed via a
    // max_mismatch regression from 1.5e-2 to 0.30 -- see journal.md
    // 2026-08-24 "orbit weight assignment" entry).
    int x = int(std::lround(orbit_xi[o](1) * k));
    int y = int(std::lround(orbit_xi[o](2) * k));
    bool interior = x > 0 && y > 0 && (x + y) < k;
    bool generic = orbit_mult[o] == int(G.size());
    orbit_r[o] = (interior && generic) ? EqualAreaGnomonicPoint(first_face_r, x, y, k)
                                        : CalcOrbitPos(o);
  }

  Vec3 north_pole(0.0, 0.0, 1.0);
  Vec3 south_pole(0.0, 0.0, -1.0);
  for (int s = 0; s < n_sites; s++) {
    int o = site_orbit[s];
    int g = site_g[s];
    r[s] = G[g] * orbit_r[o];
    r[s].normalize();
    if (AlmostEq(r[s], north_pole)) r[s] = north_pole;
    if (AlmostEq(r[s], south_pole)) r[s] = south_pole;
  }
}

/// @brief Update positions of all sites using orbits and symmetry group data
void QfeLatticeS2::UpdateOrbits() {
  // calculate the orbit positions
  std::vector<Vec3> orbit_r(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    orbit_r[o] = CalcOrbitPos(o);
  }

  // use the symmetry group data to calculate site coordinates
  Vec3 north_pole(0.0, 0.0, 1.0);
  Vec3 south_pole(0.0, 0.0, -1.0);
  for (int s = 0; s < n_sites; s++) {
    int o = site_orbit[s];
    int g = site_g[s];
    r[s] = G[g] * orbit_r[o];
    r[s].normalize();
    if (AlmostEq(r[s], north_pole)) r[s] = north_pole;
    if (AlmostEq(r[s], south_pole)) r[s] = south_pole;
  }
}

/// @brief Read symmetry group data from pre-generated files in the grp
/// directory
/// @param q polyhedron parameter
/// @param k refinement level
void QfeLatticeS2::ReadSymmetryData(int q, int k) {
  // open the symmetry group data file
  std::string grp_path = string_format("%s/elem/o3q%d.dat", GRP_DIR, q);
  FILE* grp_file = fopen(grp_path.c_str(), "r");
  assert(grp_file != nullptr);

  // read group elements
  G.clear();
  while (!feof(grp_file)) {
    GrpElemO3 g;
    g.ReadGrpElem(grp_file);
    G.push_back(g);
  }
  fclose(grp_file);

  // calculate all of the orbit positions
  std::vector<Vec3> orbit_r(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    orbit_r[o] = CalcOrbitPos(o);
  }

  // open the site group element file
  std::string g_path = string_format("%s/site_g/o3q%dk%d.dat", GRP_DIR, q, k);
  FILE* g_file = fopen(g_path.c_str(), "r");
  bool site_g_success = true;
  if (g_file != nullptr) {
    // load pre-existing symmetry data
    site_g.resize(n_sites);
    std::vector<int>::iterator it = site_g.begin();
    while (!feof(g_file)) {
      assert(it != site_g.end());
      int g;
      fscanf(g_file, "%d\n", &g);
      *it++ = g;
    }
    fclose(g_file);

    // recalculate if the file was not long enough
    if (it != site_g.end()) {
      fprintf(stderr, "Rebuilding invalid data file: %s\n", g_path.c_str());
      site_g_success = false;
    } else {
      for (int s = 0; s < n_sites; s++) {
        int o = site_orbit[s];
        int g = site_g[s];
        Vec3 r_norm = r[s].normalized();
        Vec3 gr = G[g] * orbit_r[o];
        if (!AlmostEq(r_norm, gr, 1.0e-15)) {
          site_g_success = false;
          break;
        }
      }

      if (!site_g_success) {
        fprintf(stderr, "Rebuilding invalid data file: %s\n", g_path.c_str());
      }
    }

  } else {
    site_g_success = false;
  }

  if (!site_g_success) {
    // find the group element for each site
    site_g.resize(n_sites);
    for (int s = 0; s < n_sites; s++) {
      int o = site_orbit[s];
      Vec3 r_norm = r[s].normalized();

      // find the appropriate group element
      bool found_g = false;
      for (int g = 0; g < G.size(); g++) {
        Vec3 gr = G[g] * orbit_r[o];
        if (!AlmostEq(r_norm, gr, 1.0e-15)) continue;
        site_g[s] = g;
        found_g = true;
        break;
      }
      assert(found_g);
    }

    // write the site group elements to a file
    g_file = fopen(g_path.c_str(), "w");
    for (int s = 0; s < n_sites; s++) {
      fprintf(g_file, "%d\n", site_g[s]);
    }
    fclose(g_file);
  }

  UpdateOrbits();
}

/// @brief Change the number of sites.
/// @param n_sites New number of sites
void QfeLatticeS2::ResizeSites(int n_sites) {
  QfeLattice::ResizeSites(n_sites);
  r.resize(n_sites);
  ylm.resize(n_sites);
  antipode.resize(n_sites, 0);
  site_orbit.resize(n_sites);
}

/// @brief Project all site coordinates onto a unit sphere.
void QfeLatticeS2::Inflate() {
  for (int s = 0; s < n_sites; s++) {
    r[s].normalize();
  }
}

/// @brief Identify each site's antipode, i.e. for a site with position r, find
/// the site which has position -r. A lattice with a tetrahedron base (q = 3)
/// does not have an antipode for every site.
void QfeLatticeS2::UpdateAntipodes() {
  std::unordered_map<std::string, int> antipode_map;
  for (int s = 0; s < n_sites; s++) {
    // find antipode
    Vec3 anti_r = -r[s];
    std::string key = Vec3ToString(r[s], 6);
    std::string anti_key = Vec3ToString(anti_r, 6);

    if (antipode_map.find(anti_key) != antipode_map.end()) {
      // antipode found in map
      int a = antipode_map[anti_key];
      antipode[s] = a;
      antipode[a] = s;
      antipode_map.erase(anti_key);
    } else {
      // antipode not found yet
      antipode_map[key] = s;
    }
  }

  if (antipode_map.size()) {
    // print error message if there are any unpaired sites
    fprintf(stderr, "no antipode found for %lu/%d sites\n", antipode_map.size(),
            n_sites);
    std::unordered_map<std::string, int>::iterator it;
    for (it = antipode_map.begin(); it != antipode_map.end(); it++) {
      fprintf(stderr, "%04d %s\n", it->second, it->first.c_str());
    }
  }
}

/// @brief Find the circumcenter of face
/// @param f Face id
/// @return Coordinates of face circumcenter
Vec3 QfeLatticeS2::FaceCircumcenter(int f) {
  double sq_edge_1 = EdgeSquared(faces[f].edges[0]);
  double sq_edge_2 = EdgeSquared(faces[f].edges[1]);
  double sq_edge_3 = EdgeSquared(faces[f].edges[2]);

  double w1 = sq_edge_1 * (sq_edge_2 + sq_edge_3 - sq_edge_1);
  double w2 = sq_edge_2 * (sq_edge_3 + sq_edge_1 - sq_edge_2);
  double w3 = sq_edge_3 * (sq_edge_1 + sq_edge_2 - sq_edge_3);

  Vec3 r1 = w1 * r[faces[f].sites[0]];
  Vec3 r2 = w2 * r[faces[f].sites[1]];
  Vec3 r3 = w3 * r[faces[f].sites[2]];

  return (r1 + r2 + r3) / (w1 + w2 + w3);
}

/// @brief Calculate the squared length of link
/// @param l Link id
/// @return Squared length of link
double QfeLatticeS2::EdgeSquared(int l) {
  int s_a = links[l].sites[0];
  int s_b = links[l].sites[1];
  Vec3 dr = r[s_a] - r[s_b];
  return dr.squaredNorm();
}

/// @brief Calculate the length of a link
/// @param l Link id
/// @return Link length
double QfeLatticeS2::EdgeLength(int l) { return sqrt(EdgeSquared(l)); }

/// @brief Calculate the flat area of a triangular face.
/// @param f Face id
/// @return Area of triangular face
double QfeLatticeS2::FlatArea(int f) {
  double a = EdgeLength(faces[f].edges[0]);
  double b = EdgeLength(faces[f].edges[1]);
  double c = EdgeLength(faces[f].edges[2]);
  double area = (a + b + c) * (b + c - a) * (c + a - b) * (a + b - c);
  return 0.25 * sqrt(area);
}

double QfeLatticeS2::Perimeter(int f) {
  double a = EdgeLength(faces[f].edges[0]);
  double b = EdgeLength(faces[f].edges[1]);
  double c = EdgeLength(faces[f].edges[2]);
  return a+b+c;
}

/// @brief Circumradius of a flat triangular face, R = a*b*c/(4*Area).
/// Standard formula for the circumradius of a Euclidean triangle with
/// side lengths a,b,c and area Area -- exact for `faces[f]` treated as a
/// flat triangle in the embedding space (the same "flat face" convention
/// `FlatArea`/`Perimeter` already use).
double QfeLatticeS2::Circumradius(int f) {
  double a = EdgeLength(faces[f].edges[0]);
  double b = EdgeLength(faces[f].edges[1]);
  double c = EdgeLength(faces[f].edges[2]);
  return (a * b * c) / (4.0 * FlatArea(f));
}

/// @brief Calculate FEM weights based on vertex coordinates.
void QfeLatticeS2::UpdateWeights() {
  // set site weights to zero
  for (int s = 0; s < n_sites; s++) {
    sites[s].wt = 0.0;
  }

  // loop over links to update weights
  for (int l = 0; l < n_links; l++) {
    links[l].wt = 0.0;
    for (int i = 0; i < 2; i++) {
      // find the other two edges of this face
      int f = links[l].faces[i];
      int e = 0;
      while (faces[f].edges[e] != l) e++;
      int e1 = (e + 1) % 3;
      int e2 = (e + 2) % 3;
      int l1 = faces[f].edges[e1];
      int l2 = faces[f].edges[e2];

      // find the area associated with this face
      double sq_edge = EdgeSquared(l);
      double sq_edge_1 = EdgeSquared(l1);
      double sq_edge_2 = EdgeSquared(l2);
      double half_wt = (sq_edge_1 + sq_edge_2 - sq_edge) / (8.0 * FlatArea(f));
      links[l].wt += half_wt;

      // add to the weights of the two sites connected by this link
      sites[links[l].sites[0]].wt += 0.25 * half_wt * sq_edge;
      sites[links[l].sites[1]].wt += 0.25 * half_wt * sq_edge;
    }
  }

  // normalize site weights to 1
  double site_wt_sum = 0.0;
  for (int s = 0; s < n_sites; s++) {
    site_wt_sum += sites[s].wt;
  }

  double site_wt_norm = site_wt_sum / double(n_sites);
  for (int s = 0; s < n_sites; s++) {
    sites[s].wt /= site_wt_norm;
  }

  // set face weights equal to their flat area
  double face_area_sum = 0.0;
  for (int f = 0; f < n_faces; f++) {
    double face_area = FlatArea(f);
    face_area_sum += face_area;
    faces[f].wt = face_area;
  }

  // normalize face areas
  double face_wt_norm = face_area_sum / double(n_faces);
  for (int f = 0; f < n_faces; f++) {
    faces[f].wt /= face_wt_norm;
  }
}

/// @brief Iteratively relax vertex positions on the sphere to reduce the
/// variance of triangular face areas, starting from the naive
/// flat-subdivide-then-radially-project mesh (constructor default). Each
/// site is nudged toward the centroids of adjacent faces that are larger
/// than the mean face area and away from centroids of faces smaller than
/// the mean, then re-projected onto the unit sphere. Because the update
/// rule is a fixed, purely local (combinatorial-neighbor-based) function
/// applied identically to every site every iteration, and the starting
/// mesh has exact icosahedral symmetry, the relaxed mesh remains
/// icosahedrally symmetric to floating-point precision -- no orbit/group
/// bookkeeping needs to change. See IsingS2_precision/journal.md
/// 2026-08-24 "equal-area mesh" entry for the motivation (the naive mesh's
/// non-uniform face areas are a candidate explanation for the observed
/// non-vanishing l>=3 spherical-symmetry-breaking plateau).
/// @param n_iter Number of relaxation iterations
/// @param step Step size (fraction of the area-weighted displacement
/// applied per iteration; too large can overshoot/tangle the mesh)
/// @brief Sum of face areas and the objective E = sum_f (A_f - mean)^2,
/// used by the gradient-descent version of EqualizeFaceAreas below.
static void FaceAreaEnergy(QfeLatticeS2& lattice, std::vector<double>& area,
                            double& mean, double& E) {
  area.resize(lattice.n_faces);
  double sum = 0.0;
  for (int f = 0; f < lattice.n_faces; f++) {
    area[f] = lattice.FlatArea(f);
    sum += area[f];
  }
  mean = sum / double(lattice.n_faces);
  E = 0.0;
  for (int f = 0; f < lattice.n_faces; f++) {
    double d = area[f] - mean;
    E += d * d;
  }
}

/// @brief Gradient descent (with backtracking line search) on
/// E = sum_f (FlatArea(f) - mean)^2, replacing an earlier heuristic
/// version (fixed unit-direction-toward-imbalanced-centroid update, not a
/// true gradient of any objective) that was found 2026-08-24 to have no
/// stable equilibrium at n_refine>=32 regardless of step size -- smaller
/// steps only delayed eventual drift away from the best point reached,
/// confirming the heuristic itself (not the step size) was the problem.
/// See journal.md 2026-08-24 "gradient descent" entry for the derivation
/// and comparison.
///
/// Per-vertex gradient of a single triangle's flat area A(a,b,c) =
/// 0.5*|(b-a)x(c-a)|: dA/da = 0.5*(c-b) x n_hat (n_hat = unit face
/// normal), and cyclic permutations for b, c -- standard closed-form
/// triangle-area gradient. The objective's dependence on `mean` through
/// every other face is dropped when building the per-iteration gradient
/// (a standard, cheap approximation -- each single face's sensitivity to
/// `mean` is O(1/n_faces)) but `mean` is recomputed exactly (via
/// FaceAreaEnergy) whenever E is evaluated to accept/reject a step, so the
/// backtracking decision is always exact.
///
/// After computing the gradient, it is projected onto the tangent plane
/// at each site (removing the radial component) before stepping and
/// renormalizing onto the unit sphere -- projected/Riemannian gradient
/// descent on the sphere.
/// @param n_iter Maximum number of accepted gradient steps.
/// @param step Initial step size; backtracking halves it on a rejected
/// trial step and grows it 1.2x after an accepted one, so (unlike the
/// previous version) this does not need per-resolution hand-tuning.
///
/// Early-stop criterion (added for the mesh-position cache, see
/// IsingS2_precision journal.md 2026-08-25 "iters-needed model" entry):
/// stops once E's relative decrease has been below `rel_tol` for
/// `patience` consecutive iterations, in addition to the pre-existing
/// hard stop (no improving step found at all within 30 backtracks). The
/// hard stop alone was found not to trigger within thousands of
/// iterations at n_refine>=32 -- the objective keeps making minuscule
/// improvements indefinitely (the backtracking step-size growth factor
/// 1.2x means a tiny accepted step is still "improved") long after the
/// area-uniformity gain is practically exhausted (e.g. n_refine=32's
/// std/mean plateaus around 0.0058 by iteration ~1000 but the hard stop
/// still hadn't fired by iteration 5000). `rel_tol`/`patience` give a
/// practical stopping point instead of chasing floating-point-scale
/// gradient descent.
void QfeLatticeS2::EqualizeFaceAreas(int n_iter, double step,
                                      int* iters_used, double rel_tol,
                                      int patience) {
  std::vector<double> area;
  double mean, E;
  FaceAreaEnergy(*this, area, mean, E);

  std::vector<Vec3> saved_r(n_sites);
  std::vector<Vec3> grad(n_sites);
  int iter;
  int slow_streak = 0;
  double prev_rel_decrease = std::numeric_limits<double>::infinity();
  for (iter = 0; iter < n_iter; iter++) {
    double E_before = E;
    std::fill(grad.begin(), grad.end(), Vec3::Zero());
    for (int f = 0; f < n_faces; f++) {
      int s0 = faces[f].sites[0], s1 = faces[f].sites[1], s2 = faces[f].sites[2];
      const Vec3& v0 = r[s0];
      const Vec3& v1 = r[s1];
      const Vec3& v2 = r[s2];
      Vec3 nvec = (v1 - v0).cross(v2 - v0);
      double nnorm = nvec.norm();
      if (nnorm < 1e-14) continue;
      Vec3 n_hat = nvec / nnorm;
      // dE/dA_f = 2*(A_f - mean); dA_f/dv_i = 0.5*(v_{i-1} - v_{i+1}) x
      // n_hat (cyclic; verified numerically -- the naively-plausible
      // opposite-edge order (v_{i+1}-v_{i-1}) is the *negative* of the
      // true gradient, confirmed by perturbing a flat test triangle and
      // checking area actually increased along the computed direction;
      // this was the bug behind the first version finding no improving
      // step in either sign, backtracking to full revert every iteration
      // -- see journal.md 2026-08-24 "gradient descent" entry).
      double coeff = (area[f] - mean);  // the factor of 2 is absorbed into
                                         // `step` (an overall constant).
      grad[s0] += coeff * (v1 - v2).cross(n_hat);
      grad[s1] += coeff * (v2 - v0).cross(n_hat);
      grad[s2] += coeff * (v0 - v1).cross(n_hat);
    }
    for (int s = 0; s < n_sites; s++) {
      grad[s] -= grad[s].dot(r[s]) * r[s];  // project onto tangent plane
    }

    saved_r = r;
    double trial_step = step;
    bool improved = false;
    for (int bt = 0; bt < 30; bt++) {
      for (int s = 0; s < n_sites; s++) {
        r[s] = (saved_r[s] - trial_step * grad[s]).normalized();
      }
      std::vector<double> area2;
      double mean2, E2;
      FaceAreaEnergy(*this, area2, mean2, E2);
      if (E2 < E) {
        area = area2;
        mean = mean2;
        E = E2;
        step = trial_step * 1.2;
        improved = true;
        break;
      }
      trial_step *= 0.5;
    }
    if (!improved) {
      r = saved_r;  // no improving step found within 30 backtracks -- converged
      break;
    }
    double rel_decrease = (E_before - E) / std::max(E_before, 1e-300);
    if (getenv("GD_DEBUG") != nullptr) {
      fprintf(stderr, "gd_debug iter=%d E=%.10e rel_decrease=%.6e step=%.6e\n",
              iter, E, rel_decrease, step);
    }
    // Only count an iteration toward the "stalled" streak if progress is
    // both below rel_tol *and* not still accelerating (rel_decrease still
    // growing iteration-over-iteration means the backtracking step size is
    // still in its 1.2x/iteration warm-up ramp, not stalled -- see
    // journal.md 2026-08-25 "iters-needed model, early-stop bug" entry:
    // without the rel_decrease<=prev_rel_decrease guard this fired after
    // exactly `patience` iterations at n_refine=128 even though every one
    // of those iterations was still increasing its per-step gain).
    if (rel_decrease < rel_tol && rel_decrease <= prev_rel_decrease) {
      slow_streak++;
      if (slow_streak >= patience) {
        iter++;  // count this iteration in the reported total
        break;
      }
    } else {
      slow_streak = 0;
    }
    prev_rel_decrease = rel_decrease;
  }
  if (iter > n_iter) iter = n_iter;
  if (iters_used != nullptr) *iters_used = iter;
}

/// @brief Per-face circumradius R, perimeter P, and area A (all standard
/// flat-triangle formulas -- FlatArea/Perimeter/Circumradius above), plus
/// the joint non-uniformity objective E = E_R + E_P + E_A used by
/// EqualizeCircumPerim, where E_X = <X^2>/<X>^2 - 1 (Brower/Owen
/// arXiv:2407.00459 Eq 39's E_R/E_P, extended with the same-form E_A).
/// arXiv:2407.00459 Sec 3.1/3.3/Appendix B derives that this campaign's
/// `duality`/`owen_dual` coupling rules are only exactly critical when
/// every triangle has equal circumradius AND equal perimeter -- their
/// own equal-*area*-only mesh smoother (same style as this file's
/// pre-existing EqualizeFaceAreas) is reported to fail to restore
/// spherical symmetry for exactly this reason. E_A is kept in the joint
/// objective anyway (per user direction 2026-08-25): area equalization
/// was independently shown (IsingS2_precision push7, 2026-08-25) to speed
/// up SO(3)-symmetry restoration, so there is no reason to abandon it --
/// the paper's finding is that area uniformity *alone* is insufficient,
/// not that it hurts when jointly minimized alongside R and P. See
/// IsingS2_precision/journal.md 2026-08-25 "ROOT CAUSE FOUND" entry.
static void CircumPerimAreaEnergy(QfeLatticeS2& lattice, std::vector<double>& R,
                                   std::vector<double>& P, std::vector<double>& A,
                                   double& E) {
  int F = lattice.n_faces;
  R.resize(F);
  P.resize(F);
  A.resize(F);
  double sumR = 0.0, sumR2 = 0.0, sumP = 0.0, sumP2 = 0.0, sumA = 0.0, sumA2 = 0.0;
  for (int f = 0; f < F; f++) {
    R[f] = lattice.Circumradius(f);
    P[f] = lattice.Perimeter(f);
    A[f] = lattice.FlatArea(f);
    sumR += R[f]; sumR2 += R[f] * R[f];
    sumP += P[f]; sumP2 += P[f] * P[f];
    sumA += A[f]; sumA2 += A[f] * A[f];
  }
  double Rmean = sumR / F, R2mean = sumR2 / F;
  double Pmean = sumP / F, P2mean = sumP2 / F;
  double Amean = sumA / F, A2mean = sumA2 / F;
  double E_R = R2mean / (Rmean * Rmean) - 1.0;
  double E_P = P2mean / (Pmean * Pmean) - 1.0;
  double E_A = A2mean / (Amean * Amean) - 1.0;
  E = E_R + E_P + E_A;
}

/// @brief Exact per-vertex tangent-plane gradient of `E` above via local
/// finite differences: perturbing a single vertex only changes the
/// R/P/A of its (few) incident faces, so the resulting exact new global
/// sums (and hence E) can be recomputed in O(degree) rather than O(F) by
/// subtracting the old per-face contributions and adding the new ones.
/// This is an EXACT gradient of the true global objective (unlike
/// EqualizeFaceAreas's analytic gradient, which explicitly drops each
/// face's O(1/F) sensitivity of `mean` to a single vertex as a documented
/// approximation) -- deliberately implemented via finite differences
/// rather than by hand-deriving closed-form R/P gradients, since a
/// numeric gradient cannot carry a silent sign/algebra bug the way a
/// from-memory analytic derivative could (see directive.md's
/// never-trust-an-unverified-formula rule). Central difference, step
/// `fd_h` (absolute, in embedding-space units -- safe across this
/// campaign's resolution range since edge lengths only reach ~1e-2 at
/// n_refine=128, still far above sqrt(machine epsilon)).
static double LocalEnergyWithVertexMoved(
    QfeLatticeS2& lattice, int s, const Vec3& new_pos,
    const std::vector<int>& sfaces, const std::vector<double>& R,
    const std::vector<double>& P, const std::vector<double>& A, double sumR,
    double sumR2, double sumP, double sumP2, double sumA, double sumA2,
    int F) {
  Vec3 saved = lattice.r[s];
  lattice.r[s] = new_pos;
  double newSumR = sumR, newSumR2 = sumR2;
  double newSumP = sumP, newSumP2 = sumP2;
  double newSumA = sumA, newSumA2 = sumA2;
  for (int f : sfaces) {
    newSumR -= R[f]; newSumR2 -= R[f] * R[f];
    newSumP -= P[f]; newSumP2 -= P[f] * P[f];
    newSumA -= A[f]; newSumA2 -= A[f] * A[f];
    double Rf = lattice.Circumradius(f);
    double Pf = lattice.Perimeter(f);
    double Af = lattice.FlatArea(f);
    newSumR += Rf; newSumR2 += Rf * Rf;
    newSumP += Pf; newSumP2 += Pf * Pf;
    newSumA += Af; newSumA2 += Af * Af;
  }
  lattice.r[s] = saved;
  double Rmean = newSumR / F, R2mean = newSumR2 / F;
  double Pmean = newSumP / F, P2mean = newSumP2 / F;
  double Amean = newSumA / F, A2mean = newSumA2 / F;
  double E_R = R2mean / (Rmean * Rmean) - 1.0;
  double E_P = P2mean / (Pmean * Pmean) - 1.0;
  double E_A = A2mean / (Amean * Amean) - 1.0;
  return E_R + E_P + E_A;
}

/// @brief Gradient descent (backtracking line search, same architecture as
/// EqualizeFaceAreas) on E = E_R + E_P + E_A (circumradius, perimeter, and
/// area non-uniformity jointly -- see CircumPerimAreaEnergy). Per-vertex
/// gradient is computed via exact local finite differences
/// (LocalEnergyWithVertexMoved above), projected onto the tangent plane
/// by construction (perturbation directions e1,e2 span the tangent plane
/// at each site). See IsingS2_precision/journal.md 2026-08-25 "ROOT CAUSE
/// FOUND" entry for why this mesh mode exists (arXiv:2407.00459's
/// coupling rules require circumradius+perimeter uniformity, not just
/// area uniformity, to reach the correct continuum limit).
/// @param n_iter Maximum number of accepted gradient steps.
/// @param step Initial step size; backtracking halves it on a rejected
/// trial step and grows it 1.2x after an accepted one (same convention as
/// EqualizeFaceAreas).
void QfeLatticeS2::EqualizeCircumPerim(int n_iter, double step,
                                        int* iters_used, double rel_tol,
                                        int patience) {
  const double fd_h = 1e-6;

  std::vector<std::vector<int>> site_faces(n_sites);
  for (int f = 0; f < n_faces; f++) {
    for (int e = 0; e < 3; e++) {
      site_faces[faces[f].sites[e]].push_back(f);
    }
  }

  std::vector<double> R, P, A;
  double E;
  CircumPerimAreaEnergy(*this, R, P, A, E);

  std::vector<Vec3> saved_r(n_sites);
  std::vector<Vec3> grad(n_sites);
  int iter;
  int slow_streak = 0;
  double prev_rel_decrease = std::numeric_limits<double>::infinity();
  for (iter = 0; iter < n_iter; iter++) {
    double E_before = E;

    // per-vertex running sums, needed by the local finite-difference
    // energy evaluator (avoids recomputing a full O(F) sum per component)
    double sumR = 0.0, sumR2 = 0.0, sumP = 0.0, sumP2 = 0.0, sumA = 0.0, sumA2 = 0.0;
    for (int f = 0; f < n_faces; f++) {
      sumR += R[f]; sumR2 += R[f] * R[f];
      sumP += P[f]; sumP2 += P[f] * P[f];
      sumA += A[f]; sumA2 += A[f] * A[f];
    }

    for (int s = 0; s < n_sites; s++) {
      const Vec3& x = r[s];
      // arbitrary orthonormal tangent basis at x
      Vec3 ref = (fabs(x.x()) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
      Vec3 e1 = (ref - ref.dot(x) * x).normalized();
      Vec3 e2 = x.cross(e1);

      double Ep1 = LocalEnergyWithVertexMoved(
          *this, s, (x + fd_h * e1).normalized(), site_faces[s], R, P, A,
          sumR, sumR2, sumP, sumP2, sumA, sumA2, n_faces);
      double Em1 = LocalEnergyWithVertexMoved(
          *this, s, (x - fd_h * e1).normalized(), site_faces[s], R, P, A,
          sumR, sumR2, sumP, sumP2, sumA, sumA2, n_faces);
      double Ep2 = LocalEnergyWithVertexMoved(
          *this, s, (x + fd_h * e2).normalized(), site_faces[s], R, P, A,
          sumR, sumR2, sumP, sumP2, sumA, sumA2, n_faces);
      double Em2 = LocalEnergyWithVertexMoved(
          *this, s, (x - fd_h * e2).normalized(), site_faces[s], R, P, A,
          sumR, sumR2, sumP, sumP2, sumA, sumA2, n_faces);

      double g1 = (Ep1 - Em1) / (2.0 * fd_h);
      double g2 = (Ep2 - Em2) / (2.0 * fd_h);
      grad[s] = g1 * e1 + g2 * e2;
    }

    saved_r = r;
    double trial_step = step;
    bool improved = false;
    for (int bt = 0; bt < 30; bt++) {
      for (int s = 0; s < n_sites; s++) {
        r[s] = (saved_r[s] - trial_step * grad[s]).normalized();
      }
      std::vector<double> R2v, P2v, A2v;
      double E2;
      CircumPerimAreaEnergy(*this, R2v, P2v, A2v, E2);
      if (E2 < E) {
        R = R2v; P = P2v; A = A2v;
        E = E2;
        step = trial_step * 1.2;
        improved = true;
        break;
      }
      trial_step *= 0.5;
    }
    if (!improved) {
      r = saved_r;
      break;
    }
    double rel_decrease = (E_before - E) / std::max(E_before, 1e-300);
    if (getenv("GD_DEBUG") != nullptr) {
      fprintf(stderr, "gd_rp_debug iter=%d E=%.10e rel_decrease=%.6e step=%.6e\n",
              iter, E, rel_decrease, step);
    }
    if (rel_decrease < rel_tol && rel_decrease <= prev_rel_decrease) {
      slow_streak++;
      if (slow_streak >= patience) {
        iter++;
        break;
      }
    } else {
      slow_streak = 0;
    }
    prev_rel_decrease = rel_decrease;
  }
  if (iter > n_iter) iter = n_iter;
  if (iters_used != nullptr) *iters_used = iter;
}

/// @brief Dump vertex positions only (unit vectors, one per site) to a
/// flat binary file -- the minimal state EqualizeFaceAreas's relaxation
/// changes. Used by the mesh-position cache (IsingS2_precision
/// journal.md/CLAUDE.md "equal-area mesh cache") so a converged relaxation
/// at a given (q, n_refine, step) never has to be recomputed. Format:
/// int32 n_sites header (a cheap sanity check against the caller's own
/// freshly-constructed naive mesh), then n_sites raw little-endian
/// double[3] entries in site order. Not a general lattice
/// serialization -- topology/links/faces are reconstructed by the normal
/// constructor and assumed unchanged; only positions are cached.
void QfeLatticeS2::WritePositions(FILE* file) {
  int32_t n = n_sites;
  fwrite(&n, sizeof(n), 1, file);
  fwrite(r.data(), sizeof(Vec3), n_sites, file);
}

/// @return false if the file's n_sites header doesn't match this lattice's
/// current n_sites (stale/mismatched cache entry -- caller should treat as
/// a cache miss and re-relax rather than trust a truncated/misread file).
bool QfeLatticeS2::ReadPositions(FILE* file) {
  int32_t n = 0;
  if (fread(&n, sizeof(n), 1, file) != 1) return false;
  if (n != n_sites) return false;
  size_t got = fread(r.data(), sizeof(Vec3), n_sites, file);
  return got == size_t(n_sites);
}

/// @brief Iteratively relax vertex positions on the sphere to reduce the
/// variance of each site's *dual* (Voronoi-proxy) area, as opposed to
/// EqualizeFaceAreas's triangle-face-area target. Dual area per site is
/// approximated as the barycentric proxy sum_{f adjacent to s} FlatArea(f)/3
/// (standard 1/3-per-vertex split, avoids re-deriving the cotangent-FEM
/// weight's gradient). A site with above-mean dual area is pulled toward
/// its link-neighbor centroid (shrinking its incident faces); a site with
/// below-mean dual area is pushed away (growing them). Same
/// symmetry-preservation argument as EqualizeFaceAreas applies: a fixed,
/// purely local update rule applied identically to every site preserves
/// the starting mesh's icosahedral symmetry. This is the quantity that
/// actually sets each site's effective quadrature weight (sites[s].wt from
/// UpdateWeights/OptimizeIntegrator), so it is a different equalization
/// target from EqualizeFaceAreas even though both start from the same
/// naive mesh. See journal.md 2026-08-24.
/// @param n_iter Number of relaxation iterations
/// @param step Step size (fraction of the area-weighted displacement
/// applied per iteration)
void QfeLatticeS2::EqualizeDualAreas(int n_iter, double step) {
  for (int iter = 0; iter < n_iter; iter++) {
    std::vector<double> face_area(n_faces);
    for (int f = 0; f < n_faces; f++) {
      Vec3 v0 = r[faces[f].sites[0]];
      Vec3 v1 = r[faces[f].sites[1]];
      Vec3 v2 = r[faces[f].sites[2]];
      face_area[f] = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    }

    std::vector<double> dual_area(n_sites, 0.0);
    for (int f = 0; f < n_faces; f++) {
      double share = face_area[f] / 3.0;
      for (int i = 0; i < 3; i++) {
        dual_area[faces[f].sites[i]] += share;
      }
    }
    double dual_area_mean = 0.0;
    for (int s = 0; s < n_sites; s++) dual_area_mean += dual_area[s];
    dual_area_mean /= double(n_sites);
    // characteristic edge length at this resolution -- same
    // resolution-independence reasoning as EqualizeFaceAreas.
    double length_scale = sqrt(dual_area_mean);

    std::vector<Vec3> neighbor_sum(n_sites, Vec3::Zero());
    std::vector<int> n_neighbor(n_sites, 0);
    for (int l = 0; l < n_links; l++) {
      int s0 = links[l].sites[0];
      int s1 = links[l].sites[1];
      neighbor_sum[s0] += r[s1];
      neighbor_sum[s1] += r[s0];
      n_neighbor[s0]++;
      n_neighbor[s1]++;
    }

    for (int s = 0; s < n_sites; s++) {
      if (n_neighbor[s] == 0) continue;
      Vec3 neighbor_centroid = neighbor_sum[s] / double(n_neighbor[s]);
      Vec3 dir = neighbor_centroid - r[s];
      double dir_norm = dir.norm();
      if (dir_norm > 1e-14) dir /= dir_norm;
      double rel_imbalance = (dual_area[s] - dual_area_mean) / dual_area_mean;
      Vec3 d = rel_imbalance * dir;
      r[s] = (r[s] + step * length_scale * d).normalized();
    }
  }
}

/// @brief Calculate the global effective lattice spacing
/// @return Lattice spacing a/up
double QfeLatticeS2::CalcLatticeSpacing() {
  double area_sum = 0.0;
  for (int f = 0; f < n_faces; f++) {
    area_sum += FlatArea(f);
  }
  double area_mean = area_sum / double(n_sites);
  return sqrt(area_mean);
}

/// @brief Optimize the site weights so that all linear combinations of
/// spherical harmonics invariant under the relevant symmetry group can be
/// integrated exactly up to order @p l_max.
/// @param l_max Maximum spherical harmonic eigenvalue to optimize
void QfeLatticeS2::OptimizeIntegrator(int l_max) {
  // determine which l,m combinations need to be integrated exactly
  std::vector<int> l_relevant;
  std::vector<int> m_relevant;
  int m_spacing = (q == 5) ? 5 : 4;

  for (int l = 0; l <= l_max; l++) {
    // odd l only contributes for tetrahedron
    if ((l % 2) && (q != 3)) continue;

    // number of functions at this l (overcounts for tetrahedron)
    int n_l = (l / 2) + (l / 3) + (l / q) - l + 1;

    for (int i = 0; i < n_l; i++) {
      l_relevant.push_back(l);
      int m = i * m_spacing;

      // odd ell (tetrahedron only)
      if (l % 2) m = ((2 * i + 1) * m_spacing) / 2;

      m_relevant.push_back(m);
    }
  }

  // number of relevant functions
  int n_relevant = l_relevant.size();

  // generate the rectangular matrix
  ComplexMat S = ComplexMat::Zero(n_relevant, n_distinct);
  for (int s = 0; s < n_sites; s++) {
    int id = sites[s].id;
    double theta = acos(r[s].z());
    double phi = atan2(r[s].y(), r[s].x());
    for (int i = 0; i < n_relevant; i++) {
      int l = l_relevant[i];
      int m = m_relevant[i];
      S(i, id) += boost::math::spherical_harmonic(l, m, theta, phi);
    }
  }

  // use the current weights as an initial guess
  ComplexVec x0(n_distinct);
  for (int id = 0; id < n_distinct; id++) {
    int s = distinct_first[id];
    x0(id) = sites[s].wt;
  }

  // the right hand side is the spherical harmonic orthogonality condition
  ComplexVec b = ComplexVec::Zero(n_relevant);
  b(0) = 0.28209479177387814347 * vol;  // n_sites / sqrt(4 pi)

  // compute the solution
  Eigen::LeastSquaresConjugateGradient<ComplexMat> cg;
  cg.compute(S);
  assert(cg.info() == Eigen::Success);
  ComplexVec wt = cg.solveWithGuess(b, x0);

  // apply the improved weights to the sites
  for (int s = 0; s < n_sites; s++) {
    int id = sites[s].id;
    if (s == distinct_first[id]) {
      assert(real(wt(id)) > 0.0);
      // printf("%04d %.12f %.12f\n", id, sites[s].wt, real(wt(id)));
    }
    sites[s].wt = real(wt(id));
  }
}

/// @brief Update spherical harmonic values at each site, up to a maximum l
/// eigenvalue of @p l_max
/// @param l_max Maximum spherical harmonic eigenvalue to calculate
void QfeLatticeS2::UpdateYlm(int l_max) {
  int n_ylm = ((l_max + 1) * (l_max + 2)) / 2;
  using boost::math::spherical_harmonic;

  for (int s = 0; s < n_sites; s++) {
    ylm[s].resize(n_ylm);
    double theta = acos(r[s].z());
    double phi = atan2(r[s].y(), r[s].x());

    for (int i = 0, l = 0, m = 0; i < n_ylm; i++, m++) {
      if (m > l) {
        m = 0;
        l++;
      }
      assert(i < n_ylm);
      ylm[s][i] = spherical_harmonic(l, m, theta, phi);
    }
  }
}

/// @brief Retrieve a pre-calculated spherical harmonic value.
/// @param s Site index
/// @param l Spherical harmonic eigenvalue
/// @param m Spherical harmonic eigenvalue
/// @return Spherical harmonic evaluated at site @p s
Complex QfeLatticeS2::GetYlm(int s, int l, int m) {
  int abs_m = fabs(m);
  assert(abs_m <= l);

  int i = (l * (l + 1)) / 2 + abs_m;
  assert(i < ylm[s].size());
  Complex y = ylm[s][i];

  if (m < 0) {
    y = conj(y);
    if (abs_m & 1) {
      y *= -1;
    }
  }

  return y;
}

/// @brief Calculate a spherical harmonic value (not pre-calculated)
/// @param s Site index
/// @param l Spherical harmonic eigenvalue
/// @param m Spherical harmonic eigenvalue
/// @return Spherical harmonic evaluated at site @p s
Complex QfeLatticeS2::CalcYlm(int s, int l, int m) {
  double theta = acos(r[s].z());
  double phi = atan2(r[s].y(), r[s].x());

  return boost::math::spherical_harmonic(l, m, theta, phi);
}

/// @brief Calculate the cosine of the angle between two sites. This function
/// assumes that the coordinates have been projected onto the unit sphere.
/// @param s1 1st site index
/// @param s2 2nd site index
/// @return cosine of the angle between @p s1 and @p s2
double QfeLatticeS2::CosTheta(int s1, int s2) {
  if (s1 == s2) return 1.0;
  if (antipode[s1] == s2) return -1.0;
  return r[s1].dot(r[s2]);
}

/// @brief Print the cartesian coordinates of the sites. This is helpful for
/// making plots in e.g. Mathematica.
void QfeLatticeS2::PrintCoordinates() {
  printf("{");
  for (int s = 0; s < n_sites; s++) {
    printf("{%.12f,%.12f,%.12f}", r[s].x(), r[s].y(), r[s].z());
    printf("%c\n", s == (n_sites - 1) ? '}' : ',');
  }
}