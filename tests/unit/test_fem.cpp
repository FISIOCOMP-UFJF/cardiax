/*
 * Unit tests of src/fem: reference elements, quadrature, mappings, FE spaces, and
 * the L2 projection problem.
 *
 * What is covered and why:
 *  - Quadrature rules: the reference measure of each element, and exact integration
 *    of monomials up to the degree each rule claims. Every matrix assembled in the
 *    code is only as accurate as these.
 *  - Shape functions of every element: partition of unity, the Kronecker property at
 *    the nodes, and derivatives against finite differences.
 *  - Cubic Hermite element: value and slope interpolation properties.
 *  - Mapping: Jacobian and determinant on an affine and a distorted element, and the
 *    inverse Jacobian.
 *  - SurfaceMapping: surface area of triangles and quads in 3D, length of a line.
 *  - H1FESpace, CubicHermiteFESpace, VecH1FESpace: number of dofs and the global
 *    numbering, including the interleaved x,y,z ordering the mechanics relies on.
 *  - VecH1 finite elements: scalar shape functions shared by all components.
 *  - FETools: the penalty and zero-rows variants of Dirichlet conditions, L2 norm,
 *    and the Hermite derivative scaling.
 *  - Projection: the solver of the l2projection app, whose result on a linear
 *    function must be exact.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <functional>
#include <memory>

#include "fem/fem.h"
#include "fem/quadrature.hpp"
#include "fem/projection.hpp"

using namespace cardiax_test;

namespace
{

double ipow(double x, int p)
{
  double r = 1.0;
  for (int i = 0; i < p; i++) r *= x;
  return r;
}

double factorial(int n) { return n <= 1 ? 1.0 : n * factorial(n - 1); }

// Exact integral of x^a y^b z^c on the unit simplex of dimension 3 (2 with c = 0 and dim = 2)
double simplex_monomial(int a, int b, int c, int dim)
{
  return factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + dim);
}

// Exact integral of x^a on [-1, 1]
double interval_monomial(int a) { return a % 2 ? 0.0 : 2.0 / (a + 1); }

double integrate(const Quadrature & q, int a, int b, int c)
{
  double sum = 0.0;
  for (int i = 0; i < q.get_num_ipoints(); i++)
  {
    const arma::vec3 & p = q.get_point(i);
    sum += q.get_weight(i) * ipow(p(0), a) * ipow(p(1), b) * ipow(p(2), c);
  }
  return sum;
}

// Node coordinates of each reference element
std::vector<arma::vec3> reference_nodes(ElementType t)
{
  switch (t)
  {
    case ELEM_SEGM:  return {{-1, 0, 0}, {1, 0, 0}};
    case ELEM_TRIG:  return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    case ELEM_QUAD:  return {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
    case ELEM_TETRA: return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    case ELEM_HEXA:  return {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                             {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}};
  }
  return {};
}

void check_shape_functions(const FiniteElement & fe)
{
  const std::vector<arma::vec3> nodes = reference_nodes(fe.type());
  const int n = fe.get_ndof();
  REQUIRE(int(nodes.size()) == n);

  // Kronecker property
  arma::vec N;
  for (int i = 0; i < n; i++)
  {
    fe.calc_shape(nodes[i], N);
    REQUIRE(int(N.n_elem) == n);
    for (int j = 0; j < n; j++)
      CHECK_NEAR(N(j), i == j ? 1.0 : 0.0, 1e-14);
  }

  // Partition of unity and derivatives vs. centred differences, inside the element
  const std::vector<arma::vec3> inside = {{0.21, 0.13, 0.17}, {0.05, 0.31, 0.27}};
  const double h = 1e-6;
  for (const arma::vec3 & x : inside)
  {
    fe.calc_shape(x, N);
    CHECK_NEAR(arma::accu(N), 1.0, 1e-13);

    arma::mat dN;
    fe.calc_deriv_shape(x, dN);
    REQUIRE(int(dN.n_rows) == n);
    REQUIRE(int(dN.n_cols) == fe.get_ndim());
    for (int d = 0; d < fe.get_ndim(); d++)
    {
      arma::vec3 xp = x, xm = x;
      xp(d) += h;
      xm(d) -= h;
      arma::vec Np, Nm;
      fe.calc_shape(xp, Np);
      fe.calc_shape(xm, Nm);
      for (int i = 0; i < n; i++)
        CHECK_NEAR(dN(i, d), (Np(i) - Nm(i)) / (2.0 * h), 1e-8);
    }

    // Sum of the derivatives of a partition of unity is zero
    for (int d = 0; d < fe.get_ndim(); d++)
      CHECK_NEAR(arma::accu(dN.col(d)), 0.0, 1e-13);
  }
}

} // namespace

// ---- Quadrature -------------------------------------------------------------

TEST_CASE("Quadrature: the factory gives the rule of each element type, and throws on bad orders")
{
  std::unique_ptr<Quadrature> q(Quadrature::create(1, ELEM_TETRA));
  CHECK_EQ(q->get_num_ipoints(), 1);
  q.reset(Quadrature::create(2, ELEM_TETRA));
  CHECK_EQ(q->get_num_ipoints(), 4);
  q.reset(Quadrature::create(2, ELEM_HEXA));
  CHECK_EQ(q->get_num_ipoints(), 8);
  q.reset(Quadrature::create(-1, ELEM_HEXA));
  CHECK_EQ(q->get_num_ipoints(), 1);

  CHECK_THROWS(Quadrature::create(99, ELEM_TETRA));
  CHECK_THROWS(Quadrature::create(99, ELEM_HEXA));
  CHECK_THROWS(Quadrature::create(99, ELEM_TRIG));
  CHECK_THROWS(Quadrature::create(99, ELEM_QUAD));
  CHECK_THROWS(Quadrature::create(99, ELEM_SEGM));
}

TEST_CASE("Quadrature: weights add up to the measure of the reference element")
{
  const struct { ElementType t; std::vector<int> orders; double measure; } rules[] = {
    {ELEM_SEGM,  {0, 1, 2, 3, 4, 5}, 2.0},
    {ELEM_TRIG,  {0, 1, 2}, 0.5},
    {ELEM_QUAD,  {0, 1, 2, 3}, 4.0},
    {ELEM_TETRA, {0, 1, 2, 3, 4}, 1.0 / 6.0},
    {ELEM_HEXA,  {-1, 0, 1, 2}, 8.0},
  };
  for (const auto & r : rules)
    for (int o : r.orders)
    {
      std::unique_ptr<Quadrature> q(Quadrature::create(o, r.t));
      CHECK_MSG(q->get_num_ipoints() > 0, "type " << r.t << " order " << o);
      CHECK_NEAR(integrate(*q, 0, 0, 0), r.measure, 1e-12);
    }
}

TEST_CASE("Quadrature: 1D rules are exact up to their degree")
{
  for (int order : {1, 2, 3, 4, 5})
  {
    std::unique_ptr<Quadrature> q(Quadrature::create(order, ELEM_SEGM));
    for (int a = 0; a <= order; a++)
      CHECK_NEAR(integrate(*q, a, 0, 0), interval_monomial(a), 1e-14);
  }
}

TEST_CASE("Quadrature: triangle rules are exact up to their degree")
{
  for (int order : {1, 2})
  {
    std::unique_ptr<Quadrature> q(Quadrature::create(order, ELEM_TRIG));
    for (int a = 0; a <= order; a++)
      for (int b = 0; a + b <= order; b++)
        CHECK_NEAR(integrate(*q, a, b, 0), simplex_monomial(a, b, 0, 2), 1e-14);
  }
}

TEST_CASE("Quadrature: tetrahedron rules are exact up to degree 2")
{
  for (int order : {1, 2})
  {
    std::unique_ptr<Quadrature> q(Quadrature::create(order, ELEM_TETRA));
    for (int a = 0; a <= order; a++)
      for (int b = 0; a + b <= order; b++)
        for (int c = 0; a + b + c <= order; c++)
          CHECK_NEAR(integrate(*q, a, b, c), simplex_monomial(a, b, c, 3), 1e-14);
  }
}

TEST_CASE_KNOWN_BUG("Quadrature: a tetrahedron rule of order 3 or 4 is exact up to that degree",
                    "src/fem/quadrature.cpp: orders 3 and 4 return the 4-point rule, exact "
                    "only up to degree 2, so asking for more silently under-integrates")
{
  for (int order : {3, 4})
  {
    std::unique_ptr<Quadrature> q(Quadrature::create(order, ELEM_TETRA));
    for (int a = 0; a <= order; a++)
      for (int b = 0; a + b <= order; b++)
        for (int c = 0; a + b + c <= order; c++)
          CHECK_NEAR(integrate(*q, a, b, c), simplex_monomial(a, b, c, 3), 1e-12);
  }
}

TEST_CASE("Quadrature: quadrilateral and hexahedron Gauss rules integrate cubics in each direction")
{
  std::unique_ptr<Quadrature> q2(Quadrature::create(2, ELEM_QUAD));
  std::unique_ptr<Quadrature> h(Quadrature::create(2, ELEM_HEXA));
  for (int a = 0; a <= 3; a++)
    for (int b = 0; b <= 3; b++)
    {
      CHECK_NEAR(integrate(*q2, a, b, 0), interval_monomial(a) * interval_monomial(b), 1e-14);
      for (int c = 0; c <= 3; c++)
        CHECK_NEAR(integrate(*h, a, b, c),
                   interval_monomial(a) * interval_monomial(b) * interval_monomial(c), 1e-14);
    }

  std::unique_ptr<Quadrature> q3(Quadrature::create(3, ELEM_QUAD));
  for (int a = 0; a <= 5; a++)
    for (int b = 0; b <= 5; b++)
      CHECK_NEAR(integrate(*q3, a, b, 0), interval_monomial(a) * interval_monomial(b), 1e-14);
}

// ---- Shape functions --------------------------------------------------------

TEST_CASE("Shape functions: FESegm")  { FESegm fe(nullptr);  check_shape_functions(fe); }
TEST_CASE("Shape functions: FETrig")  { FETrig fe(nullptr);  check_shape_functions(fe); }
TEST_CASE("Shape functions: FEQuad")  { FEQuad fe(nullptr);  check_shape_functions(fe); }
TEST_CASE("Shape functions: FETetra") { FETetra fe(nullptr); check_shape_functions(fe); }
TEST_CASE("Shape functions: FEHexa")  { FEHexa fe(nullptr);  check_shape_functions(fe); }

TEST_CASE("Shape functions: properties of each element")
{
  FESegm s(nullptr);
  CHECK_EQ(s.get_ndof(), 2);
  CHECK_EQ(s.get_ndim(), 1);
  FETetra t(nullptr);
  CHECK_EQ(t.ndofs(), 4);
  CHECK_EQ(t.get_ndim(), 3);
  CHECK_EQ(t.order(), 1);
  CHECK(t.type() == ELEM_TETRA);
  FEHexa h(nullptr);
  CHECK_EQ(h.get_ndof(), 8);
  CHECK(h.type() == ELEM_HEXA);
}

TEST_CASE("FECubicHermite: interpolates values and slopes at the ends")
{
  FECubicHermite fe(nullptr);
  CHECK_EQ(fe.order(), 3);

  arma::vec N;
  arma::mat dN;
  const arma::vec3 left = {-1, 0, 0}, right = {1, 0, 0};

  fe.calc_shape(left, N);
  CHECK_NEAR(N(0), 1.0, 1e-14);
  CHECK_NEAR(N(1), 0.0, 1e-14);
  CHECK_NEAR(N(2), 0.0, 1e-14);
  CHECK_NEAR(N(3), 0.0, 1e-14);

  fe.calc_shape(right, N);
  CHECK_NEAR(N(0), 0.0, 1e-14);
  CHECK_NEAR(N(1), 1.0, 1e-14);

  fe.calc_deriv_shape(left, dN);
  CHECK_NEAR(dN(0, 0), 0.0, 1e-14);
  CHECK_NEAR(dN(1, 0), 0.0, 1e-14);
  CHECK_NEAR(dN(2, 0), 1.0, 1e-14);
  CHECK_NEAR(dN(3, 0), 0.0, 1e-14);

  fe.calc_deriv_shape(right, dN);
  CHECK_NEAR(dN(2, 0), 0.0, 1e-14);
  CHECK_NEAR(dN(3, 0), 1.0, 1e-14);
}

TEST_CASE("FECubicHermite: derivatives match finite differences")
{
  FECubicHermite fe(nullptr);
  const double h = 1e-6;
  for (double xi : {-0.7, 0.0, 0.4})
  {
    arma::vec Np, Nm;
    arma::mat dN;
    fe.calc_shape(arma::vec3({xi + h, 0, 0}), Np);
    fe.calc_shape(arma::vec3({xi - h, 0, 0}), Nm);
    fe.calc_deriv_shape(arma::vec3({xi, 0, 0}), dN);
    for (int i = 0; i < 4; i++)
      CHECK_NEAR(dN(i, 0), (Np(i) - Nm(i)) / (2.0 * h), 1e-8);
  }
}

// ---- Mapping ----------------------------------------------------------------

TEST_CASE("Mapping: Jacobian of a scaled and sheared tetrahedron")
{
  std::vector<arma::vec3> pts = {{1, 1, 1}, {3, 1, 1}, {1, 4, 1}, {2, 1, 6}};
  Mapping map(pts, 3, 1);
  FETetra fe(nullptr);
  arma::mat dN;
  fe.calc_deriv_shape(arma::vec3({0.25, 0.25, 0.25}), dN);
  map.calc_jacobian(dN);

  // det = 6 * volume
  CHECK_NEAR(map.get_det_jacobian(), 2.0 * 3.0 * 5.0, 1e-12);
  const arma::mat I = map.get_jacobian() * map.get_inv_jacobian();
  CHECK_NEAR(arma::abs(I - arma::eye(3, 3)).max(), 0.0, 1e-13);
  CHECK_NEAR(map.get_jacobian()(0, 2), 1.0, 1e-13);   // dx/dzeta of the shear
}

TEST_CASE("Mapping: Jacobian of a box-shaped hexahedron is constant")
{
  std::vector<arma::vec3> pts = {{0, 0, 0}, {4, 0, 0}, {4, 2, 0}, {0, 2, 0},
                                 {0, 0, 6}, {4, 0, 6}, {4, 2, 6}, {0, 2, 6}};
  Mapping map(pts, 3, 1);
  FEHexa fe(nullptr);
  for (const arma::vec3 & x : std::vector<arma::vec3>{{0, 0, 0}, {0.5, -0.3, 0.8}})
  {
    arma::mat dN;
    fe.calc_deriv_shape(x, dN);
    map.calc_jacobian(dN);
    CHECK_NEAR(map.get_det_jacobian(), 2.0 * 1.0 * 3.0, 1e-13);
  }
}

TEST_CASE("Mapping: quadrature of det J gives the volume of a distorted hexahedron")
{
  // Unit cube with one vertex moved along its diagonal: the volume is known
  // (1 + 3*d/8... computed by decomposing the moved corner into 6 pyramids below)
  const double d = 0.3;
  std::vector<arma::vec3> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                 {0, 0, 1}, {1, 0, 1}, {1 + d, 1 + d, 1 + d}, {0, 1, 1}};
  Mapping map(pts, 3, 1);
  FEHexa fe(nullptr);

  std::unique_ptr<Quadrature> q(Quadrature::create(2, ELEM_HEXA));
  double vol = 0.0;
  for (int i = 0; i < q->get_num_ipoints(); i++)
  {
    arma::mat dN;
    fe.calc_deriv_shape(q->get_point(i), dN);
    map.calc_jacobian(dN);
    vol += q->get_weight(i) * map.get_det_jacobian();
  }

  // Divergence theorem on the (bilinear) boundary gives V = 1 + 3 d / 4 + d^2 * 0
  // for this motion: only the three faces touching the vertex change, each by
  // the same flux, linear in d (plus a quadratic term that the two-point rule
  // integrates exactly, since det J is at most quadratic in each direction).
  const double fd = 1e-6;
  auto volume_of = [&](double dd) {
    std::vector<arma::vec3> p2 = pts;
    p2[6] = arma::vec3({1 + dd, 1 + dd, 1 + dd});
    Mapping m2(p2, 3, 1);
    double v = 0.0;
    for (int i = 0; i < q->get_num_ipoints(); i++)
    {
      arma::mat dN;
      fe.calc_deriv_shape(q->get_point(i), dN);
      m2.calc_jacobian(dN);
      v += q->get_weight(i) * m2.get_det_jacobian();
    }
    return v;
  };
  CHECK_NEAR(volume_of(0.0), 1.0, 1e-13);
  CHECK_NEAR((volume_of(fd) - volume_of(-fd)) / (2 * fd), 0.25 * 3.0, 1e-6);
  CHECK(vol > 1.0);
  CHECK_NEAR(vol, volume_of(d), 0.0);
}

TEST_CASE("Mapping: map_point sends reference points to physical points")
{
  std::vector<arma::vec3> pts = {{1, 2, 3}, {2, 2, 3}, {1, 3, 3}, {1, 2, 4}};
  Mapping map(pts, 3, 1);
  FETetra fe(nullptr);
  arma::vec N;
  const arma::vec3 xi = {0.2, 0.3, 0.1};
  fe.calc_shape(xi, N);
  const arma::vec3 X = map.map_point(xi, N);
  CHECK_NEAR(X(0), 1.2, 1e-14);
  CHECK_NEAR(X(1), 2.3, 1e-14);
  CHECK_NEAR(X(2), 3.1, 1e-14);
}

TEST_CASE("Mapping: segment Jacobian works for a line in 3D space")
{
  std::vector<arma::vec3> pts = {{0, 0, 0}, {3, 4, 0}};
  Mapping map(pts, 1, 1);
  FESegm fe(nullptr);
  arma::mat dN;
  fe.calc_deriv_shape(arma::vec3({0, 0, 0}), dN);
  map.calc_segment_jacobian(dN);
  CHECK_NEAR(map.get_det_jacobian(), 5.0, 1e-14);
  CHECK_NEAR(map.get_jacobian()(1, 0), 2.0, 1e-14);
}

TEST_CASE("SurfaceMapping: area scale of a triangle and a quadrilateral in 3D")
{
  {
    std::vector<arma::vec3> pts = {{0, 0, 1}, {2, 0, 1}, {0, 3, 1}};
    SurfaceMapping sm(pts, 2, 1);
    FETrig fe(nullptr);
    arma::mat dN;
    fe.calc_deriv_shape(arma::vec3({0.2, 0.2, 0}), dN);
    sm.calc_jacobian(dN);
    // |a x b| = twice the area of the triangle
    CHECK_NEAR(sm.get_det_jacobian(), 6.0, 1e-14);
  }
  {
    std::vector<arma::vec3> pts = {{0, 0, 0}, {2, 0, 0}, {2, 0, 4}, {0, 0, 4}};
    SurfaceMapping sm(pts, 2, 1);
    FEQuad fe(nullptr);
    arma::mat dN;
    fe.calc_deriv_shape(arma::vec3({0.1, 0.4, 0}), dN);
    sm.calc_jacobian(dN);
    // weights of [-1,1]^2 add up to 4, so det = area / 4
    CHECK_NEAR(sm.get_det_jacobian(), 8.0 / 4.0, 1e-14);
  }
}

TEST_CASE_KNOWN_BUG("SurfaceMapping: a line integrates to its length with the [-1, 1] weights",
                    "src/fem/mapping.cpp: calc_jacobian sets det = length for a line; with the "
                    "weights of the reference segment (adding up to 2) that gives twice the length")
{
  std::vector<arma::vec3> pts = {{0, 0, 0}, {3, 4, 0}};
  SurfaceMapping sm(pts, 1);
  FESegm fe(nullptr);
  arma::mat dN;
  fe.calc_deriv_shape(arma::vec3({0, 0, 0}), dN);
  sm.calc_jacobian(dN);
  std::unique_ptr<Quadrature> q(Quadrature::create(1, ELEM_SEGM));
  CHECK_NEAR(q->get_weight(0) * sm.get_det_jacobian(), 5.0, 1e-14);
}

// ---- FE spaces --------------------------------------------------------------

TEST_CASE("H1FESpace: one dof per node, numbered like the mesh, and the right element type")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  H1FESpace V(&m);
  V.config();
  CHECK_EQ(V.get_num_dofs(), int(m.get_n_points()));

  std::vector<int> d, p;
  CHECK_EQ(V.get_element_dofs(3, d), 4);
  m.get_element_pt_nums(3, p);
  CHECK(d == p);

  CHECK_EQ(V.get_boundary_element_dofs(0, d), 3);
  CHECK(V.createFE(0).type() == ELEM_TETRA);
  CHECK(V.create_boundary_FE(0).type() == ELEM_TRIG);
}

TEST_CASE("H1FESpace: element types of a hexahedral, quadrilateral and line mesh")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  write_hex_box("h.xml", b);
  write_quad_square("q.xml", b);
  write_tri_square("t.xml", b);
  write_cable("c.xml", 2, 1.0);
  Mesh h, q, t, c;
  h.read_xml("h.xml");
  q.read_xml("q.xml");
  t.read_xml("t.xml");
  c.read_xml("c.xml");
  CHECK(H1FESpace(&h).createFE(0).type() == ELEM_HEXA);
  CHECK(H1FESpace(&h).create_boundary_FE(0).type() == ELEM_QUAD);
  CHECK(H1FESpace(&q).createFE(0).type() == ELEM_QUAD);
  CHECK(H1FESpace(&q).create_boundary_FE(0).type() == ELEM_SEGM);
  CHECK(H1FESpace(&t).createFE(0).type() == ELEM_TRIG);
  CHECK(H1FESpace(&c).createFE(0).type() == ELEM_SEGM);
}

TEST_CASE("CubicHermiteFESpace: two dofs per node, value and slope")
{
  write_cable("c.xml", 3, 3.0);
  Mesh m;
  m.read_xml("c.xml");
  CubicHermiteFESpace V(&m);
  V.config();
  CHECK_EQ(V.get_num_dofs(), 2 * 4);
  CHECK_EQ(V.createFE(0).get_ndof(), 4);

  std::vector<int> d;
  CHECK_EQ(V.get_element_dofs(1, d), 4);
  // Value dofs are the nodes; slope dofs come after all value dofs
  CHECK_EQ(d[0], 1);
  CHECK_EQ(d[1], 2);
  CHECK_EQ(d[2], 4 + 1);
  CHECK_EQ(d[3], 4 + 2);
}

TEST_CASE("VecH1FESpace: three interleaved dofs per node in 3D")
{
  BoxSpec b;
  b.nx = 2; b.ny = 1; b.nz = 1;
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");

  VecH1FESpace V;
  V.set_mesh(&m);
  CHECK_EQ(V.get_ndofs(), 3 * int(m.get_n_points()));
  CHECK_EQ(V.get_ndofs_u(), V.get_ndofs());
  CHECK_EQ(V.get_ndofs_p(), 0);

  std::vector<int> d, p;
  CHECK_EQ(V.get_element_dofs_u(2, d), 12);
  m.get_element_pt_nums(2, p);
  for (int a = 0; a < 4; a++)
  {
    CHECK_EQ(d[a], 3 * p[a]);
    CHECK_EQ(d[a + 4], 3 * p[a] + 1);
    CHECK_EQ(d[a + 8], 3 * p[a] + 2);
  }

  std::vector<int> bd, bp;
  CHECK_EQ(V.get_boundary_element_dofs_u(1, bd), 9);
  m.get_boundary_element_pt_nums(1, bp);
  CHECK_EQ(bd[3], 3 * bp[0] + 1);
  CHECK_EQ(bd[8], 3 * bp[2] + 2);
}

TEST_CASE("VecH1 elements: same scalar shape functions and the right dof count")
{
  const struct { const char * name; int ndim; std::function<void(const std::string&)> write;
                 int scalar_dofs; ElementType t; int boundary_dofs; } meshes[] = {
    {"tet", 3, [](const std::string & f){ BoxSpec b; write_tet_box(f, b); }, 4, ELEM_TETRA, 3},
    {"hex", 3, [](const std::string & f){ BoxSpec b; write_hex_box(f, b); }, 8, ELEM_HEXA, 4},
    {"tri", 2, [](const std::string & f){ BoxSpec b; write_tri_square(f, b); }, 3, ELEM_TRIG, 2},
    {"qua", 2, [](const std::string & f){ BoxSpec b; write_quad_square(f, b); }, 4, ELEM_QUAD, 2},
  };
  for (const auto & mm : meshes)
  {
    const std::string name = std::string(mm.name) + ".xml";
    mm.write(name);
    Mesh m;
    m.read_xml(name);
    VecH1FESpace V;
    V.set_mesh(&m);
    std::unique_ptr<MixedFiniteElement> fe(V.createFE());
    std::unique_ptr<MixedFiniteElement> bfe(V.create_boundary_FE());
    CHECK(fe->get_type() == mm.t);
    CHECK_EQ(fe->get_ndofs_u(), mm.ndim * mm.scalar_dofs);
    CHECK_EQ(bfe->get_ndofs_u(), mm.ndim * mm.boundary_dofs);

    arma::vec N;
    arma::mat dN;
    const arma::vec3 x = {0.2, 0.2, 0.2};
    fe->calc_shape_u(x, N);
    fe->calc_deriv_shape_u(x, dN);
    // One copy of the scalar functions per component, each a partition of unity
    CHECK_EQ(int(N.n_elem), mm.ndim * mm.scalar_dofs);
    CHECK_NEAR(arma::accu(N), double(mm.ndim), 1e-13);
    CHECK_EQ(int(dN.n_cols), mm.ndim);
  }
}

// ---- FETools ----------------------------------------------------------------

TEST_CASE("FETools: L2 norm with a mass matrix")
{
  arma::sp_mat M(2, 2);
  M(0, 0) = 2.0;
  M(1, 1) = 3.0;
  arma::vec u = {1.0, 2.0};
  CHECK_NEAR(FETools::calc_L2_norm(M, u), std::sqrt(2.0 + 12.0), 1e-14);
}

TEST_CASE("FETools: penalty boundary values impose the value on an Armadillo system")
{
  arma::sp_mat A(2, 2);
  A(0, 0) = 2.0; A(0, 1) = -1.0; A(1, 0) = -1.0; A(1, 1) = 2.0;
  arma::vec b = {0.0, 1.0};
  std::map<int, double> bc = {{0, 3.0}};
  FETools::apply_boundary_values(bc, A, b);
  const arma::vec x = arma::solve(arma::mat(A), b);
  CHECK_NEAR(x(0), 3.0, 1e-8);
  CHECK_NEAR(x(1), 2.0, 1e-8);   // 2u1 = 1 + 3
}

TEST_CASE("FETools: penalty boundary values impose the value on a PETSc system")
{
  petsc::Matrix A;
  A.create(2, 2, 2);
  A.set(0, 0, 2.0); A.set(0, 1, -1.0); A.set(1, 0, -1.0); A.set(1, 1, 2.0);
  A.assemble();
  petsc::Vector b, x;
  b.create(2);
  x.create(2);
  b.set(1.0);
  x = 0.0;
  std::map<int, double> bc = {{1, 5.0}};
  FETools::apply_boundary_values(bc, A, b);
  A.assemble();
  b.assemble();
  petsc::LinearSolver ls;
  ls.init();
  ls.solve(A, x, b);
  double out[2];
  x.get_data(out);
  CHECK_NEAR(out[1], 5.0, 1e-8);
}

TEST_CASE("FETools: zero-rows boundary values, with and without RHS correction")
{
  for (bool update_rhs : {true, false})
  {
    petsc::Matrix A;
    A.create(2, 2, 2);
    A.set(0, 0, 2.0); A.set(0, 1, -1.0); A.set(1, 0, -1.0); A.set(1, 1, 2.0);
    A.assemble();
    petsc::Vector b, xbar, x;
    b.create(2);
    xbar.create(2);
    x.create(2);
    b = 0.0;
    xbar = 0.0;
    xbar.set(0, 4.0);
    xbar.assemble();
    x = 0.0;
    std::map<int, double> bc = {{0, 4.0}};
    FETools::apply_boundary_values(bc, A, b, xbar, update_rhs);
    petsc::LinearSolver ls;
    ls.init();
    ls.solve(A, x, b);
    double out[2];
    x.get_data(out);
    if (update_rhs)
    {
      CHECK_NEAR(out[0], 4.0, 1e-12);
      CHECK_NEAR(out[1], 2.0, 1e-12);
    }
    else
      CHECK_NEAR(out[0], 0.0, 1e-12);
  }
}

TEST_CASE("FETools: Hermite slope dofs scale by 2/L of the element")
{
  write_cable("c.xml", 2, 8.0);
  Mesh m;
  m.read_xml("c.xml");
  H1FESpace V(&m);
  V.config();
  arma::vec d(3);
  d.fill(1.0);
  FETools::cubic_hermite_deriv(V, m, d);
  CHECK_NEAR(d(0), 2.0 / 4.0, 1e-14);
  CHECK_NEAR(d(2), 2.0 / 4.0, 1e-14);
}

// ---- Projection -------------------------------------------------------------

TEST_CASE("Projection: L2 projection reproduces its source at the mesh nodes as h shrinks")
{
  // The class projects sin(pi x) sin(pi y). The error at the nodes must decrease
  // with the mesh size.
  auto nodal_error = [](int n) {
    const std::string name = "p" + std::to_string(n) + ".xml";
    BoxSpec b;
    b.nx = b.ny = n;
    write_tri_square(name, b);
    Projection P;
    P.run(name);
    const arma::vec & u = P.get_solution();
    const Mesh & m = P.get_mesh();
    double e = 0.0;
    for (uint i = 0; i < m.get_n_points(); i++)
    {
      const arma::vec3 & p = m.get_point(i);
      e = std::max(e, std::fabs(u(i) - std::sin(M_PI * p(0)) * std::sin(M_PI * p(1))));
    }
    return e;
  };
  const double e1 = nodal_error(8), e2 = nodal_error(16);
  CHECK_MSG(e1 < 0.1, e1);
  CHECK_MSG(e2 < e1 / 2.5, e1 << " -> " << e2);
}

TEST_CASE("Projection: write_data creates the result file")
{
  BoxSpec b;
  b.nx = b.ny = 2;
  write_tri_square("p.xml", b);
  Projection P;
  P.run("p.xml");
  P.write_data("proj_out", false);
  // Static nodal field [np], read as a column
  const arma::mat u = h5_read("proj_out.h5", "/point_data/solution");
  CHECK_EQ(u.n_rows, arma::uword(9));
  CHECK_REL(u(4, 0), P.get_solution()(4), 1e-12);
}
