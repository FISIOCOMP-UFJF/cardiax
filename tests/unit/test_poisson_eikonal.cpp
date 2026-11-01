/*
 * Unit tests of the scalar PDE solvers in src/pdes: Poisson, Laplace and Eikonal.
 *
 * What is covered and why:
 *  - Element matrices of the Poisson problem, reached through a subclass that exposes
 *    them: the stiffness of the reference tetrahedron has a known value, and every
 *    element stiffness has the constant vector in its null space and is symmetric. The
 *    mass matrix adds up to the volume.
 *  - Solutions the discretisation represents exactly: a linear solution with Dirichlet
 *    and Neumann conditions. This covers boundary-condition assembly on tets, triangles,
 *    quads and hexahedra.
 *  - Rate of convergence against a manufactured solution, the check that catches an
 *    assembly mistake that a single fixed solution could hide.
 *  - Reading conditions and source term from TOML, the only input path of the app.
 *  - Laplace: a linear field between two Dirichlet walls.
 *  - Eikonal: times along mesh edges equal distance/velocity, roots keep their time,
 *    several roots activate the nearest one first, and anisotropy follows the fibres.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include "pdes/poisson.hpp"
#include "pdes/laplace.hpp"
#include "pdes/eikonal.hpp"
#include "util/expr_function.hpp"

using namespace cardiax_test;

namespace
{

struct OpenPoisson : public Poisson
{
  using Poisson::calc_elmat_poisson;
  using Poisson::calc_elmat_mass;
  using Poisson::fespace;
  using Poisson::mesh;

  void load(const std::string & file)
  {
    mesh.read_xml(file);
    fespace.set_mesh(&mesh);
    fespace.config();
  }
};

toml::table poisson_cfg(const std::string & body)
{
  return toml::parse(body);
}

// Solves and returns the L2 error
double solve_l2(const std::string & mesh, const std::string & cfg, const std::string & exact)
{
  Poisson p;
  p.config(poisson_cfg(cfg));
  p.run(mesh);
  ExprScalarFunction ex(exact);
  return p.calc_l2_error(ex);
}

const char * kLinearCfg =
  "[physical]\n"
  "source = \"0\"\n"
  "dirichlet = [ { marker = 1, value = 0.0 } ]\n"
  "neumann   = [ { marker = 2, value = 1.0 }, { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 },\n"
  "              { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n";

const char * kParabolaCfg =
  "[physical]\n"
  "source = \"2\"\n"
  "dirichlet = [ { marker = 1, value = 0.0 }, { marker = 2, value = 0.0 } ]\n"
  "neumann   = [ { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 },\n"
  "              { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n";

} // namespace

// ---- Element matrices -------------------------------------------------------

TEST_CASE("Poisson: stiffness of the reference tetrahedron")
{
  write_text("t.xml",
    "<mesh celltype=\"tetrahedron\" dim=\"3\"><nodes size=\"4\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"0\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"0\" z=\"1\"/>"
    "</nodes><elements size=\"1\"><element id=\"0\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\"/>"
    "</elements></mesh>");
  OpenPoisson p;
  p.load("t.xml");
  FiniteElement & fe = p.fespace.createFE(0);
  arma::mat K(4, 4);
  p.calc_elmat_poisson(0, fe, K);

  const arma::mat expected = (1.0 / 6.0) * arma::mat({{3, -1, -1, -1},
                                                      {-1, 1, 0, 0},
                                                      {-1, 0, 1, 0},
                                                      {-1, 0, 0, 1}});
  CHECK_NEAR(arma::abs(K - expected).max(), 0.0, 1e-14);
  delete &fe;
}

TEST_CASE("Poisson: mass of the reference tetrahedron")
{
  write_text("t.xml",
    "<mesh celltype=\"tetrahedron\" dim=\"3\"><nodes size=\"4\">"
    "<node id=\"0\" x=\"0\" y=\"0\" z=\"0\"/><node id=\"1\" x=\"1\" y=\"0\" z=\"0\"/>"
    "<node id=\"2\" x=\"0\" y=\"1\" z=\"0\"/><node id=\"3\" x=\"0\" y=\"0\" z=\"1\"/>"
    "</nodes><elements size=\"1\"><element id=\"0\" v0=\"0\" v1=\"1\" v2=\"2\" v3=\"3\"/>"
    "</elements></mesh>");
  OpenPoisson p;
  p.load("t.xml");
  FiniteElement & fe = p.fespace.createFE(0);
  arma::mat M(4, 4);
  p.calc_elmat_mass(0, fe, M);

  // Exact P1 mass: V/20 * (1 + delta_ij), with V = 1/6
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      CHECK_NEAR(M(i, j), (1.0 / 6.0) / 20.0 * (i == j ? 2.0 : 1.0), 1e-14);
  delete &fe;
}

TEST_CASE("Poisson: element matrices of every element type are symmetric, with constants in the kernel")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.Lx = 2.0; b.Ly = 3.0; b.Lz = 5.0;
  write_tet_box("tet.xml", b);
  write_hex_box("hex.xml", b);
  write_tri_square("tri.xml", b);
  write_quad_square("qua.xml", b);

  for (const char * name : {"tet.xml", "hex.xml", "tri.xml", "qua.xml"})
  {
    OpenPoisson p;
    p.load(name);
    FiniteElement & fe = p.fespace.createFE(0);
    const int n = fe.get_ndof();
    double mass = 0.0;
    for (int e = 0; e < p.mesh.get_n_elements(); e++)
    {
      arma::mat K(n, n), M(n, n);
      p.calc_elmat_poisson(e, fe, K);
      p.calc_elmat_mass(e, fe, M);
      CHECK_MSG(arma::abs(K - K.t()).max() < 1e-12, name);
      CHECK_MSG(arma::abs(M - M.t()).max() < 1e-12, name);
      CHECK_MSG(arma::abs(K * arma::ones<arma::vec>(n)).max() < 1e-12, name);
      mass += arma::accu(M);
    }
    const double domain = (p.mesh.get_n_dim() == 3) ? 2.0 * 3.0 * 5.0 : 2.0 * 3.0;
    CHECK_MSG(std::fabs(mass - domain) < 1e-11, name << " " << mass);
    delete &fe;
  }
}

// ---- Exact solutions --------------------------------------------------------

TEST_CASE("Poisson: linear solution with Dirichlet and Neumann, tetrahedra")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 3;
  write_tet_box("m.xml", b);
  CHECK_NEAR(solve_l2("m.xml", kLinearCfg, "x"), 0.0, 1e-7);
}

TEST_CASE("Poisson: linear solution with Dirichlet and Neumann, hexahedra")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 3;
  write_hex_box("m.xml", b);
  CHECK_NEAR(solve_l2("m.xml", kLinearCfg, "x"), 0.0, 1e-7);
}

TEST_CASE_KNOWN_BUG("Poisson: linear solution with Dirichlet and Neumann, triangles",
                    "src/fem/mapping.cpp: the SurfaceMapping of a line returns its length as det J, "
                    "but the reference segment is [-1, 1] (det J = length/2): every boundary "
                    "integral in 2D is twice too large, and here u = 2x comes out")
{
  BoxSpec b;
  b.nx = b.ny = 4;
  write_tri_square("m.xml", b);
  CHECK_NEAR(solve_l2("m.xml", kLinearCfg, "x"), 0.0, 1e-7);
}

TEST_CASE_KNOWN_BUG("Poisson: linear solution with Dirichlet and Neumann, quadrilaterals",
                    "src/fem/mapping.cpp: the SurfaceMapping of a line returns its length as det J, "
                    "but the reference segment is [-1, 1] (det J = length/2): every boundary "
                    "integral in 2D is twice too large, and here u = 2x comes out")
{
  BoxSpec b;
  b.nx = b.ny = 4;
  write_quad_square("m.xml", b);
  CHECK_NEAR(solve_l2("m.xml", kLinearCfg, "x"), 0.0, 1e-7);
}

TEST_CASE("Poisson: a nonzero Dirichlet value is imposed")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  write_tet_box("m.xml", b);
  const std::string cfg =
    "[physical]\nsource = \"0\"\n"
    "dirichlet = [ { marker = 1, value = 3.0 } ]\n"
    "neumann = [ { marker = 2, value = 0.0 }, { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 },\n"
    "            { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n";
  CHECK_NEAR(solve_l2("m.xml", cfg, "3"), 0.0, 1e-5);
}

TEST_CASE("Poisson: fixed nodes pin the solution at the nodes")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  write_tet_box("m.xml", b);
  const int corner = box_node(b, 0, 0, 0);
  const std::string cfg =
    "[physical]\nsource = \"0\"\n"
    "fix_node = [ { marker = " + std::to_string(corner) + ", value = 2.5 } ]\n"
    "neumann = [ { marker = 1, value = 0.0 }, { marker = 2, value = 0.0 }, { marker = 3, value = 0.0 },\n"
    "            { marker = 4, value = 0.0 }, { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n";
  Poisson p;
  p.config(poisson_cfg(cfg));
  p.run("m.xml");
  CHECK_NEAR(p.get_solution()(corner), 2.5, 1e-6);
  // No other load: the constant solution
  CHECK_NEAR(p.get_solution()(box_node(b, 2, 2, 2)), 2.5, 1e-6);
}

// ---- Convergence ------------------------------------------------------------

TEST_CASE("Poisson: second order in L2 for a parabola, on triangles")
{
  auto err = [](int n) {
    const std::string name = "t" + std::to_string(n) + ".xml";
    BoxSpec b;
    b.nx = n; b.ny = 2;
    write_tri_square(name, b);
    return solve_l2(name, kParabolaCfg, "x*(1-x)");
  };
  const double e1 = err(8), e2 = err(16);
  CHECK_MSG(e2 < e1, e1 << " " << e2);
  CHECK_NEAR(std::log2(e1 / e2), 2.0, 0.2);
}

TEST_CASE("Poisson: second order in L2 for a parabola, on tetrahedra")
{
  auto err = [](int n) {
    const std::string name = "m" + std::to_string(n) + ".xml";
    BoxSpec b;
    b.nx = n; b.ny = b.nz = 2;
    write_tet_box(name, b);
    return solve_l2(name, kParabolaCfg, "x*(1-x)");
  };
  const double e1 = err(8), e2 = err(16);
  CHECK_NEAR(std::log2(e1 / e2), 2.0, 0.2);
}

TEST_CASE("Poisson: calc_mass_matrix integrates to the volume, with the error norm it defines")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.Lx = 2.0;
  write_tet_box("m.xml", b);
  Poisson p;
  p.config(poisson_cfg(kLinearCfg));
  p.run("m.xml");
  arma::sp_mat M;
  p.calc_mass_matrix(M);
  CHECK_NEAR(arma::accu(M), 2.0, 1e-12);

  arma::vec one(p.get_mesh().get_n_points(), arma::fill::ones);
  CHECK_NEAR(FETools::calc_L2_norm(M, one), std::sqrt(2.0), 1e-12);
}

TEST_CASE("Poisson: write_data stores the solution")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  write_tet_box("m.xml", b);
  Poisson p;
  p.config(poisson_cfg(kLinearCfg));
  p.run("m.xml");
  p.write_data("poisson_out");
  const arma::mat vm = h5_read("poisson_out.h5", "/point_data/vm");
  CHECK_NEAR(vm(0, box_node(b, 2, 0, 0)), p.get_solution()(box_node(b, 2, 0, 0)), 0.0);
}

// ---- Laplace ----------------------------------------------------------------

TEST_CASE("Laplace: linear field between two walls")
{
  BoxSpec b;
  b.nx = 4; b.ny = b.nz = 2;
  write_tet_box("m.xml", b);
  const std::string cfg =
    "[physical]\n"
    "dirichlet = [ { marker = 1, value = 1.0 }, { marker = 2, value = 3.0 } ]\n"
    "neumann = [ { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 },\n"
    "            { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n";
  Laplace p;
  p.config(poisson_cfg(cfg));
  p.run("m.xml");
  ExprScalarFunction exact("1 + 2*x");
  CHECK_NEAR(p.calc_l2_error(exact), 0.0, 1e-6);
}

// ---- Eikonal ----------------------------------------------------------------

namespace
{

// A cube of side 10, 5 elements per direction, fibres along x
BoxSpec eikonal_box()
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 5;
  b.Lx = b.Ly = b.Lz = 10.0;
  return b;
}

toml::table eikonal_cfg(const std::string & type, double vf, double vs, double vn,
                        const std::string & roots)
{
  std::ostringstream os;
  os << "[problem]\nmesh = \"m.xml\"\n"
     << "[physical]\nconductivity_type = \"" << type << "\"\n"
     << "vel_f = " << vf << "\nvel_s = " << vs << "\nvel_n = " << vn << "\n"
     << "[activation]\nsource = \"inline\"\n" << roots
     << "[output]\nfilename = \"eik\"\n";
  return toml::parse(os.str());
}

const char * kCornerRoot = "[[activation.root]]\nid = 0\ntime = 0.0\n";

} // namespace

TEST_CASE("Eikonal: along a mesh edge, activation time is distance over velocity")
{
  const BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_ISOTROPIC", 2.0, 2.0, 2.0, kCornerRoot));
  e.solve();
  const arma::vec & lat = e.get_lat();
  REQUIRE(int(lat.n_elem) == box_num_nodes(b));

  CHECK_NEAR(lat(0), 0.0, 0.0);
  for (int i = 1; i <= 5; i++)
    CHECK_NEAR(lat(box_node(b, i, 0, 0)), i * 2.0 / 2.0, 1e-12);
  CHECK_NEAR(lat(box_node(b, 0, 5, 0)), 5.0, 1e-12);
  CHECK_NEAR(lat(box_node(b, 0, 0, 5)), 5.0, 1e-12);

  // The Kuhn diagonals of this mesh start at the root: exact along the body diagonal
  CHECK_NEAR(lat(box_node(b, 5, 5, 5)), std::sqrt(3.0) * 10.0 / 2.0, 1e-12);
}

TEST_CASE("Eikonal: isotropic activation has the symmetries of the problem")
{
  const BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_ISOTROPIC", 1.0, 1.0, 1.0, kCornerRoot));
  e.solve();
  const arma::vec & lat = e.get_lat();
  for (int i = 0; i <= 5; i++)
    for (int j = 0; j <= 5; j++)
    {
      CHECK_NEAR(lat(box_node(b, i, j, 0)), lat(box_node(b, j, i, 0)), 1e-12);
      CHECK_NEAR(lat(box_node(b, i, 0, j)), lat(box_node(b, j, 0, i)), 1e-12);
    }
  CHECK(lat.min() >= 0.0);
  CHECK(lat.is_finite());
}

TEST_CASE("Eikonal: activation is faster along the fibres")
{
  BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_ORTHOTROPIC", 3.0, 1.0, 0.5, kCornerRoot));
  e.solve();
  const arma::vec & lat = e.get_lat();
  CHECK_NEAR(lat(box_node(b, 5, 0, 0)), 10.0 / 3.0, 1e-12);
  CHECK_NEAR(lat(box_node(b, 0, 5, 0)), 10.0 / 1.0, 1e-12);
  CHECK_NEAR(lat(box_node(b, 0, 0, 5)), 10.0 / 0.5, 1e-12);
}

TEST_CASE("Eikonal: rotating the fibres rotates the fast direction")
{
  BoxSpec b = eikonal_box();
  b.fiber = arma::vec3({0.0, 1.0, 0.0});
  b.sheet = arma::vec3({1.0, 0.0, 0.0});
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_ORTHOTROPIC", 3.0, 1.0, 0.5, kCornerRoot));
  e.solve();
  const arma::vec & lat = e.get_lat();
  CHECK_NEAR(lat(box_node(b, 0, 5, 0)), 10.0 / 3.0, 1e-12);
  CHECK_NEAR(lat(box_node(b, 5, 0, 0)), 10.0 / 1.0, 1e-12);
}

TEST_CASE("Eikonal: transversely isotropic reads the sheet velocity for both cross directions")
{
  BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_TRANSVERSE", 3.0, 1.0, 99.0, kCornerRoot));
  e.solve();
  CHECK_NEAR(e.get_lat()(box_node(b, 0, 0, 5)), 10.0, 1e-12);
}

TEST_CASE("Eikonal: the root times are honoured, and the earliest wavefront wins")
{
  const BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  const int far_corner = box_node(b, 5, 0, 0);
  std::ostringstream roots;
  roots << "[[activation.root]]\nid = 0\ntime = 1.0\n"
        << "[[activation.root]]\nid = " << far_corner << "\ntime = 2.0\n";
  Eikonal e;
  e.setup(eikonal_cfg("M_ISOTROPIC", 1.0, 1.0, 1.0, roots.str()));
  e.solve();
  const arma::vec & lat = e.get_lat();
  CHECK_NEAR(lat(0), 1.0, 0.0);
  CHECK_NEAR(lat(far_corner), 2.0, 0.0);
  CHECK_EQ(e.get_root_nodes().size(), std::size_t(2));

  // Between the roots, the two fronts meet at 1 + x = 2 + (10 - x), x = 5.5
  CHECK_NEAR(lat(box_node(b, 2, 0, 0)), 5.0, 1e-12);   // x = 4, from root 0
  CHECK_NEAR(lat(box_node(b, 3, 0, 0)), 6.0, 1e-12);   // x = 6, from the other root
}

TEST_CASE("Eikonal: the output files hold the result")
{
  const BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_ISOTROPIC", 1.0, 1.0, 1.0, kCornerRoot));
  e.solve();
  CHECK(file_exists("eik_lat.h5"));
  CHECK(file_exists("eik_lat.xmf"));
  const std::string text = read_text("eik.txt");
  CHECK(text.find("<node id=\"0\" lat=\"0\"") != std::string::npos);
}

TEST_CASE("Eikonal: a config with no mesh is rejected, and no roots skips the solve")
{
  toml::table empty;
  Eikonal bad;
  CHECK_THROWS(bad.setup(empty));

  const BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Eikonal e;
  e.setup(eikonal_cfg("M_ISOTROPIC", 1.0, 1.0, 1.0, ""));
  e.solve();
  CHECK_EQ(e.get_lat().n_elem, arma::uword(0));
}

TEST_CASE("Eikonal: it can use a mesh owned by the caller, with roots set directly")
{
  // This is the path of the coupled solvers
  const BoxSpec b = eikonal_box();
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  Eikonal e;
  e.set_mesh(&m);
  e.set_conductivity(M_ISOTROPIC);
  e.set_velocities(2.0, 2.0, 2.0);
  e.set_root_nodes({0}, {0.5});
  e.solve();
  CHECK_NEAR(e.get_lat()(box_node(b, 5, 0, 0)), 0.5 + 5.0, 1e-12);
}
