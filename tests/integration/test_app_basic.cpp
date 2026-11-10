/*
 * Integration tests of the simple solvers, each one run as the real executable on a
 * small problem with a known answer:
 *
 *  - poisson: manufactured solution from a TOML file, with the L2 error it prints, and
 *    the Laplace variant selected in the same file.
 *  - l2projection: runs and writes the projected function.
 *  - elasticity: uniaxial tension of a bar under a traction, against Hooke's law.
 *  - eikonal: the cube and cable cases in tests/data/eikonal, and a generated cube with a
 *    known activation time; the output files must hold what the solver reports.
 *  - Command line: each solver stops with a usage message when called with no input.
 *
 * The executables take the PETSc direct solver through -ksp_type/-pc_type, the same
 * options a user would pass, so that this does not need a GPU.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <filesystem>

using namespace cardiax_test;

namespace
{
const char * kPetsc = "-linsolver petsc -ksp_type preonly -pc_type lu";
}

TEST_CASE("poisson: reports a small L2 error for a problem with a known solution")
{
  BoxSpec b;
  b.nx = 16; b.ny = 2; b.nz = 2;
  std::filesystem::create_directories("poisson");
  write_tet_box("poisson/m.xml", b);
  write_text("poisson/c.toml",
    "[problem]\ntype = \"poisson\"\nmesh = \"m.xml\"\n"
    "[physical]\nsource = \"2\"\n"
    "dirichlet = [ { marker = 1, value = 0.0 }, { marker = 2, value = 0.0 } ]\n"
    "neumann = [ { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 }, { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n"
    "[output]\nfilename = \"out\"\n"
    "[validation]\nexact = \"x*(1-x)\"\n");

  const AppResult r = run_app("poisson", std::string("-f c.toml ") + kPetsc, "poisson");
  CHECK_EQ(r.exit_code, 0);
  const double err = number_after(r.log, "L2 error norm:");
  CHECK_MSG(err < 2e-3, err);
  const arma::mat u = h5_read("poisson/out.h5", "/point_data/solution");
  CHECK_NEAR(u(box_node(b, 8, 1, 1), 0), 0.25, 2e-3);
}

TEST_CASE("poisson: without an exact solution it says so and still solves")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  std::filesystem::create_directories("poisson2");
  write_tet_box("poisson2/m.xml", b);
  write_text("poisson2/c.toml",
    "[problem]\nmesh = \"m.xml\"\n"
    "[physical]\ndirichlet = [ { marker = 1, value = 0.0 } ]\n"
    "neumann = [ { marker = 2, value = 0.0 }, { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 }, { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n");
  const AppResult r = run_app("poisson", std::string("-f c.toml ") + kPetsc, "poisson2");
  CHECK_EQ(r.exit_code, 0);
  CHECK(r.log.find("skipping L2 error") != std::string::npos);
  CHECK(std::filesystem::exists("poisson2/output.h5"));
}

TEST_CASE("poisson: a malformed TOML file fails with an error")
{
  std::filesystem::create_directories("poisson3");
  write_text("poisson3/c.toml", "[problem\nmesh = ");
  const AppResult r = run_app("poisson", "-f c.toml", "poisson3");
  CHECK(r.exit_code != 0);
  CHECK(r.log.find("ERROR parsing") != std::string::npos);
}

TEST_CASE_KNOWN_BUG("poisson: a boundary marker with no condition is an error",
                    "src/pdes/poisson.cpp, coeff_robin_mat/coeff_robin_vec: prints a message and "
                    "calls exit(0), so a misconfigured run looks successful to a script")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  std::filesystem::create_directories("poisson4");
  write_tet_box("poisson4/m.xml", b);
  write_text("poisson4/c.toml",
    "[problem]\nmesh = \"m.xml\"\n[physical]\ndirichlet = [ { marker = 1, value = 0.0 } ]\n");
  const AppResult r = run_app("poisson", std::string("-f c.toml ") + kPetsc, "poisson4");
  CHECK(r.exit_code != 0);
}

TEST_CASE("poisson: the Laplace problem is selected in the file")
{
  BoxSpec b;
  b.nx = 4; b.ny = b.nz = 2;
  std::filesystem::create_directories("laplace");
  write_tet_box("laplace/m.xml", b);
  write_text("laplace/c.toml",
    "[problem]\ntype = \"laplace\"\nmesh = \"m.xml\"\n"
    "[physical]\n"
    "dirichlet = [ { marker = 1, value = 1.0 }, { marker = 2, value = 3.0 } ]\n"
    "neumann = [ { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 }, { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n"
    "[validation]\nexact = \"1 + 2*x\"\n");
  const AppResult r = run_app("poisson", std::string("-f c.toml ") + kPetsc, "laplace");
  CHECK_EQ(r.exit_code, 0);
  const double err = number_after(r.log, "L2 error norm:");
  CHECK_MSG(err < 1e-5, err);
}

TEST_CASE("l2projection: writes the projection of its source function")
{
  BoxSpec b;
  b.nx = b.ny = 16;
  std::filesystem::create_directories("l2");
  write_tri_square("l2/m.xml", b);
  const AppResult r = run_app("l2projection", std::string("-f m.xml -o proj ") + kPetsc, "l2");
  CHECK_EQ(r.exit_code, 0);
  REQUIRE(std::filesystem::exists("l2/proj.h5"));
  const arma::mat u = h5_read("l2/proj.h5", "/point_data/solution");
  // sin(pi x) sin(pi y) at the centre of the unit square
  CHECK_NEAR(u(8 + 8 * 17, 0), 1.0, 0.02);
}

TEST_CASE_KNOWN_BUG("elasticity: uniaxial tension of a bar follows Hooke's law",
                    "src/pdes/elasticity.cpp, under #define DOBRAR: assemble_system reads K_l with "
                    "MatGetRowIJ without assembling it (the run aborts with a garbage column index), "
                    "duplicates the system, and numbers fixed dofs as dir*np + node while the space "
                    "uses ndim*node + dir")
{
  const double E = 1000.0, nu = 0.3, traction = 10.0, L = 4.0;
  BoxSpec b;
  b.nx = 4; b.ny = b.nz = 2;
  b.Lx = L;

  // Symmetry conditions: x = 0 slides in the plane, and the lateral motion is free
  std::ostringstream fixed;
  for (int k = 0; k <= b.nz; k++)
    for (int j = 0; j <= b.ny; j++)
      fixed << "<node id=\"" << box_node(b, 0, j, k) << "\" direction=\"0\" value=\"0\"/>\n";
  fixed << "<node id=\"" << box_node(b, 0, 0, 0) << "\" direction=\"1\" value=\"0\"/>\n"
        << "<node id=\"" << box_node(b, 0, 0, 0) << "\" direction=\"2\" value=\"0\"/>\n"
        << "<node id=\"" << box_node(b, 0, 2, 0) << "\" direction=\"2\" value=\"0\"/>\n";

  std::ostringstream prefix;
  prefix << "<elasticity type=\"THREE_DIM\"><parameters><young>" << E << "</young>"
         << "<poisson>" << nu << "</poisson></parameters>"
         << "<neumann><node id=\"1\" marker=\"2\" t0=\"" << traction << "\" t1=\"0\" t2=\"0\"/></neumann>"
         << "<prescribed_displacement>" << fixed.str() << "</prescribed_displacement></elasticity>\n";
  b.prefix_xml = prefix.str();

  std::filesystem::create_directories("elas");
  write_hex_box("elas/m.xml", b);
  const AppResult r = run_app("elasticity", std::string("-f m.xml -o res ") + kPetsc, "elas");
  CHECK_EQ(r.exit_code, 0);
  REQUIRE(std::filesystem::exists("elas/res.h5"));

  const arma::mat u = h5_read("elas/res.h5", "/point_data/displacements");
  const int tip = box_node(b, 4, 1, 1);
  CHECK_REL(u(0, 3 * tip + 0), traction * L / E, 1e-6);
  const int side = box_node(b, 2, 2, 1);
  CHECK_REL(u(0, 3 * side + 1) - u(0, 3 * box_node(b, 2, 0, 1) + 1),
            -nu * traction / E * 1.0, 1e-6);
}

TEST_CASE("eikonal: generated cube, isotropic")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 4;
  b.Lx = b.Ly = b.Lz = 8.0;
  std::filesystem::create_directories("eik");
  write_tet_box("eik/m.xml", b);
  write_text("eik/c.toml",
    "[problem]\nmesh = \"m.xml\"\n"
    "[physical]\nconductivity_type = \"M_ISOTROPIC\"\nvel_f = 2.0\n"
    "[activation]\nsource = \"inline\"\n[[activation.root]]\nid = 0\ntime = 1.0\n"
    "[output]\nfilename = \"res\"\n");
  const AppResult r = run_app("eikonal", "-f c.toml", "eik");
  CHECK_EQ(r.exit_code, 0);
  CHECK_NEAR(number_after(r.log, "Earliest activation (min):"), 1.0, 1e-9);
  CHECK(std::filesystem::exists("eik/res_lat.h5"));
  CHECK(std::filesystem::exists("eik/res_lat.xmf"));
  const std::string list = read_text("eik/res.txt");
  CHECK(list.find("<node id=\"4\" lat=\"5\"") != std::string::npos);  // 1 + 8/2
}

TEST_CASE("eikonal: repository cases run and keep their activation range")
{
  const struct { const char * toml; const char * mesh; const char * key; } cases[] = {
    {"tests/data/eikonal/input_cube.toml",  "tests/data/eikonal/cube.xml",            "eikonal.cube.latmax"},
    {"tests/data/eikonal/input_cable.toml", "tests/data/eikonal/cable_10x_fiber.xml", "eikonal.cable.latmax"},
    {"tests/data/eikonal/input_donut.toml", "tests/data/eikonal/donut.xml",           "eikonal.donut.latmax"},
  };
  for (const auto & c : cases)
  {
    const std::string dir = std::string("repo_") + c.key;
    const std::string toml = copy_input(c.toml, dir);
    copy_input(c.mesh, dir);
    const AppResult r = run_app("eikonal", "-f " + toml, dir);
    CHECK_MSG(r.exit_code == 0, c.toml);
    const double lat = number_after(r.log, "Latest activation (max):");
    CHECK_MSG(lat > 0.0, c.toml << " " << lat);
    CHECK_GOLDEN(c.key, lat, 1e-9);
  }
}

TEST_CASE("Command line: each solver reports its usage when given no input")
{
  for (const char * app : {"poisson", "l2projection", "elasticity", "nonlinearelas",
                           "monodomain", "bidomain", "cardiac_cycle", "electromech", "eikonal"})
  {
    const AppResult r = run_app(app, "", std::string("usage_") + app);
    CHECK_MSG(r.log.find("sage") != std::string::npos, app << ": " << r.log.substr(0, 200));
  }
}
