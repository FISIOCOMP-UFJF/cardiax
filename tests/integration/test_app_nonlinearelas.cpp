/*
 * Integration tests of the nonlinearelas executable: quasi-static mechanics with Newton.
 *
 * What is covered and why:
 *  - Uniaxial stretch of a free cube with prescribed displacement: homogeneous, and against
 *    the closed form of the compressible neo-Hookean law.
 *  - The nearly incompressible limit: the tissue volume is preserved.
 *  - Updated and Total Lagrangian agree on the same problem.
 *  - Convergence of Newton: the number of iterations is bounded. A broken tangent would
 *    still reach the solution through line search, but with many more iterations.
 *  - Following pressure and passive cases in tests/data/passive; active cases in
 *    tests/data/active; the number of increments must not change the converged state.
 *  - Land et al. 2015 benchmark problem 2 (examples/prob2_12x27x2_k300.xml), a pressurised
 *    thin ellipsoid: the deformed state against the regression baseline.
 *  - PV history and volume files, checkpoint written with -save_state, and restart from
 *    the middle of the load ramp.
 *  - Springs (Robin condition) on the boundary.
 *
 * nonlinearelas reads a TOML file (-f): each test writes it with write_nonlinearelas_toml,
 * and the mesh XML it names keeps only the boundary-condition topology.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <filesystem>

using namespace cardiax_test;
namespace fs = std::filesystem;

namespace
{

const char * kArgs = " -num_threads 2 -linsolver petsc -ksp_type preonly -pc_type lu";

// The mesh XML keeps only the topology of the boundary conditions; the values go in the
// TOML (see write_input)
std::string elasticity_block(const std::string & body)
{
  return "<elasticity type=\"THREE_DIM\">" + body + "</elasticity>\n";
}

// Writes <dir>/<base>.toml for the mesh <dir>/<mesh> and returns "<base>.toml"
std::string write_input(const std::string & dir, const std::string & material,
                        const std::string & coefs, int ninc,
                        const std::vector<std::pair<int, double>> & pressure = {},
                        const std::string & mesh = "m.xml")
{
  NonlinearElasInput in;
  in.mesh = mesh;
  in.material = material;
  in.coefs = coefs;
  in.ninc = ninc;
  in.pressure = pressure;
  const std::string toml = mesh.substr(0, mesh.rfind(".xml")) + ".toml";
  write_nonlinearelas_toml(dir + "/" + toml, in);
  return toml;
}

const char * kGuccioneStiff = "2.0, 8.0, 2.0, 4.0, 0.0, 1000.0";
const char * kGuccione      = "2.0, 8.0, 2.0, 4.0, 0.0, 200.0";

// Symmetric stretch of a unit cube along x by a total displacement d: x = 0 fixed in x,
// x = L moved in x, three nodes hold the rigid modes.
std::string uniaxial_bc(const BoxSpec & b, double d)
{
  std::ostringstream os;
  os << "<prescribed_displacement>";
  for (int k = 0; k <= b.nz; k++)
    for (int j = 0; j <= b.ny; j++)
    {
      os << "<node id=\"" << box_node(b, 0, j, k) << "\" direction=\"0\" value=\"0\"/>";
      os << "<node id=\"" << box_node(b, b.nx, j, k) << "\" direction=\"0\" value=\"" << d << "\"/>";
    }
  os << "<node id=\"" << box_node(b, 0, 0, 0) << "\" direction=\"1\" value=\"0\"/>"
     << "<node id=\"" << box_node(b, 0, 0, 0) << "\" direction=\"2\" value=\"0\"/>"
     << "<node id=\"" << box_node(b, 0, b.ny, 0) << "\" direction=\"2\" value=\"0\"/>"
     << "</prescribed_displacement>";
  return os.str();
}

// Runs and returns the final displacement field, with rows = nodes, cols = xyz
arma::mat final_displacement(const std::string & dir, const std::string & base)
{
  const arma::mat u = h5_read(dir + "/" + base + "_output_nl.h5", "/point_data/displacements");
  const int np = u.n_cols / 3;
  // The file has one row per increment; pick the last row with a nonzero entry
  arma::uword row = 0;
  for (arma::uword i = 0; i < u.n_rows; i++)
    if (arma::abs(u.row(i)).max() > 0.0) row = i;
  arma::mat out(np, 3);
  for (int i = 0; i < np; i++)
    for (int d = 0; d < 3; d++)
      out(i, d) = u(row, 3 * i + d);
  return out;
}

int newton_iterations(const std::string & log)
{
  int n = 0;
  const std::string key = "Newton-LS step";
  for (std::size_t p = log.find(key); p != std::string::npos; p = log.find(key, p + key.size()))
    n++;
  return n;
}

} // namespace

TEST_CASE("nonlinearelas: a homogeneous uniaxial stretch stays homogeneous, in both formulations")
{
  const double d = 0.2;
  arma::mat results[2];
  const char * forms[2] = {"ul", "tl"};
  for (int f = 0; f < 2; f++)
  {
    const std::string dir = std::string("uni_") + forms[f];
    fs::create_directories(dir);
    BoxSpec b;
    b.nx = b.ny = b.nz = 2;
    b.prefix_xml = elasticity_block(uniaxial_bc(b, d));
    write_hex_box(dir + "/m.xml", b);
    const std::string toml = write_input(dir, "Guccione", kGuccioneStiff, 4);

    const AppResult r = run_app("nonlinearelas", "-f " + toml + " -s " + forms[f] + kArgs, dir);
    REQUIRE(r.exit_code == 0);
    results[f] = final_displacement(dir, "m");

    // Linear displacement in x, and the same lateral contraction at both ends
    CHECK_NEAR(results[f](box_node(b, 2, 1, 1), 0), d, 1e-6);
    CHECK_NEAR(results[f](box_node(b, 1, 1, 1), 0), d / 2.0, 1e-6);
    const double lat_mid = results[f](box_node(b, 1, 2, 1), 1) - results[f](box_node(b, 1, 0, 1), 1);
    const double lat_end = results[f](box_node(b, 2, 2, 1), 1) - results[f](box_node(b, 2, 0, 1), 1);
    CHECK_NEAR(lat_mid, lat_end, 1e-6);
    CHECK_MSG(lat_mid < 0.0, "a stretched bar contracts laterally: " << lat_mid);
  }
  // Same law, but UL uses mean dilatation for the volumetric part and TL does not
  CHECK_NEAR(arma::abs(results[0] - results[1]).max(), 0.0, 5e-4);
}

// Was a known bug until allocate_Ta stopped reading an empty vector (5668b24)
TEST_CASE("nonlinearelas: zero lateral stress gives the closed form of the compressible neo-Hookean")
{
  // Bonet & Wood: psi = mu/2 (I1 - 3) - mu ln J + lam/2 (ln J)^2. With free lateral faces
  // the transverse stretch t solves mu (t^2 - 1) + lam ln(L t^2) = 0.
  const double mu = 1.0, lam = 2.0, L = 1.2;
  double t = 1.0;
  for (int i = 0; i < 100; i++)
  {
    const double f = mu * (t * t - 1.0) + lam * std::log(L * t * t);
    const double df = 2.0 * mu * t + 2.0 * lam / t;
    t -= f / df;
  }

  fs::create_directories("closed");
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block(uniaxial_bc(b, L - 1.0));
  write_hex_box("closed/m.xml", b);
  write_input("closed", "CompressibleNeoHookean", "1.0, 2.0", 4);
  const AppResult r = run_app("nonlinearelas", std::string("-f m.toml -s tl") + kArgs, "closed");
  REQUIRE(r.exit_code == 0);
  const arma::mat u = final_displacement("closed", "m");
  const double width = 1.0 + u(box_node(b, 1, 2, 1), 1) - u(box_node(b, 1, 0, 1), 1);
  CHECK_NEAR(width, t, 1e-4);
}

TEST_CASE("nonlinearelas: a nearly incompressible stretch preserves the volume")
{
  fs::create_directories("incomp");
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block(uniaxial_bc(b, 0.1));
  write_hex_box("incomp/m.xml", b);
  write_input("incomp", "Guccione", kGuccioneStiff, 4);
  const AppResult r = run_app("nonlinearelas", std::string("-f m.toml") + kArgs, "incomp");
  REQUIRE(r.exit_code == 0);
  // <base>_volumes.dat holds V_LV, V_RV and the tissue volume at the end, one per line
  const double v = read_table("incomp/m_volumes.dat")(2, 0);
  CHECK_NEAR(v, 1.0, 5e-3);
}

TEST_CASE("nonlinearelas: Newton converges in few iterations per increment")
{
  fs::create_directories("newton");
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block(uniaxial_bc(b, 0.1));
  write_tet_box("newton/m.xml", b);
  write_input("newton", "Guccione", kGuccione, 5);
  const AppResult r = run_app("nonlinearelas", std::string("-f m.toml") + kArgs, "newton");
  REQUIRE(r.exit_code == 0);
  const int its = newton_iterations(r.log);
  CHECK_MSG(its > 0 && its <= 5 * 8, its << " iterations in 5 increments");
}

TEST_CASE("nonlinearelas: the passive cases of tests/data/passive run to the same volume per material")
{
  // The values are those of the <parameters> and <pressure> blocks of each XML
  const struct {
    const char * file; const char * key; const char * material; const char * coefs;
    std::vector<std::pair<int, double>> pressure;
  } cases[] = {
    {"tests/data/passive/teste_unitario_guccione.xml",  "passive.guccione",
     "Guccione", "0.650, 6.62, 3.65, 2.65, 100, 200.0", {{10, 4.0}}},
    {"tests/data/passive/teste_unitario_holzapfel.xml", "passive.holzapfel",
     "HolzapfelOgden", "0.15, 0.11685, 0.372, 0.410, 6.0, 11.83425, 5.16, 11.3, 0, 333.0",
     {{50, 4.0}, {10, 0.0}}},
  };
  for (const auto & c : cases)
  {
    const std::string name = copy_input(c.file, c.key);
    const std::string toml = write_input(c.key, c.material, c.coefs, 10, c.pressure, name);
    const AppResult r = run_app("nonlinearelas", "-f " + toml + kArgs, c.key);
    CHECK_MSG(r.exit_code == 0, c.file);
    const std::string base = name.substr(0, name.rfind(".xml"));
    const arma::mat u = final_displacement(c.key, base);
    CHECK_GOLDEN(std::string(c.key) + ".unorm", arma::norm(arma::vectorise(u)), 1e-6);
  }
}

TEST_CASE("nonlinearelas: active stress cases, and the increments do not change the converged state")
{
  double norm_of[2] = {0.0, 0.0};
  const char * dirs[2] = {"ninc1", "ninc10"};
  const int ninc[2] = {1, 25};   // the <ninc> of each XML (ninc10 holds 25)
  for (int i = 0; i < 2; i++)
  {
    const std::string dir = std::string("active_") + dirs[i];
    const std::string name = copy_input(std::string("tests/data/active/") + dirs[i] +
                                        "/teste_unitario_guccione.xml", dir);
    const std::string toml = write_input(dir, "Guccione", "0.650, 6.62, 3.65, 2.65, 100, 200.0",
                                         ninc[i], {{10, 0.0}}, name);
    const AppResult r = run_app("nonlinearelas", "-f " + toml + kArgs, dir);
    REQUIRE(r.exit_code == 0);
    const arma::mat u = final_displacement(dir, name.substr(0, name.rfind(".xml")));
    norm_of[i] = arma::norm(arma::vectorise(u));
    CHECK(u.is_finite());
    CHECK(norm_of[i] > 0.0);
  }
  CHECK_REL(norm_of[0], norm_of[1], 1e-4);
  CHECK_GOLDEN("active.guccione.unorm", norm_of[1], 1e-6);
}

TEST_CASE("nonlinearelas: Land et al. benchmark 2, pressurised ellipsoid")
{
  const std::string name = copy_input("examples/prob2_12x27x2_k300.xml", "prob2");
  const std::string toml = write_input("prob2", "Guccione", "10, 1, 1, 1, 0, 300", 100,
                                       {{11, 10.0}}, name);
  const AppResult r = run_app("nonlinearelas", "-f " + toml + kArgs, "prob2");
  REQUIRE(r.exit_code == 0);
  const arma::mat u = final_displacement("prob2", name.substr(0, name.rfind(".xml")));
  CHECK(u.is_finite());

  const double apex_like = arma::abs(u.col(0)).max();
  CHECK_GOLDEN("prob2.umax", apex_like, 1e-5);
  CHECK_GOLDEN("prob2.unorm", arma::norm(arma::vectorise(u)), 1e-5);
}

TEST_CASE("nonlinearelas: PV history and volume files of a pressurised box")
{
  fs::create_directories("pv");
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  std::ostringstream fixed;
  fixed << "<prescribed_displacement>";
  for (int j = 0; j <= b.ny; j++)
    for (int i = 0; i <= b.nx; i++)
      for (int d = 0; d < 3; d++)
        fixed << "<node id=\"" << box_node(b, i, j, b.nz) << "\" direction=\"" << d << "\" value=\"0\"/>";
  fixed << "</prescribed_displacement>";
  b.prefix_xml = elasticity_block(fixed.str());
  write_tet_box("pv/m.xml", b);
  write_input("pv", "Guccione", kGuccione, 5,
              {{30, 1.0}, {2, 0.0}, {3, 0.0}, {4, 0.0}, {5, 0.0}, {6, 0.0}});

  // Turn x = 0 into the 'LV' surface, marker 30
  std::string xml = read_text("pv/m.xml");
  const std::string from = "marker=\"1\" ", to = "marker=\"30\" ";
  for (std::size_t p = xml.find(from); p != std::string::npos; p = xml.find(from, p + to.size()))
    xml.replace(p, from.size(), to);
  write_text("pv/m.xml", xml);

  const AppResult r = run_app("nonlinearelas", std::string("-f m.toml") + kArgs, "pv");
  REQUIRE(r.exit_code == 0);

  const arma::mat h = read_table("pv/m_pv_history.dat");
  REQUIRE(h.n_rows == 6);
  CHECK_NEAR(h(0, 1), 0.0, 0.0);
  CHECK_NEAR(h(5, 1), 1.0, 1e-9);

  const arma::mat v = read_table("pv/m_volumes.dat");
  CHECK_EQ(v.n_rows, arma::uword(3));
}

TEST_CASE("nonlinearelas: a spring boundary restrains the motion")
{
  // Every boundary marker needs a load or a spring: the loaded end, the fixed end, and four
  // lateral faces that are either free (zero pressure) or held by normal springs.
  // `springs` is an XML <spring> block; `lateral` the pressures of the lateral faces.
  using Pressures = std::vector<std::pair<int, double>>;
  auto tip = [](const std::string & dir, const std::string & springs, const Pressures & lateral) {
    fs::create_directories(dir);
    BoxSpec b;
    b.nx = b.ny = b.nz = 1;
    std::ostringstream fixed;
    fixed << "<prescribed_displacement>";
    for (int k = 0; k <= 1; k++)
      for (int j = 0; j <= 1; j++)
        for (int d = 0; d < 3; d++)
          fixed << "<node id=\"" << box_node(b, 0, j, k) << "\" direction=\"" << d << "\" value=\"0\"/>";
    fixed << "</prescribed_displacement>";
    b.prefix_xml = elasticity_block(springs + fixed.str());
    write_tet_box(dir + "/m.xml", b);
    Pressures pressure = {{2, -0.1}, {1, 0.0}};
    pressure.insert(pressure.end(), lateral.begin(), lateral.end());
    write_input(dir, "Guccione", kGuccione, 2, pressure);
    const AppResult r = run_app("nonlinearelas", std::string("-f m.toml") + kArgs, dir);
    REQUIRE(r.exit_code == 0);
    const arma::mat u = final_displacement(dir, "m");
    return u(box_node(b, 1, 0, 0), 0);
  };
  const std::string springs =
    "<spring><node id=\"1\" marker=\"3\" value=\"10.0\"/><node id=\"2\" marker=\"4\" value=\"10.0\"/>"
    "<node id=\"3\" marker=\"5\" value=\"10.0\"/><node id=\"4\" marker=\"6\" value=\"10.0\"/></spring>";
  const double free_tip = tip("spring0", "", {{3, 0.0}, {4, 0.0}, {5, 0.0}, {6, 0.0}});
  const double held_tip = tip("spring1", springs, {});
  CHECK_MSG(std::fabs(free_tip) > 1e-4, free_tip);
  CHECK_MSG(std::fabs(held_tip) < std::fabs(free_tip), free_tip << " vs " << held_tip);
}

TEST_CASE("nonlinearelas: an unknown formulation or a missing mesh is an error")
{
  fs::create_directories("bad");
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  b.prefix_xml = elasticity_block("");
  write_tet_box("bad/m.xml", b);
  write_input("bad", "Guccione", kGuccione, 1);
  write_input("bad", "Guccione", kGuccione, 1, {}, "nope.xml");   // a TOML whose mesh is missing
  CHECK(run_app("nonlinearelas", std::string("-f m.toml -s xx") + kArgs, "bad").exit_code != 0);
  CHECK(run_app("nonlinearelas", std::string("-f nope.toml") + kArgs, "bad").exit_code != 0);
  CHECK(run_app("nonlinearelas", std::string("-f missing.toml") + kArgs, "bad").exit_code != 0);
}

TEST_CASE("nonlinearelas: restart from a checkpoint in the middle of the load ramp")
{
  fs::create_directories("chk");
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  b.prefix_xml = elasticity_block(uniaxial_bc(b, 0.2));
  write_hex_box("chk/m.xml", b);
  write_input("chk", "Guccione", kGuccione, 4);

  const AppResult full = run_app("nonlinearelas", std::string("-f m.toml -save_state 1") + kArgs, "chk");
  REQUIRE(full.exit_code == 0);
  const arma::mat u_full = final_displacement("chk", "m");

  // A checkpoint per increment; restart from the first one found
  std::vector<std::string> checkpoints;
  for (const auto & e : fs::directory_iterator("chk"))
    if (e.path().extension() == ".h5" && e.path().filename().string() != "m_output_nl.h5")
      checkpoints.push_back(e.path().filename().string());
  REQUIRE(!checkpoints.empty());
  std::sort(checkpoints.begin(), checkpoints.end());

  const AppResult again =
      run_app("nonlinearelas", "-f m.toml -restore " + checkpoints.front() + kArgs, "chk");
  REQUIRE(again.exit_code == 0);
  const arma::mat u_rest = final_displacement("chk", "m");
  CHECK_NEAR(arma::abs(u_full - u_rest).max(), 0.0, 1e-5);
}
