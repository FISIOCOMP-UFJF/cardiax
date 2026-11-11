/*
 * Integration tests of the electrophysiology executables, on small tissue meshes.
 *
 * What is covered and why:
 *  - monodomain, TOML input: a stimulus at one end of a slab gives an action potential
 *    that travels to the other end, at a conduction velocity in the physiological range
 *    and matching the regression baseline. Faster along the fibres than across them.
 *  - Crank-Nicolson and implicit Euler (the numerics.theta option) give close velocities.
 *  - No stimulus means no activation.
 *  - The same answer with one and two OpenMP threads (also for ToRORdLand, whose
 *    result depends on the cell type, with no cell types in the mesh).
 *  - Different ionic models run through the same path.
 *  - Checkpoint: -save_state writes a file, -restore continues the run and arrives at the
 *    same final state as the uninterrupted run.
 *  - Pseudo-ECG: the output file exists with one column per lead.
 *  - Purkinje network (-ep purkinje) and tissue coupled with a Purkinje network
 *    (-ep coupled_purkinje): the wave crosses the network and, when coupled, reaches
 *    the tissue through the terminal junctions.
 *  - bidomain: the wave propagates like the monodomain one, with a finite extracellular
 *    potential.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <filesystem>

using namespace cardiax_test;
namespace fs = std::filesystem;

namespace
{

const char * kLU = " -linsolver petsc -ksp_type preonly -pc_type lu";

// Slab of 2 x 0.1 x 0.1 cm, in um. Short enough to run in seconds.
BoxSpec slab(int nx = 40)
{
  BoxSpec b;
  b.nx = nx; b.ny = 2; b.nz = 2;
  b.Lx = 20000.0; b.Ly = 1000.0; b.Lz = 1000.0;
  return b;
}

std::string mono_toml(const std::string & cell, double t_end, double dt, double theta,
                      const std::string & cond = "isotropic",
                      const std::string & stim_region = "max = [1000.0, 1000.0, 1000.0]",
                      const std::string & extra = "")
{
  std::ostringstream os;
  os << "[problem]\nmesh = \"m.xml\"\nmodel = \"monodomain\"\n"
     << "[time]\ndt = " << dt << "\nt_end = " << t_end << "\ndt_output = 1.0\n"
     << "[cell_model]\nname = \"" << cell << "\"\n"
     << "[numerics]\node_solver = \"ExplicitEuler\"\ntheta = " << theta << "\n"
     << "[physical]\nconductivity_type = \"" << cond << "\"\n"
     << "sigma_f = 0.0005\nsigma_s = 0.0001\nsigma_n = 0.0001\nsurface_to_volume = 0.14\n"
     << "[stimuli]\ntissue = [ { start = 0.0, duration = 2.0, value = -52.0, "
     << "min = [0.0, 0.0, 0.0], " << stim_region << " } ]\n"
     << extra;
  return os.str();
}

// First output step at which each node crosses the threshold (NaN if never)
arma::vec activation_step(const arma::mat & vm, double threshold = -20.0)
{
  arma::vec at(vm.n_cols);
  at.fill(arma::datum::nan);
  for (arma::uword n = 0; n < vm.n_cols; n++)
    for (arma::uword s = 0; s < vm.n_rows; s++)
      if (vm(s, n) > threshold) { at(n) = static_cast<double>(s); break; }
  return at;
}

arma::mat run_monodomain(const std::string & dir, const BoxSpec & b, const std::string & toml,
                         const std::string & args = "", int threads = 2)
{
  fs::create_directories(dir);
  write_hex_box(dir + "/m.xml", b);
  write_text(dir + "/c.toml", toml);
  const AppResult r = run_app("monodomain",
      "-f c.toml -num_threads " + std::to_string(threads) + " " + args + kLU, dir);
  CHECK_MSG(r.exit_code == 0, r.log.substr(r.log.size() > 800 ? r.log.size() - 800 : 0));
  return h5_read(dir + "/m_output.h5", "/point_data/vm");
}

// Conduction velocity [um/ms] between 1/4 and 3/4 of the slab
double conduction_velocity(const BoxSpec & b, const arma::mat & vm, double dt_output)
{
  const arma::vec at = activation_step(vm);
  const int i1 = b.nx / 4, i2 = 3 * b.nx / 4;
  const double t1 = at(box_node(b, i1, 1, 1)) * dt_output;
  const double t2 = at(box_node(b, i2, 1, 1)) * dt_output;
  return (i2 - i1) * (b.Lx / b.nx) / (t2 - t1);
}

} // namespace

TEST_CASE("monodomain: the action potential crosses the slab at a physiological velocity")
{
  const BoxSpec b = slab();
  const arma::mat vm = run_monodomain("cv", b, mono_toml("TT2", 40.0, 0.02, 0.5));

  CHECK(vm.is_finite());
  CHECK_MSG(vm.max() > 10.0, "no upstroke: " << vm.max());
  CHECK_MSG(vm.max() < 80.0, "potential out of range: " << vm.max());

  const arma::vec at = activation_step(vm);
  CHECK_MSG(at.is_finite(), "some nodes were not activated");

  // Monotonic front along x
  for (int i = 1; i <= b.nx; i++)
    CHECK(at(box_node(b, i, 1, 1)) >= at(box_node(b, i - 1, 1, 1)));

  const double cv = conduction_velocity(b, vm, 1.0);
  CHECK_MSG(cv > 200.0 && cv < 1500.0, "CV " << cv << " um/ms");
  CHECK_GOLDEN("monodomain.TT2.cv", cv, 1e-6);
  CHECK_GOLDEN("monodomain.TT2.vm_last", vm(vm.n_rows - 1, box_node(b, b.nx, 1, 1)), 1e-6);
}

TEST_CASE("monodomain: no stimulus, no activation")
{
  const BoxSpec b = slab(10);
  std::string toml = mono_toml("TT2", 5.0, 0.02, 0.5);
  toml = toml.substr(0, toml.find("[stimuli]"));
  const arma::mat vm = run_monodomain("rest", b, toml);
  CHECK_MSG(vm.max() < -70.0, vm.max());
}

TEST_CASE("monodomain: faster along the fibres than across them")
{
  auto cv_with_fibre = [](const std::string & dir, const arma::vec3 & f, const arma::vec3 & s) {
    BoxSpec b = slab();
    b.fiber = f;
    b.sheet = s;
    const arma::mat vm = run_monodomain(dir, b, mono_toml("TT2", 40.0, 0.02, 0.5, "orthotropic"));
    return conduction_velocity(b, vm, 1.0);
  };
  const double along = cv_with_fibre("along", arma::vec3({1, 0, 0}), arma::vec3({0, 1, 0}));
  const double across = cv_with_fibre("across", arma::vec3({0, 1, 0}), arma::vec3({0, 0, 1}));
  CHECK_MSG(across > 0.0 && along > 1.5 * across, along << " vs " << across);
}

TEST_CASE_KNOWN_BUG("monodomain: implicit Euler (theta = 1) and Crank-Nicolson agree on the velocity",
                    "src/cardiac/monodomain.cpp, assemble_matrices: the right-hand matrix uses "
                    "M - theta*K where the theta scheme has M - (1 - theta)*K. Only theta = 0.5 is "
                    "right; theta = 1 is Crank-Nicolson with twice the conductivity, which gives "
                    "the velocity times sqrt(2)")
{
  const BoxSpec b = slab();
  auto fine_output = [](std::string toml) {
    toml.replace(toml.find("dt_output = 1.0"), 15, "dt_output = 0.1");
    return toml;
  };
  const double cn = conduction_velocity(b, run_monodomain("cn", b, fine_output(mono_toml("TT2", 40.0, 0.02, 0.5))), 0.1);
  const double ie = conduction_velocity(b, run_monodomain("ie", b, fine_output(mono_toml("TT2", 40.0, 0.02, 1.0))), 0.1);
  CHECK_REL(cn, ie, 0.1);
}

TEST_CASE("monodomain: the result does not depend on the number of threads")
{
  const BoxSpec b = slab(10);
  const std::string toml = mono_toml("TT2", 10.0, 0.02, 0.5);
  const arma::mat a = run_monodomain("t1", b, toml, "", 1);
  const arma::mat c = run_monodomain("t2", b, toml, "", 2);
  CHECK_NEAR(arma::abs(a - c).max(), 0.0, 1e-8);
}

// ToRORdLand depends on the cell type (endo/mid/epi). With no types in the mesh,
// CellModel::type used to be left uninitialized, so each per-thread copy of the cell
// model ran with whatever was in its memory and the result changed with the
// thread -> node mapping.
TEST_CASE("monodomain: ToRORdLand without cell types does not depend on the number of threads")
{
  const BoxSpec b = slab(10);
  const std::string toml = mono_toml("ToRORdLand", 10.0, 0.02, 0.5);
  const arma::mat a = run_monodomain("tor_t1", b, toml, "", 1);
  const arma::mat c = run_monodomain("tor_t4", b, toml, "", 4);
  CHECK(a.is_finite());
  CHECK_NEAR(arma::abs(a - c).max(), 0.0, 1e-8);
}

TEST_CASE("monodomain: the command line overrides the file")
{
  const BoxSpec b = slab(10);
  const arma::mat vm = run_monodomain("cli", b, mono_toml("TT2", 40.0, 0.02, 0.5), "-t 5");
  CHECK_MSG(vm.n_rows <= 6, "output rows " << vm.n_rows);
}

TEST_CASE("monodomain: other ionic models run through the same path")
{
  for (const char * cell : {"TNNP", "LR1"})
  {
    const BoxSpec b = slab(10);
    const std::string dir = std::string("cell_") + cell;
    const arma::mat vm = run_monodomain(dir, b, mono_toml(cell, 5.0, 0.01, 0.5));
    CHECK_MSG(vm.is_finite(), cell);
    CHECK_MSG(vm.max() > 0.0, cell << " peak " << vm.max());
  }
}

TEST_CASE("monodomain: ToRORd-Land stays finite and activates")
{
  const BoxSpec b = slab(10);
  std::string toml = mono_toml("ToRORdLand", 30.0, 0.02, 0.5);
  toml.replace(toml.find("value = -52.0"), 13, "value = -53.0");
  const arma::mat vm = run_monodomain("cell_torord", b, toml);
  CHECK(vm.is_finite());
  CHECK_MSG(vm.max() > 0.0 && vm.max() < 90.0, vm.max());
}

TEST_CASE("monodomain: an unknown model or cell name is an error")
{
  const BoxSpec b = slab(4);
  fs::create_directories("badcell");
  write_hex_box("badcell/m.xml", b);
  write_text("badcell/c.toml", mono_toml("NoSuchCell", 1.0, 0.1, 0.5));
  CHECK(run_app("monodomain", std::string("-f c.toml") + kLU, "badcell").exit_code != 0);

  write_text("badcell/nomesh.toml", "[time]\ndt = 0.1\n");
  CHECK(run_app("monodomain", std::string("-f nomesh.toml") + kLU, "badcell").exit_code != 0);
}

TEST_CASE("monodomain: restarting from a checkpoint reaches the same state")
{
  const BoxSpec b = slab(10);
  const std::string toml = mono_toml("TT2", 10.0, 0.02, 0.5);

  run_monodomain("restart_full", b, toml, "-save_state 5");
  REQUIRE(fs::exists("restart_full/checkpoint_t_5.00ms.h5"));
  REQUIRE(fs::exists("restart_full/checkpoint_t_10.00ms.h5"));

  fs::create_directories("restart_rest");
  fs::copy_file("restart_full/checkpoint_t_5.00ms.h5", "restart_rest/start.h5",
                fs::copy_options::overwrite_existing);
  fs::rename("restart_full/checkpoint_t_10.00ms.h5", "restart_full/final.h5");
  run_monodomain("restart_rest", b, toml, "-save_state 5 -restore start.h5");
  REQUIRE(fs::exists("restart_rest/checkpoint_t_10.00ms.h5"));

  std::string why;
  CHECK_MSG(h5_equal("restart_full/final.h5", "restart_rest/checkpoint_t_10.00ms.h5", 1e-10, &why), why);
}

TEST_CASE("monodomain: the pseudo-ECG writes one signal per lead")
{
  const BoxSpec b = slab(10);
  const std::string extra =
    "[ecg]\nsigma_b = 20.0\n"
    "leads = [ { x = -5000.0, y = 500.0, z = 500.0, name = \"a\" },\n"
    "          { x = 25000.0, y = 500.0, z = 500.0, name = \"b\" } ]\n";
  run_monodomain("ecg", b, mono_toml("TT2", 10.0, 0.02, 0.5, "isotropic",
                                     "max = [1000.0, 1000.0, 1000.0]", extra));
  std::string ecg;
  for (const auto & e : fs::directory_iterator("ecg"))
    if (e.path().filename().string().find("ecg") != std::string::npos)
      ecg = e.path().string();
  REQUIRE(!ecg.empty());
  const arma::mat table = read_table(ecg);
  REQUIRE(table.n_rows >= 5);
  CHECK(table.is_finite());
  // The time column and one column for each lead
  CHECK_EQ(table.n_cols, arma::uword(3));
  CHECK(arma::abs(table.col(1)).max() > 0.0);
  // Leads on opposite sides of a travelling wave see opposite signs
  CHECK(arma::dot(table.col(1), table.col(2)) < 0.0);
}

// ---- Purkinje ---------------------------------------------------------------

namespace
{

// A line that forks, with elements of about 100 um as the conductivity of the Purkinje model
// expects: from x = 0 to 1000 along y = z = 500, then two branches that end at the tissue
// nodes (2000, 0, 500) and (2000, 1000, 500) of a slab(20).
void write_purkinje_tree(const std::string & path)
{
  std::vector<arma::vec3> pts;
  for (int i = 0; i <= 10; i++) pts.push_back(arma::vec3({100.0 * i, 500.0, 500.0}));
  for (int i = 1; i <= 10; i++) pts.push_back(arma::vec3({1000.0 + 100.0 * i, 500.0 - 50.0 * i, 500.0}));
  for (int i = 1; i <= 10; i++) pts.push_back(arma::vec3({1000.0 + 100.0 * i, 500.0 + 50.0 * i, 500.0}));

  std::ostringstream os;
  os << "<?xml version=\"1.0\"?>\n<mesh celltype=\"line\" dim=\"1\"><nodes size=\"" << pts.size() << "\">\n";
  for (std::size_t i = 0; i < pts.size(); i++)
    os << "<node id=\"" << i << "\" x=\"" << pts[i](0) << "\" y=\"" << pts[i](1)
       << "\" z=\"" << pts[i](2) << "\"/>\n";
  os << "</nodes><elements size=\"30\">\n";
  int e = 0;
  auto seg = [&](int a, int b) { os << "<element id=\"" << e++ << "\" v0=\"" << a << "\" v1=\"" << b << "\"/>\n"; };
  for (int i = 0; i < 10; i++) seg(i, i + 1);
  seg(10, 11);
  for (int i = 11; i < 20; i++) seg(i, i + 1);
  seg(10, 21);
  for (int i = 21; i < 30; i++) seg(i, i + 1);
  os << "</elements></mesh>\n";
  write_text(path, os.str());
}

} // namespace

TEST_CASE("monodomain -ep purkinje: the wave crosses the network")
{
  fs::create_directories("pk");
  write_purkinje_tree("pk/t.xml");
  write_text("pk/c.toml",
    "[problem]\nmesh = \"t.xml\"\nmodel = \"purkinje\"\n"
    "[time]\ndt = 0.02\nt_end = 30.0\ndt_output = 1.0\n"
    "[cell_model]\nname = \"TT2\"\n"
    "[numerics]\node_solver = \"ExplicitEuler\"\n"
    "[stimuli]\npurkinje = [ { start = 0.0, duration = 2.0, value = -500.0, "
    "min = [0.0, 0.0, 0.0], max = [300.0, 1000.0, 1000.0] } ]\n");
  const AppResult r = run_app("monodomain", std::string("-f c.toml -num_threads 1") + kLU, "pk");
  CHECK_MSG(r.exit_code == 0, r.log.substr(r.log.size() > 800 ? r.log.size() - 800 : 0));
  REQUIRE(fs::exists("pk/t.xml_output.h5"));
  const arma::mat vm = h5_read("pk/t.xml_output.h5", "/point_data/vm");
  CHECK(vm.is_finite());
  // Terminal nodes of both branches activate
  const arma::vec at = activation_step(vm);
  CHECK_MSG(at.is_finite(), "network not fully activated");
}

TEST_CASE("monodomain -ep coupled_purkinje: the tissue activates through the network")
{
  fs::create_directories("pkt");
  write_purkinje_tree("pkt/t.xml");
  BoxSpec b = slab(20);
  write_hex_box("pkt/m.xml", b);
  write_text("pkt/c.toml",
    "[problem]\nmesh = \"m.xml\"\npurkinje_mesh = \"t.xml\"\nmodel = \"coupled_purkinje\"\n"
    "[time]\ndt = 0.02\nt_end = 40.0\ndt_output = 1.0\n"
    "[cell_model]\nname = \"TT2\"\n"
    "[numerics]\node_solver = \"ExplicitEuler\"\ntheta = 0.5\n"
    "[physical]\nconductivity_type = \"isotropic\"\nsigma_f = 0.0005\nsigma_s = 0.0005\nsigma_n = 0.0005\n"
    "[stimuli]\npurkinje = [ { start = 0.0, duration = 2.0, value = -500.0, "
    "min = [0.0, 0.0, 0.0], max = [300.0, 1000.0, 1000.0] } ]\n");
  const AppResult r = run_app("monodomain", std::string("-f c.toml -num_threads 1") + kLU, "pkt");
  CHECK_MSG(r.exit_code == 0, r.log.substr(r.log.size() > 800 ? r.log.size() - 800 : 0));

  REQUIRE(fs::exists("pkt/m_output.h5"));
  const arma::mat vm = h5_read("pkt/m_output.h5", "/point_data/vm");
  CHECK(vm.is_finite());
  // The stimulus is only in the network, and the end of the tissue gets activated
  CHECK_MSG(vm.col(box_node(b, b.nx, 1, 1)).max() > 0.0,
            "the tissue never activated: " << vm.max());
}

// ---- Bidomain ---------------------------------------------------------------

TEST_CASE_KNOWN_BUG("bidomain: propagation and a finite extracellular potential",
                    "the bidomain executable aborts at the first output step with an HDF5 'no output "
                    "buffer' error, thrown as HighFive::DataSetException from the writer")
{
  BoxSpec b = slab(10);
  b.Lx = 5000.0;
  b.extra_mesh_xml = "";
  fs::create_directories("bido");
  write_hex_box("bido/m.xml", b);
  std::string xml = read_text("bido/m.xml");
  xml += "<electrophysiology><stimuli number=\"1\">"
         "<stim start=\"0.0\" duration=\"2.0\" value=\"-52.0\" x0=\"0\" x1=\"500\" "
         "y0=\"0\" y1=\"1000\" z0=\"0\" z1=\"1000\"/></stimuli></electrophysiology>\n";
  write_text("bido/m.xml", xml);

  const AppResult r = run_app("bidomain",
      // The elliptic matrix has the constants in its kernel: a direct solver hits a zero pivot
      "-f m.xml -dt 0.02 -t 15 -pr 1 -pa 1 -c TT2 -m ExplicitEuler -linsolver petsc -ksp_type cg -pc_type jacobi",
      "bido");
  CHECK_MSG(r.exit_code == 0, r.log.substr(r.log.size() > 800 ? r.log.size() - 800 : 0));

  REQUIRE(fs::exists("bido/m_output.h5"));
  const arma::mat vm = h5_read("bido/m_output.h5", "/point_data/vm");
  CHECK(vm.is_finite());
  CHECK_MSG(vm.col(box_node(b, b.nx, 1, 1)).max() > 0.0, "no propagation: " << vm.max());
  if (h5_has("bido/m_output.h5", "/point_data/ve"))
    CHECK(h5_read("bido/m_output.h5", "/point_data/ve").is_finite());
  CHECK_GOLDEN("bidomain.TT2.vm_end",
               vm(vm.n_rows - 1, box_node(b, b.nx, 1, 1)), 1e-5);
}
