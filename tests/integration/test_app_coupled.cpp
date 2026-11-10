/*
 * Integration tests of the coupled problems: electromech and cardiac_cycle.
 *
 * They run on small meshes from the repository, for a short time, with the PETSc direct
 * solver. The biventricular mesh of tests/data/teste_pvloop_save_state has 702 nodes.
 *
 * Regression baselines use -ta_curve, which tests mechanics and circulation by themselves.
 *
 * What is covered and why:
 *  - electromech: the cook's beam case runs and deforms.
 *  - cardiac_cycle with each circulation model:
 *      legacy (tabulated curve): the prescribed pressure reaches the cavity,
 *      -circ 1 (Regazzoni): total blood volume conserved, the 3D volume follows the 0D one,
 *      -pw 1 (piecewise five-phase): the phases advance in order, the pressure ramp
 *      of PRELOAD is prescribed, the isovolumic contraction keeps its volume and the
 *      filling follows the prescribed inflow.
 *  - Ta unit: the mechanics receive Pa. The ceiling -t only reduces beats.
 *  - -ta_curve replaces the active tension of the cell model.
 *  - Mechano-electric feedback (-lam, -lamrate): the stretch reaches the cell.
 *  - The flag table with its provenance. The layout of the output files is not tested.
 *  - Checkpoint: for each circulation model, a run restarted from -restore ends in the
 *    same state as the uninterrupted one (the check of tests/data/teste_pvloop_save_state,
 *    made automatic).
 *  - Invalid flag combinations are rejected.
 *  - Regression baselines of the volumes and pressures reached, for each model.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <filesystem>

using namespace cardiax_test;
namespace fs = std::filesystem;

namespace
{

const char * kMesh = "tests/data/teste_pvloop_save_state/malha_descarregada_ldrb_ms.xml";
const char * kMeshName = "malha_descarregada_ldrb_ms.xml";
const char * kCommon = " -num_threads 2 -linsolver petsc -ksp_type preonly -pc_type lu"
                       " -c ToRORdLand -dt 0.02 -dt_mech 1 -logc 1000 ";

AppResult cycle(const std::string & dir, const std::string & flags)
{
  copy_input(kMesh, dir);
  return run_app("cardiac_cycle", std::string("-f ") + kMeshName + kCommon + flags, dir);
}

// Flags are read at their first occurrence, so this one leaves the choice of cell model
// and of coupling step to the caller
AppResult cycle_bare(const std::string & dir, const std::string & flags)
{
  copy_input(kMesh, dir);
  return run_app("cardiac_cycle",
                 std::string("-f ") + kMeshName +
                 " -num_threads 2 -linsolver petsc -ksp_type preonly -pc_type lu -logc 1000 " + flags, dir);
}

// One thread: the element loops add to the residual in whatever order the threads arrive, so
// two runs with several threads differ by roundoff, and only one thread can be compared exactly
AppResult cycle_serial(const std::string & dir, const std::string & flags)
{
  copy_input(kMesh, dir);
  return run_app("cardiac_cycle",
                 std::string("-f ") + kMeshName +
                 " -num_threads 1 -linsolver petsc -ksp_type preonly -pc_type lu -c ToRORdLand"
                 " -dt 0.02 -dt_mech 1 -logc 1000 " + flags, dir);
}

std::string tail(const std::string & log, std::size_t n = 1200)
{
  return log.substr(log.size() > n ? log.size() - n : 0);
}

// The mechanics writer holds the frames (fields of cardiac_cycle)
std::string frames(const std::string & dir)
{
  const std::string base = kMeshName;
  return dir + "/" + base.substr(0, base.rfind(".xml")) + "_output_nl.h5";
}

// A smooth twitch of 15 kPa over one period, zero at both ends. A jump in Ta when the table
// wraps around would invert elements in one step.
void write_ta_curve(const std::string & dir, double period = 100.0)
{
  fs::create_directories(dir);
  std::ostringstream os;
  for (int i = 0; i <= 40; i++)
  {
    const double t = period * i / 40.0;
    os << t << " " << 15000.0 * std::pow(std::sin(M_PI * t / period), 2) << "\n";
  }
  write_text(dir + "/ta.txt", os.str());
}

std::string history(const std::string & dir, const std::string & suffix)
{
  return dir + "/" + kMeshName + suffix;
}

// The checkpoints of the full run and of the restarted one must hold the same state.
// The output frame bookkeeping is not physics, so it is not compared.
void check_checkpoints_match(const std::string & full, const std::string & rest,
                             const std::vector<std::string> & times)
{
  for (const std::string & t : times)
  {
    const std::string name = "/checkpoint_t_" + t + "ms.h5";
    REQUIRE(fs::exists(full + name));
    REQUIRE(fs::exists(rest + name));
    std::string why;
    CHECK_MSG(h5_equal(full + name, rest + name, 1e-12, &why), "t = " << t << ": " << why);
  }
}

} // namespace

// ---- electromech ------------------------------------------------------------

TEST_CASE("electromech: the cook's beam case runs and deforms")
{
  const std::string name = copy_input("tests/data/mechanics/cooks_modificado_electromech.xml", "em");
  const AppResult r = run_app("electromech",
      "-f " + name + " -c ToRORdLand -m ExplicitEuler -ep mono -dt 0.02 -dt_mech 1.0 -t 30 "
      "-num_threads 2 -linsolver petsc -ksp_type preonly -pc_type lu", "em");
  CHECK_MSG(r.exit_code == 0, tail(r.log));

  const std::string stem = "em/" + name.substr(0, name.rfind(".xml"));
  REQUIRE(fs::exists(stem + "_output_nl.h5"));
  const arma::mat u = h5_read(stem + "_output_nl.h5", "/point_data/displacements");
  CHECK(u.is_finite());
  CHECK_MSG(arma::abs(u).max() > 0.0, "no displacement");
}

TEST_CASE("electromech: an unknown EP model is an error")
{
  const std::string name = copy_input("tests/data/mechanics/cooks_modificado_electromech.xml", "em_bad");
  const AppResult r = run_app("electromech", "-f " + name + " -ep nope -t 1", "em_bad");
  CHECK(r.exit_code != 0);
}

// ---- cardiac_cycle, legacy --------------------------------------------------

TEST_CASE("cardiac_cycle legacy: the curve drives the cavity and the volume responds")
{
  write_ta_curve("legacy");
  const AppResult r = cycle("legacy", "-t 30 -ta_curve ta.txt");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  CHECK(r.log.find("flags desta rodada") != std::string::npos);

  const arma::mat h = read_table(history("legacy", "_pvloop.txt"));
  REQUIRE(h.n_rows >= 25);
  CHECK(h.is_finite());
  CHECK_GOLDEN("cycle.legacy.vlv_last", h(h.n_rows - 1, 2), 1e-5);
  CHECK_GOLDEN("cycle.legacy.plv_last", h(h.n_rows - 1, 1), 1e-5);
}

// ---- cardiac_cycle, Regazzoni -----------------------------------------------

TEST_CASE("cardiac_cycle -circ 1: blood volume is conserved, and the volumes are finite")
{
  write_ta_curve("circ", 40.0);
  const AppResult r = cycle("circ", "-circ 1 -rr 40 -beats 1 -ta_curve ta.txt");
  CHECK_MSG(r.exit_code == 0, tail(r.log));

  const arma::mat h = read_table(history("circ", "_circulation.dat"));
  REQUIRE(h.n_rows == 40);
  CHECK(h.is_finite());

  const arma::vec vtot = h.col(h.n_cols - 1);
  CHECK_MSG(vtot.max() - vtot.min() < 1e-3 * vtot(0), "volume drift " << vtot.max() - vtot.min());

  // V_LV and p_LV: the mesh volume was measured, positive
  CHECK(h(0, 1) > 0.0);
  CHECK_GOLDEN("cycle.circ.vlv_last", h(h.n_rows - 1, 1), 1e-4);
  CHECK_GOLDEN("cycle.circ.plv_last", h(h.n_rows - 1, 2), 1e-4);

  CHECK(r.log.find("Total 3D solves") != std::string::npos);
  const double solves = number_after(r.log, "Total 3D solves:");
  CHECK_MSG(solves >= 40.0, solves);
}

TEST_CASE("cardiac_cycle -circ 1: -t is a ceiling that reduces the beats")
{
  const AppResult r = cycle("ceiling", "-circ 1 -rr 20 -beats 3 -t 45");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  CHECK(r.log.find("reduzido para 2") != std::string::npos);
  const arma::mat h = read_table(history("ceiling", "_circulation.dat"));
  CHECK_EQ(h.n_rows, arma::uword(40));
}

TEST_CASE("cardiac_cycle: the result does not depend on the number of threads")
{
  copy_input(kMesh, "th1");
  copy_input(kMesh, "th2");
  const std::string f = std::string("-f ") + kMeshName +
      " -linsolver petsc -ksp_type preonly -pc_type lu -c ToRORdLand -dt 0.02 -dt_mech 1 "
      "-logc 1000 -circ 1 -rr 10 -beats 1 -num_threads ";
  REQUIRE(run_app("cardiac_cycle", f + "1", "th1").exit_code == 0);
  REQUIRE(run_app("cardiac_cycle", f + "2", "th2").exit_code == 0);
  const arma::mat a = read_table(history("th1", "_circulation.dat"));
  const arma::mat b = read_table(history("th2", "_circulation.dat"));
  REQUIRE(a.n_rows == b.n_rows);
  // Pressure comes from a Newton solve with a stopping tolerance
  CHECK_MSG(arma::abs(a - b).max() < 1e-2, "max gap " << arma::abs(a - b).max());
}

// ---- cardiac_cycle, piecewise -----------------------------------------------

TEST_CASE("cardiac_cycle -pw 1: phases go in order, and the volume constraints are met")
{
  // Short period: PRELOAD ramp up to 15 ms, then IVC. The small mesh
  // only reaches ~1300 Pa, so the arterial pressures are lowered for the valves to open
  write_ta_curve("pw", 40.0);
  const AppResult r = cycle("pw", "-pw 1 -pwinit 10 -pwedp 5 -rr 40 -beats 1"
                                  " -pwpej 1150 -pwpej_rv 450 -ta_curve ta.txt");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  CHECK(r.log.find("passo(s) terminaram sem satisfazer") == std::string::npos);

  // Columns: t, then per cavity: phase, P, V, Q, P_c (phase under which the step was solved)
  const arma::mat h = read_table(history("pw", "_piecewise.dat"));
  REQUIRE(h.n_rows == 40);
  CHECK(h.is_finite());

  // The phase of each cavity never goes back within the beat
  for (int k : {1, 6})
  {
    for (arma::uword i = 1; i < h.n_rows; i++)
      CHECK_MSG(h(i, k) >= h(i - 1, k), "phase went from " << h(i - 1, k) << " to " << h(i, k));
    CHECK_NEAR(h(0, k), 0.0, 0.0);
    CHECK_MSG(h(h.n_rows - 1, k) == 4.0, "did not reach FILLING: " << h(h.n_rows - 1, k));
  }

  // PRELOAD is a prescribed pressure ramp: 0 -> 500 Pa at t_zero = 10 ms (LV defaults)
  CHECK_NEAR(h(9, 2), 500.0, 1e-6);
  CHECK_MSG(h(9, 2) > h(0, 2), "pressure did not grow in the ramp");

  // Isovolumic contraction keeps the end-diastolic volume, within -pwtol (1e-3 mL)
  for (int k : {1, 6})
  {
    double v_ed = 0.0;
    for (arma::uword i = 0; i < h.n_rows; i++)
    {
      if (h(i, k) == 0.0) v_ed = h(i, k + 2);
      if (h(i, k) == 1.0)
        CHECK_MSG(std::fabs(h(i, k + 2) - v_ed) < 1.1e-3,
                  "t = " << h(i, 0) << ": V moved " << h(i, k + 2) - v_ed << " mL in IVC");
    }
  }
  // FILLING prescribes the inflow (defaults 0.104 and 0.057 mL/ms)
  const arma::uword n = h.n_rows;
  CHECK_NEAR(h(n - 1, 3) - h(n - 2, 3), 0.104, 1.1e-3);
  CHECK_NEAR(h(n - 1, 8) - h(n - 2, 8), 0.057, 1.1e-3);
  CHECK_GOLDEN("cycle.pw.vlv_last", h(n - 1, 3), 1e-5);
}

// ---- Active tension ---------------------------------------------------------

TEST_CASE("cardiac_cycle: -ta_curve prescribes the active tension in Pa")
{
  fs::create_directories("tac");
  write_text("tac/ta.txt", "0 0\n10 5000\n20 10000\n");
  const AppResult r = cycle("tac", "-t 20 -ta_curve ta.txt -prc 1");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  CHECK(r.log.find("PRESCRIBED") != std::string::npos);

  const std::string h5 = frames("tac");
  REQUIRE(fs::exists(h5));
  const arma::mat ta = h5_read(h5, "/point_data/active_stress");
  // The table goes 500 Pa per ms at every node, one frame per ms (the writer allocates more
  // frames than the run fills, so only the first ones are read)
  REQUIRE(ta.n_rows >= 21);
  for (int k = 0; k <= 20; k++)
  {
    const arma::vec row = arma::vectorise(ta.row(k));
    CHECK_NEAR(row.max(), 500.0 * k, 1e-6);
    CHECK_NEAR(row.min(), 500.0 * k, 1e-6);
  }
}

TEST_CASE("cardiac_cycle: -tascale scales the tension from the cell model")
{
  auto peak = [](const std::string & dir, const std::string & flags) {
    fs::create_directories(dir);
    write_text(dir + "/ta.txt", "0 1000\n20 1000\n");
    const AppResult r = cycle(dir, "-t 10 -ta_curve ta.txt -prc 1 " + flags);
    REQUIRE(r.exit_code == 0);
    const std::string h5 = frames(dir);
    return h5_read(h5, "/point_data/active_stress").max();
  };
  CHECK_REL(peak("sc2", "-tascale 2"), 2.0 * peak("sc1", ""), 1e-9);
}

// ---- Mechano-electric feedback ----------------------------------------------

TEST_CASE("cardiac_cycle: with feedback on, the stretch reaches the cell model")
{
  const AppResult off = cycle("mef_off", "-circ 1 -rr 20 -beats 1 -prc 1");
  const AppResult on  = cycle("mef_on", "-circ 1 -rr 20 -beats 1 -prc 1 -lam 1 -lamrate 1");
  CHECK_MSG(off.exit_code == 0, tail(off.log));
  CHECK_MSG(on.exit_code == 0, tail(on.log));
  CHECK(on.log.find("Mechano-electric feedback ON") != std::string::npos);
  CHECK(off.log.find("Mechano-electric feedback off") != std::string::npos);

  const std::string h5 = frames("mef_on");
  REQUIRE(fs::exists(h5));
  const arma::mat lam = h5_read(h5, "/point_data/lambda_f");
  CHECK(lam.is_finite());
  CHECK_MSG(arma::abs(lam.row(lam.n_rows - 2) - 1.0).max() > 1e-6, "lambda stayed at 1");
  CHECK_MSG(lam.max() <= 1.3 && lam.max() > 0.0, "out of the clip: " << lam.max());

  const arma::mat a = read_table(history("mef_off", "_circulation.dat"));
  const arma::mat b = read_table(history("mef_on", "_circulation.dat"));
  REQUIRE(a.n_rows == b.n_rows);
  CHECK(b.is_finite());
}

TEST_CASE("cardiac_cycle: feedback on a model that ignores stretch is turned off, with a message")
{
  const AppResult r = cycle_bare("mef_kerck", "-t 5 -c Kerkoff2003 -lam 1");
  CHECK(r.log.find("coupling disabled") != std::string::npos);
}

// ---- Cell models ------------------------------------------------------------

TEST_CASE("cardiac_cycle: ToRORd-Land gives a real membrane potential and twitch")
{
  // This mesh has no activation times, so its cells would never be stimulated: all nodes
  // activate at t = 0 here. A Ta of 1 kPa is reached well before the end of the run.
  copy_input(kMesh, "torord");
  std::ostringstream lat;
  lat << "<eikonal>\n";
  for (int i = 0; i < 702; i++) lat << "<node id=\"" << i << "\" lat=\"0.0\"/>\n";
  lat << "</eikonal>\n</mesh>";
  std::string xml = read_text(std::string("torord/") + kMeshName);
  xml.replace(xml.rfind("</mesh>"), 7, lat.str());
  write_text(std::string("torord/") + kMeshName, xml);
  const AppResult r = run_app("cardiac_cycle",
      std::string("-f ") + kMeshName + kCommon + "-circ 1 -rr 100 -beats 1 -prc 1", "torord");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  const arma::mat vm = h5_read(frames("torord"), "/point_data/vm");
  CHECK(vm.is_finite());
  CHECK_MSG(vm.max() > 0.0 && vm.max() < 90.0, "Vm peak " << vm.max());
  const arma::mat ta = h5_read(frames("torord"), "/point_data/active_stress");
  CHECK(ta.is_finite());
  CHECK_MSG(ta.max() > 1000.0, "Ta peak [Pa] " << ta.max());
}

TEST_CASE("cardiac_cycle: Kerckhoffs model, activated by the LAT")
{
  copy_input(kMesh, "kerck");
  const AppResult r = run_app("cardiac_cycle",
      std::string("-f ") + kMeshName +
      " -num_threads 2 -linsolver petsc -ksp_type preonly -pc_type lu -c Kerkoff2003 -t 10 -logc 1000",
      "kerck");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
}

// ---- Checkpoint -------------------------------------------------------------

TEST_CASE("cardiac_cycle: restart continues the legacy run exactly")
{
  REQUIRE(cycle_serial("chk_legacy_full", "-t 20 -save_state 5").exit_code == 0);
  fs::create_directories("chk_legacy_rest");
  fs::copy_file("chk_legacy_full/checkpoint_t_10.00ms.h5", "chk_legacy_rest/start.h5",
                fs::copy_options::overwrite_existing);
  const AppResult r = cycle_serial("chk_legacy_rest", "-t 20 -save_state 5 -restore start.h5");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  check_checkpoints_match("chk_legacy_full", "chk_legacy_rest", {"15.00", "20.00"});
}

TEST_CASE("cardiac_cycle: restart continues the -circ 1 run exactly")
{
  const std::string flags = "-circ 1 -rr 20 -beats 1 -save_state 5";
  REQUIRE(cycle_serial("chk_circ_full", flags).exit_code == 0);
  fs::create_directories("chk_circ_rest");
  fs::copy_file("chk_circ_full/checkpoint_t_10.00ms.h5", "chk_circ_rest/start.h5",
                fs::copy_options::overwrite_existing);
  const AppResult r = cycle_serial("chk_circ_rest", flags + " -restore start.h5");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  check_checkpoints_match("chk_circ_full", "chk_circ_rest", {"15.00", "20.00"});
}

TEST_CASE("cardiac_cycle: restart continues the -circ 1 run with feedback exactly")
{
  const std::string flags = "-circ 1 -rr 20 -beats 1 -lam 1 -lamrate 1 -save_state 5";
  REQUIRE(cycle_serial("chk_mef_full", flags).exit_code == 0);
  fs::create_directories("chk_mef_rest");
  fs::copy_file("chk_mef_full/checkpoint_t_10.00ms.h5", "chk_mef_rest/start.h5",
                fs::copy_options::overwrite_existing);
  const AppResult r = cycle_serial("chk_mef_rest", flags + " -restore start.h5");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  check_checkpoints_match("chk_mef_full", "chk_mef_rest", {"15.00", "20.00"});
}

TEST_CASE("cardiac_cycle: restart continues the -pw 1 run exactly")
{
  // Ta from a table: with the cell model's Ta and a period this short, the mechanics diverges
  const std::string flags = "-pw 1 -pwinit 10 -pwedp 5 -rr 40 -beats 1 -pwpej 1150 -pwpej_rv 450"
                            " -save_state 5 -ta_curve ta.txt";
  write_ta_curve("chk_pw_full", 40.0);
  write_ta_curve("chk_pw_rest", 40.0);
  REQUIRE(cycle_serial("chk_pw_full", flags).exit_code == 0);
  // Cut in the isovolumic contraction: the restart carries the Newton state and the
  // Windkessel through the ejection
  fs::copy_file("chk_pw_full/checkpoint_t_20.00ms.h5", "chk_pw_rest/start.h5",
                fs::copy_options::overwrite_existing);
  const AppResult r = cycle_serial("chk_pw_rest", flags + " -restore start.h5");
  CHECK_MSG(r.exit_code == 0, tail(r.log));
  check_checkpoints_match("chk_pw_full", "chk_pw_rest", {"25.00", "30.00", "35.00", "40.00"});
}

// ---- Invalid input ----------------------------------------------------------

TEST_CASE("cardiac_cycle: invalid flag combinations are rejected")
{
  CHECK(cycle("bad1", "-circ 1 -pw 1").exit_code != 0);
  CHECK(cycle_bare("bad2", "-c ToRORdLand -circ 1 -dt_mech 0").exit_code != 0);

  copy_input(kMesh, "bad3");
  const AppResult r = run_app("cardiac_cycle",
      std::string("-f ") + kMeshName + " -c TT2 -t 1 -linsolver petsc", "bad3");
  CHECK(r.exit_code != 0);
  CHECK(r.log.find("sem mapeamento") != std::string::npos);
}

TEST_CASE("cardiac_cycle: -h prints the help and leaves the output directory alone")
{
  fs::create_directories("help/output");
  write_text("help/output/keep.txt", "x");
  const AppResult r = run_app("cardiac_cycle", "-h", "help");
  CHECK_EQ(r.exit_code, 0);
  CHECK(r.log.find("Usage: cardiac_cycle") != std::string::npos);
  CHECK(fs::exists("help/output/keep.txt"));
}
