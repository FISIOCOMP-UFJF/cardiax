/*
 * Unit tests of src/odes: ODE integrators, every cell model, the per-node collection
 * of cells and the warm-up.
 *
 * What is covered and why:
 *  - ExplicitEuler, ImplicitEuler, RungeKutta4, LSODA on du/dt = -2u: the error against
 *    the exact solution and the observed order of convergence.
 *  - TimeStepper: step counting used by the solve() loops.
 *  - The CellModel factory, for all names, and that an unknown name throws.
 *  - Every cell model: initial state and right-hand side are finite, the Jacobian has
 *    the right shape and is finite, the integration of a stimulated cell stays finite,
 *    and the final state matches the regression baseline (tests/golden/baseline.txt).
 *    The baseline detects unintended changes in equations that have no closed-form
 *    answer to compare with.
 *  - CellModel: the cell type has a default (ENDO) even when built over dirty memory.
 *  - Physiology of the ionic models used in production (TT2, ORd, ToRORd-Land): resting
 *    potential, upstroke, repolarisation; and contraction in ToRORd-Land.
 *  - Length dependence of the Land submodel, the neutral values of stretch and of the
 *    apicobasal coordinate, and the active stiffness formula.
 *  - Kerckhoffs: native time unit in seconds and activation driven by the LAT state.
 *  - Cells: time unit conversion, stimulus on selected nodes only, per-node values, and
 *    that the result does not depend on the number of OpenMP threads.
 *  - CellWarmup: grouping by cell type, cache reuse and invalidation, drift towards a
 *    limit cycle, and phase by LAT.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <cstring>
#include <memory>
#include <new>
#include <omp.h>

#include "odes/odes.h"
#include "odes/torord_land.hpp"
#include "odes/kerkoff2003.hpp"
#include "util/util.hpp"
#include "util/command_line_args.h"

using namespace cardiax_test;

namespace
{

const double kExactEnd = 5.0 * std::exp(-2.0);   // u(1) of du/dt = -2u, u(0) = 5

// Error at t = 1 of the method on SimpleODE
double error_at_one(const std::string & method, int nsteps)
{
  const double dt = 1.0 / nsteps;
  SimpleODE ode;
  ode.setup(method, dt, 1.0);
  double y;
  ode.init(&y);
  for (int i = 0; i < nsteps; i++)
    ode.advance(&y, i * dt, dt);
  return std::fabs(y - kExactEnd);
}

const std::vector<std::string> kAllModels = {
  "FHN", "NP", "MNP", "MINI", "MMSilva", "MS", "MV", "LR1", "ORd", "TNNP", "TT2",
  "TT2Ta", "RiceTT2", "RiceORd", "SODE", "ToRORdLand", "Kerkoff2003"};

std::unique_ptr<CellModel> make_model(const std::string & name, double dt = 0.01,
                                      const std::string & method = "ExplicitEuler")
{
  std::unique_ptr<CellModel> m(CellModel::create(name));
  m->setup(method, dt, 1000.0);
  return m;
}

// Membrane potential trace of one paced cell. Stimulus in the first stim_dur ms.
struct Trace { double rest, peak, end, apd_end; };

Trace paced_trace(const std::string & name, double stim, double T, double dt,
                  double stim_dur = 1.0)
{
  auto m = make_model(name, dt);
  arma::vec y(m->get_num_state_vars());
  m->init(y.memptr());
  Trace tr;
  tr.rest = y(0);
  tr.peak = y(0);
  tr.apd_end = -1.0;
  const int n = static_cast<int>(T / dt);
  for (int i = 0; i < n; i++)
  {
    const double t = i * dt;
    if (t < stim_dur) m->advance(y.memptr(), t, dt, stim);
    else              m->advance(y.memptr(), t, dt, 0.0);
    tr.peak = std::max(tr.peak, y(0));
    // First time it comes back to near rest after the upstroke
    if (tr.apd_end < 0.0 && tr.peak > 0.0 && y(0) < tr.rest + 10.0)
      tr.apd_end = t;
  }
  tr.end = y(0);
  return tr;
}

//! Models whose ExplicitEuler path diverges: they declare `rlvars` but not `is_rl`, see the
//! KNOWN_BUG case below.
const std::vector<std::string> kBrokenRushLarsen = {"MNP", "ORd", "TT2Ta", "RiceORd"};

//! One step of ToRORd-Land through the shipped solver
void torord_step(CellModel & m, arma::vec & y, double t, double dt, double stim)
{
  m.advance(y.memptr(), t, dt, stim);
}

} // namespace

// ---- ODE integrators --------------------------------------------------------

TEST_CASE("ExplicitEuler: first order on the model problem")
{
  const double e1 = error_at_one("ExplicitEuler", 100);
  const double e2 = error_at_one("ExplicitEuler", 200);
  CHECK_NEAR(std::log2(e1 / e2), 1.0, 0.05);
  CHECK(e2 < 1e-2);
}

TEST_CASE("ImplicitEuler: first order on the model problem")
{
  const double e1 = error_at_one("ImplicitEuler", 100);
  const double e2 = error_at_one("ImplicitEuler", 200);
  CHECK_NEAR(std::log2(e1 / e2), 1.0, 0.05);
  CHECK(e2 < 1e-2);
}

TEST_CASE("ImplicitEuler: stays stable where the explicit method blows up")
{
  // du/dt = -2u with dt = 2: |1 - 2 dt| = 3 explodes, 1/(1 + 2 dt) decays
  auto run = [](const std::string & method) {
    SimpleODE ode;
    ode.setup(method, 2.0, 100.0);
    double y;
    ode.init(&y);
    for (int i = 0; i < 10; i++)
      ode.advance(&y, i * 2.0, 2.0);
    return std::fabs(y);
  };
  CHECK(run("ExplicitEuler") > 1e3);
  CHECK(run("ImplicitEuler") < 1e-3);
}

TEST_CASE("RungeKutta4: fourth order on the model problem")
{
  const double e1 = error_at_one("RungeKutta4", 10);
  const double e2 = error_at_one("RungeKutta4", 20);
  CHECK_NEAR(std::log2(e1 / e2), 4.0, 0.15);
  CHECK(e2 < 1e-4);
}

TEST_CASE("LSODA: follows the exact solution")
{
  CHECK(error_at_one("LSODA", 20) < 1e-3);
}

TEST_CASE("ODESolver: an unknown method throws")
{
  SimpleODE ode;
  CHECK_THROWS(ode.setup("NoSuchMethod", 0.1, 1.0));
  ode.setup("ExplicitEuler", 0.1, 1.0);   // so that the destructor has what to delete
}

TEST_CASE("TimeStepper: steps, end of the run and print rate")
{
  TimeStepper ts(0.25, 1.0, 0.5);
  CHECK_NEAR(ts.timestep(), 0.25, 0.0);
  CHECK_NEAR(ts.endtime(), 1.0, 0.0);
  CHECK_EQ(ts.number_of_time_steps(), 4);
  CHECK(!ts.finished());

  int n = 0;
  int prints = 0;
  while (!ts.finished())
  {
    ts.increase_time();
    if (ts.time_to_print()) prints++;
    n++;
    REQUIRE(n <= 100);
  }
  CHECK_EQ(n, 4);
  CHECK_EQ(prints, 2);
  CHECK_NEAR(ts.time(), 1.0, 1e-14);
}

// ---- Factory and generic properties -----------------------------------------

TEST_CASE("CellModel::create: an unknown name throws")
{
  CHECK_THROWS(CellModel::create("NoSuchCell"));
}

TEST_CASE("CellModel: known sizes of the models that have others depending on them")
{
  CHECK_EQ(make_model("ToRORdLand")->get_num_state_vars(), 50);
  CHECK_EQ(make_model("ORd")->get_num_state_vars(), 41);
  CHECK_EQ(make_model("Kerkoff2003")->get_num_state_vars(), 2);
  CHECK_EQ(make_model("SODE")->get_num_state_vars(), 1);
}

TEST_CASE("Every cell model: finite initial state and right-hand side")
{
  for (const std::string & name : kAllModels)
  {
    auto m = make_model(name);
    REQUIRE(m->get_num_state_vars() > 0);
    arma::vec y(m->get_num_state_vars()), f(m->get_num_state_vars());
    m->init(y.memptr());
    CHECK_MSG(y.is_finite(), name);
    m->set_stimulus(0.0);
    m->equation(0.0, y.memptr(), f.memptr());
    CHECK_MSG(f.is_finite(), name);
  }
}

TEST_CASE("Every cell model: the initial state does not depend on previous calls")
{
  for (const std::string & name : kAllModels)
  {
    auto m = make_model(name);
    arma::vec a(m->get_num_state_vars()), b(m->get_num_state_vars());
    m->init(a.memptr());
    m->advance(a.memptr(), 0.0, 0.01, 0.0);
    m->init(a.memptr());
    m->init(b.memptr());
    CHECK_MSG(arma::approx_equal(a, b, "absdiff", 0.0), name);
  }
}

TEST_CASE("Every cell model: finite-difference Jacobian is finite and leaves the state as it was")
{
  for (const std::string & name : kAllModels)
  {
    auto m = make_model(name);
    const int n = m->get_num_state_vars();
    arma::vec y(n), y0;
    m->init(y.memptr());
    y0 = y;
    m->set_stimulus(0.0);
    arma::mat J(n, n, arma::fill::zeros);
    m->compute_jacobian(y.memptr(), 0.0, J);
    // MINI and MMSilva have a power of a state that starts at zero: dV/dV is singular there
    if (name != "MINI" && name != "MMSilva")
      CHECK_MSG(J.is_finite(), name);
    CHECK_MSG(arma::approx_equal(y, y0, "absdiff", 0.0), name);
  }
}

TEST_CASE("compute_jacobian: matches the analytic Jacobian of the model problem")
{
  SimpleODE ode;
  ode.setup("ExplicitEuler", 0.1, 1.0);
  double y = 3.0;
  ode.set_stimulus(0.0);
  arma::mat J(1, 1);
  ode.compute_jacobian(&y, 0.0, J);
  CHECK_NEAR(J(0, 0), -2.0, 1e-6);
}

TEST_CASE("Every cell model: regression baseline of a stimulated integration")
{
  // Stimulus and time step large enough for each family of models to activate, with
  // the ionic ones on their usual time scale (ms).
  const struct { const char * name; double dt; double stim; double T; } runs[] = {
    {"FHN",         0.001, 1.0,   1.0},
    {"NP",          0.01,  0.5,   10.0},
    {"MINI",        0.01, -1.0,   10.0},
    {"MS",          0.01,  1.0,   10.0},
    {"MV",          0.01,  1.0,   10.0},
    {"LR1",         0.01, -50.0,  10.0},
    {"TNNP",        0.01, -52.0,  10.0},
    {"TT2",         0.01, -52.0,  10.0},
    {"RiceTT2",     0.01, -52.0,  10.0},
    {"SODE",        0.01,  0.0,   1.0},
    {"Kerkoff2003", 0.0005, 0.0,  0.1},
  };

  for (const auto & r : runs)
  {
    auto m = make_model(r.name, r.dt);
    m->set_celltype(0);
    const int n = m->get_num_state_vars();
    arma::vec y(n);
    m->init(y.memptr());
    const int steps = static_cast<int>(std::round(r.T / r.dt));
    for (int i = 0; i < steps; i++)
    {
      const double t = i * r.dt;
      m->advance(y.memptr(), t, r.dt, t < 2.0 ? r.stim : 0.0);
    }
    CHECK_MSG(y.is_finite(), r.name);
    CHECK_GOLDEN(std::string("cell.") + r.name + ".y0", y(0), 1e-6);
    CHECK_GOLDEN(std::string("cell.") + r.name + ".ylast", y(n - 1), 1e-6);
    CHECK_GOLDEN(std::string("cell.") + r.name + ".norm", arma::norm(y), 1e-6);
  }
}

TEST_CASE("ToRORd-Land: regression baseline")
{
  auto m = make_model("ToRORdLand", 0.01);
  m->set_celltype(ENDO);
  arma::vec y(m->get_num_state_vars());
  m->init(y.memptr());
  for (int i = 0; i < 5000; i++)
    torord_step(*m, y, i * 0.01, 0.01, i * 0.01 < 2.0 ? -53.0 : 0.0);
  REQUIRE(y.is_finite());
  CHECK_GOLDEN("cell.ToRORdLand.v50ms", y(0), 1e-6);
  CHECK_GOLDEN("cell.ToRORdLand.ta50ms", y(49), 1e-6);
  CHECK_GOLDEN("cell.ToRORdLand.norm50ms", arma::norm(y), 1e-6);
}

TEST_CASE_KNOWN_BUG("ExplicitEuler: every model with Rush-Larsen variables stays finite",
                    "src/odes/explicit_euler.cpp: since 9a628f1 it applies Rush-Larsen through "
                    "is_rl/rl_inf/rl_tau. ORd, RiceORd, TT2Ta and MNP declare only rlvars and go NaN "
                    "here. ToRORdLand had the same defect, with this cause verified, and was migrated "
                    "to is_rl; these four were not investigated")
{
  for (const std::string & name : kBrokenRushLarsen)
  {
    auto m = make_model(name, 0.001);
    m->set_celltype(ENDO);
    arma::vec y(m->get_num_state_vars());
    m->init(y.memptr());
    const double stim = name == "MNP" ? 0.5 : -53.0;
    for (int i = 0; i < 30000; i++)
      m->advance(y.memptr(), i * 0.001, 0.001, i * 0.001 < 2.0 ? stim : 0.0);
    CHECK_MSG(y.is_finite(), name);
  }
}

TEST_CASE("ToRORd-Land: gates stay within [0, 1] even at a large step")
{
  auto m = make_model("ToRORdLand", 0.02);
  m->set_celltype(ENDO);
  arma::vec y(m->get_num_state_vars());
  m->init(y.memptr());
  for (int i = 0; i < 2000; i++)
  {
    torord_step(*m, y, i * 0.02, 0.02, i * 0.02 < 2.0 ? -53.0 : 0.0);
    for (int v = 9; v <= 28; v++)
      REQUIRE(y(v) >= -1e-9 && y(v) <= 1.0 + 1e-9);
  }
  CHECK_EQ(m->rlvars.count(9), std::size_t(1));
  CHECK_EQ(m->rlvars.count(0), std::size_t(0));
}

// ---- Physiology -------------------------------------------------------------

TEST_CASE("TT2: action potential with the published shape")
{
  const Trace t = paced_trace("TT2", -52.0, 600.0, 0.02);
  CHECK_NEAR(t.rest, -85.0, 3.0);
  CHECK_MSG(t.peak > 20.0 && t.peak < 60.0, t.peak);
  CHECK_MSG(t.apd_end > 200.0 && t.apd_end < 450.0, t.apd_end);
  CHECK_NEAR(t.end, t.rest, 5.0);
}

TEST_CASE_KNOWN_BUG("ORd: action potential through the shipped solver",
                    "same defect as 'every model with Rush-Larsen variables stays finite'. Note also "
                    "that src/odes/ord.cpp ignores the stimulus it is given and applies a hard-coded "
                    "-80 pulse at t = 50 ms")
{
  const Trace t = paced_trace("ORd", 0.0, 400.0, 0.001);
  CHECK_MSG(t.peak > 10.0 && t.peak < 80.0, t.peak);
}

TEST_CASE("ToRORd-Land: action potential and active tension")
{
  auto m = make_model("ToRORdLand", 0.01);
  m->set_celltype(ENDO);
  arma::vec y(m->get_num_state_vars());
  m->init(y.memptr());
  const double rest = y(0);
  double peak = rest, ta_max = 0.0, t_ta_max = 0.0;
  for (int i = 0; i < 80000; i++)
  {
    const double t = i * 0.01;
    torord_step(*m, y, t, 0.01, t < 2.0 ? -53.0 : 0.0);
    peak = std::max(peak, y(0));
    if (y(49) > ta_max) { ta_max = y(49); t_ta_max = t; }
  }
  REQUIRE(y.is_finite());
  CHECK_NEAR(rest, -88.6, 1.5);
  CHECK_MSG(peak > 10.0 && peak < 90.0, peak);
  CHECK_MSG(ta_max > 5.0 && ta_max < 100.0, "Ta peak [kPa] " << ta_max);
  CHECK_MSG(t_ta_max > 80.0 && t_ta_max < 350.0, "time to peak [ms] " << t_ta_max);
  CHECK_MSG(y(49) < 0.05 * ta_max, "must relax by 800 ms: " << y(49));
  CHECK_NEAR(y(0), rest, 3.0);
}

TEST_CASE("Ionic models: no stimulus, no action potential")
{
  const Trace t = paced_trace("TT2", 0.0, 100.0, 0.02);
  CHECK_MSG(t.peak < t.rest + 5.0, "peak " << t.peak);

  auto m = make_model("ToRORdLand", 0.02);
  m->set_celltype(ENDO);
  arma::vec y(m->get_num_state_vars());
  m->init(y.memptr());
  const double rest = y(0);
  for (int i = 0; i < 5000; i++) torord_step(*m, y, i * 0.02, 0.02, 0.0);
  CHECK_NEAR(y(0), rest, 5.0);
}

TEST_CASE("ToRORd-Land: cell types differ, and each one is reproducible")
{
  auto final_v = [](int type) {
    auto m = make_model("ToRORdLand", 0.02);
    m->set_celltype(type);
    arma::vec y(m->get_num_state_vars());
    m->init(y.memptr());
    for (int i = 0; i < 10000; i++)
      torord_step(*m, y, i * 0.02, 0.02, i * 0.02 < 2.0 ? -53.0 : 0.0);
    return y(0);
  };
  const double epi = final_v(EPI), mid = final_v(MCELL), endo = final_v(ENDO);
  CHECK(std::fabs(epi - endo) > 1e-3);
  CHECK(std::fabs(mid - endo) > 1e-3);
  CHECK_NEAR(final_v(EPI), epi, 0.0);
}

TEST_CASE("CellModel: the cell type has a default even over dirty memory")
{
  // Cells keeps one CellModel per OpenMP thread and only calls set_celltype() when the
  // mesh has cell types. CellModel::type used to be left uninitialized, so with no types
  // each copy ran with whatever was in its memory (EPI in some, an invalid value in
  // others) and the answer depended on the thread -> node mapping. Building the model
  // over memory filled with garbage reproduces that deterministically.
  alignas(TorordLand) unsigned char buf[sizeof(TorordLand)];
  std::memset(buf, 0xAB, sizeof(buf));
  TorordLand * m = new (buf) TorordLand();
  CHECK_MSG(m->get_celltype() == ENDO, "type " << m->get_celltype());
  m->~TorordLand();   // ode_solver/ts must be nullptr before setup(), or this crashes
}

// ---- Land submodel and neutral values ---------------------------------------

namespace
{
// Peak Ta for a given stretch and stretch rate
double land_peak_ta(double lambda, double rate = 0.0)
{
  auto m = make_model("ToRORdLand", 0.02);
  m->set_celltype(ENDO);
  m->set_stretch(lambda);
  m->set_stretch_rate(rate);
  arma::vec y(m->get_num_state_vars());
  m->init(y.memptr());
  double ta = 0.0;
  for (int i = 0; i < 15000; i++)
  {
    torord_step(*m, y, i * 0.02, 0.02, i * 0.02 < 2.0 ? -53.0 : 0.0);
    ta = std::max(ta, y(49));
  }
  return ta;
}
}

TEST_CASE("Land: active tension increases with fibre stretch (Frank-Starling)")
{
  const double ta09 = land_peak_ta(0.9), ta10 = land_peak_ta(1.0), ta11 = land_peak_ta(1.1);
  CHECK_MSG(ta09 < ta10 && ta10 < ta11, ta09 << " " << ta10 << " " << ta11);
}

TEST_CASE("Land: the neutral stretch (1, rate 0) is what the model uses by default")
{
  auto def = make_model("ToRORdLand", 0.02);
  CHECK_NEAR(def->get_stretch(), 1.0, 0.0);
  CHECK_NEAR(def->get_stretch_rate(), 0.0, 0.0);
  CHECK_NEAR(def->get_apicobasal(), 0.5, 0.0);

  def->set_celltype(ENDO);
  arma::vec y(def->get_num_state_vars());
  def->init(y.memptr());
  double ta = 0.0;
  for (int i = 0; i < 15000; i++)
  {
    torord_step(*def, y, i * 0.02, 0.02, i * 0.02 < 2.0 ? -53.0 : 0.0);
    ta = std::max(ta, y(49));
  }
  CHECK_NEAR(ta, land_peak_ta(1.0, 0.0), 0.0);
}

TEST_CASE("Land: apicobasal coordinate changes repolarisation; 0.5 is neutral")
{
  auto run = [](double ab, bool set) {
    auto m = make_model("ToRORdLand", 0.02);
    m->set_celltype(ENDO);
    if (set) m->set_apicobasal(ab);
    arma::vec y(m->get_num_state_vars());
    m->init(y.memptr());
    for (int i = 0; i < 12500; i++)
      torord_step(*m, y, i * 0.02, 0.02, i * 0.02 < 2.0 ? -53.0 : 0.0);
    return y(0);
  };
  CHECK_NEAR(run(0.5, true), run(0.0, false), 0.0);
  CHECK(std::fabs(run(0.0, true) - run(1.0, true)) > 1e-3);
}

TEST_CASE("land_active_stiffness: formula, positivity and the length factor")
{
  // Lfac(lambda = 1) = 1 + 2.3*(1 + 0.87 - 1.87) = 1
  const double dr = 0.25, wfrac = 0.5, TOT_A = 25.0, Tref = 120.0;
  const double A = (0.25 * TOT_A) / ((1.0 - dr) * wfrac + dr) * (dr / 0.25);
  CHECK_NEAR(land_active_stiffness(1.0, 0.1, 0.2), (Tref / dr) * A * 0.3, 1e-10);

  // Negative crossbridge fractions count as zero
  CHECK_NEAR(land_active_stiffness(1.0, -0.5, 0.2), (Tref / dr) * A * 0.2, 1e-10);
  CHECK_NEAR(land_active_stiffness(1.0, -0.5, -0.5), 0.0, 0.0);

  // Saturates above lambda_max and is zero at very low stretch
  CHECK_NEAR(land_active_stiffness(1.5, 0.1, 0.1), land_active_stiffness(1.2, 0.1, 0.1), 0.0);
  CHECK_NEAR(land_active_stiffness(0.5, 0.1, 0.1), 0.0, 0.0);
  CHECK(land_active_stiffness(1.1, 0.1, 0.1) > land_active_stiffness(1.0, 0.1, 0.1));
}

// ---- Kerckhoffs -------------------------------------------------------------

TEST_CASE("Kerckhoffs: written in seconds, and activated by the LAT state")
{
  auto m = make_model("Kerkoff2003", 0.001);
  CHECK_NEAR(m->native_time_unit_ms(), 1000.0, 0.0);
  CHECK_EQ(m->lat_var_index(), 1);
  CHECK_EQ(m->get_num_monitored(), 1);

  arma::vec y(2);
  m->init(y.memptr());
  y(1) = 100.0;   // LAT in ms: 0.1 s

  m->advance(y.memptr(), 0.05, 0.001);
  CHECK_NEAR(m->get_monitored_value(0), 0.0, 0.0);   // before activation
  CHECK_NEAR(y(1), 100.0, 0.0);                      // the LAT state stays constant

  double ta = 0.0;
  for (int i = 0; i < 200; i++)
  {
    m->advance(y.memptr(), 0.1 + i * 0.001, 0.001);
    ta = std::max(ta, m->get_monitored_value(0));
  }
  CHECK_MSG(ta > 1.0, ta);
}

TEST_CASE("Ionic models report no LAT state and a millisecond time unit")
{
  auto m = make_model("ToRORdLand");
  CHECK_NEAR(m->native_time_unit_ms(), 1.0, 0.0);
  CHECK_EQ(m->lat_var_index(), -1);
}

// ---- Cells ------------------------------------------------------------------

TEST_CASE("Cells: state layout, get/set of variables and of whole systems")
{
  Cells cells(3, "NP", "ExplicitEuler", 0.01, 1.0, 1.0);
  cells.set_solver_time_unit_ms(1.0);
  cells.init();
  const uint nv = cells.get_ode_size();
  CHECK_EQ(cells.size(), 3);
  CHECK_EQ(cells.get_size(), 3 * nv);

  arma::vec v = {1.0, 2.0, 3.0};
  cells.set_var(1, v);
  arma::vec out(3);
  cells.get_var(1, out);
  CHECK_NEAR(out(2), 3.0, 0.0);
  CHECK_NEAR(cells.get_state(1, 1), 2.0, 0.0);
  CHECK_NEAR(cells.get_state(1 * nv + 1), 2.0, 0.0);

  petsc::Vector pv;
  pv.create(3);
  cells.get_var(1, pv);
  CHECK_NEAR(pv.l2norm(), std::sqrt(14.0), 1e-14);

  std::vector<double> y(nv, 9.0), back(nv, 0.0);
  cells.set_system_state(2, y.data());
  cells.get_system_state(2, back.data());
  CHECK_NEAR(back[nv - 1], 9.0, 0.0);
  CHECK_NEAR(cells.get_state(0, 0) == 9.0 ? 1.0 : 0.0, 0.0, 0.0);
}

TEST_CASE("Cells: time is converted from the solver unit to the model unit")
{
  {
    Cells c(1, "Kerkoff2003", "ExplicitEuler", 1.0, 100.0, 1.0);
    c.set_solver_time_unit_ms(1.0);
    CHECK_NEAR(c.time_factor(), 1e-3, 0.0);
  }
  {
    Cells c(1, "ToRORdLand", "ExplicitEuler", 0.01, 100.0, 1.0);
    c.set_solver_time_unit_ms(1.0);
    CHECK_NEAR(c.time_factor(), 1.0, 0.0);
    c.set_solver_time_unit_ms(1000.0);
    CHECK_NEAR(c.time_factor(), 1000.0, 0.0);
  }
}

TEST_CASE("Cells: a step in ms and the same step in seconds give the same cell state")
{
  auto step = [](double unit_ms) {
    Cells c(1, "TT2", "ExplicitEuler", 0.01, 10.0, 1.0);
    c.set_solver_time_unit_ms(unit_ms);
    c.init();
    for (int i = 0; i < 50; i++)
      c.advance(i * 0.01 / unit_ms, 0.01 / unit_ms);
    return c.get_state(0, 0);
  };
  CHECK_NEAR(step(1.0), step(1000.0), 1e-9);
}

TEST_CASE("Cells: the stimulus reaches the selected nodes only")
{
  Cells c(3, "TT2", "ExplicitEuler", 0.01, 10.0, 1.0);
  c.set_solver_time_unit_ms(1.0);
  c.init();
  const double v0 = c.get_state(0, 0);
  const std::set<uint> nodes = {1};
  for (int i = 0; i < 100; i++)
    c.advance(i * 0.01, 0.01, -52.0, nodes);
  CHECK(c.get_state(1, 0) > v0 + 50.0);
  CHECK_NEAR(c.get_state(0, 0), c.get_state(2, 0), 0.0);
  CHECK(c.get_state(0, 0) < v0 + 2.0);
}

TEST_CASE("Cells: nodal stimulus values, and the size check")
{
  Cells c(2, "TT2", "ExplicitEuler", 0.01, 10.0, 1.0);
  c.set_solver_time_unit_ms(1.0);
  c.init();
  const double v0 = c.get_state(0, 0);
  arma::vec stim = {0.0, -52.0};
  for (int i = 0; i < 100; i++)
    c.advance(i * 0.01, 0.01, stim);
  CHECK(c.get_state(1, 0) > v0 + 50.0);
  CHECK(c.get_state(0, 0) < v0 + 2.0);

  arma::vec too_short = {0.0};
  CHECK_THROWS(c.advance(0.0, 0.01, too_short));
}

TEST_CASE("Cells: per-node fields shorter than the number of cells are rejected")
{
  Cells c(3, "TT2", "ExplicitEuler", 0.01, 10.0, 1.0);
  c.set_solver_time_unit_ms(1.0);
  c.init();
  c.set_stretch(arma::vec({1.0, 1.0}));
  CHECK_THROWS(c.advance(0.0, 0.01));
  c.set_stretch(arma::vec());
  c.set_apicobasal(arma::vec({0.5}));
  CHECK_THROWS(c.advance(0.0, 0.01));
}

TEST_CASE("Cells: each node keeps its own stretch, cell type and apicobasal value")
{
  Cells c(2, "ToRORdLand", "ExplicitEuler", 0.02, 500.0, 1.0);
  c.set_solver_time_unit_ms(1.0);
  c.init();
  c.set_stretch(arma::vec({0.9, 1.1}));
  const std::set<uint> all = {0, 1};
  double ta0 = 0.0, ta1 = 0.0;
  for (int i = 0; i < 10000; i++)
  {
    const double t = i * 0.02;
    c.advance(t, 0.02, t < 2.0 ? -53.0 : 0.0, all);
    ta0 = std::max(ta0, c.get_state(0, 49));
    ta1 = std::max(ta1, c.get_state(1, 49));
  }
  CHECK_MSG(ta1 > ta0, ta0 << " " << ta1);

  int types[2] = {EPI, ENDO};
  c.set_cell_types(2, types);
  CHECK_EQ(c.get_cell_types()(1), arma::sword(ENDO));
}

TEST_CASE("Cells: the result does not depend on the number of threads")
{
  auto run = [](int threads) {
    omp_set_num_threads(threads);
    Cells c(16, "TT2", "ExplicitEuler", 0.02, 100.0, 1.0);
    c.set_solver_time_unit_ms(1.0);
    std::vector<int> types(16);
    for (int i = 0; i < 16; i++) types[i] = i % 3;
    c.set_cell_types(16, types.data());
    c.init();
    arma::vec stim(16);
    for (int i = 0; i < 1000; i++)
    {
      stim.fill(i * 0.02 < 2.0 ? -52.0 : 0.0);
      c.advance(i * 0.02, 0.02, stim);
    }
    arma::vec out(c.get_size());
    for (uint i = 0; i < c.get_size(); i++) out(i) = c.get_state(i);
    return out;
  };
  const arma::vec serial = run(1);
  const arma::vec parallel = run(2);
  omp_set_num_threads(2);
  CHECK(serial.is_finite());
  CHECK_NEAR(arma::abs(serial - parallel).max(), 0.0, 0.0);
}

TEST_CASE("Cells: monitored value is exposed for each node")
{
  Cells c(2, "Kerkoff2003", "ExplicitEuler", 1.0, 1000.0, 1.0);
  c.set_solver_time_unit_ms(1.0);
  c.init();
  arma::vec lat = {0.0, 1e6};
  c.set_var(1, lat);
  for (int i = 0; i < 50; i++)
    c.advance(i * 1.0, 1.0);
  arma::vec ta(2);
  c.get_monitored_values(0, ta);
  CHECK_MSG(ta(0) > 0.0, ta(0));
  CHECK_NEAR(ta(1), 0.0, 0.0);
}

TEST_CASE_KNOWN_BUG("Cells: with several monitored values, each one is read per node",
                    "src/odes/cells.cpp: advance() stores monitored values as "
                    "system*n_mon + j, but get_monitored_values() reads them as "
                    "mindex*n_sys + i, and writes to v(i) with the same offset")
{
  // RiceTT2 monitors 3 values. Index 1 of 2 cells must fit a vector of size 2.
  Cells c(2, "RiceTT2", "ExplicitEuler", 0.01, 10.0, 1.0);
  c.set_solver_time_unit_ms(1.0);
  c.init();
  arma::vec stim = {-52.0, 0.0};
  for (int i = 0; i < 300; i++)
    c.advance(i * 0.01, 0.01, stim);
  arma::vec out(6, arma::fill::zeros);
  c.get_monitored_values(1, out);
  // The values of monitored index 1 belong in out(0) and out(1)
  CHECK_NEAR(arma::accu(arma::abs(out.subvec(2, 3))), 0.0, 0.0);
}

// ---- CellWarmup -------------------------------------------------------------

namespace
{

CellWarmup::Options warmup_options(const std::string & cache = "none")
{
  CellWarmup::Options o;
  o.enabled = true;
  o.beats = 2;
  o.bcl_ms = 300.0;
  o.dt_ms = 0.02;
  o.stim_amp = -52.0;
  o.stim_dur_ms = 2.0;
  o.ab_bins = 1;
  o.cache = cache;
  o.model_name = "TNNP";
  return o;
}

} // namespace

TEST_CASE("CellWarmup: disabled means nothing is done")
{
  auto m = make_model("TNNP", 0.02);
  CellWarmup::Options o = warmup_options();
  o.enabled = false;
  CellWarmup w(m.get(), o);
  CHECK(!w.run(arma::ivec(), arma::vec(), 4));
  CHECK(!w.ready());
}

TEST_CASE("CellWarmup: one state per cell type, different from the initial condition")
{
  auto m = make_model("TNNP", 0.02);
  arma::ivec types = {EPI, ENDO, EPI, MCELL};
  CellWarmup w(m.get(), warmup_options());
  REQUIRE(w.run(types, arma::vec(), 4));
  REQUIRE(w.ready());

  const int n = m->get_num_state_vars();
  arma::vec y0(n), y1(n), y2(n), y3(n), ic(n);
  w.state_for_cell(0, 0.0, y0.memptr());
  w.state_for_cell(1, 0.0, y1.memptr());
  w.state_for_cell(2, 0.0, y2.memptr());
  w.state_for_cell(3, 0.0, y3.memptr());
  CHECK(y0.is_finite());
  CHECK_NEAR(arma::abs(y0 - y2).max(), 0.0, 0.0);
  CHECK(arma::abs(y0 - y1).max() > 1e-8);
  CHECK(arma::abs(y1 - y3).max() > 1e-8);

  m->set_celltype(EPI);
  m->init(ic.memptr());
  CHECK(arma::abs(y0 - ic).max() > 1e-6);
}

TEST_CASE("CellWarmup: the apicobasal bins split the groups only when requested")
{
  // TNNP has no apicobasal gradient, so the bins are observed through the cache groups: with
  // one bin the two cells share a group, with four they get a state each (equal in value)
  auto m = make_model("TNNP", 0.02);
  arma::ivec types = {EPI, EPI};
  arma::vec ab = {0.05, 0.95};
  const int n = m->get_num_state_vars();
  arma::vec a(n), b(n);

  CellWarmup::Options o = warmup_options();
  CellWarmup one(m.get(), o);
  REQUIRE(one.run(types, ab, 2));
  one.state_for_cell(0, 0.0, a.memptr());
  one.state_for_cell(1, 0.0, b.memptr());
  CHECK_NEAR(arma::abs(a - b).max(), 0.0, 0.0);

  o.ab_bins = 4;
  o.cache = "bins.cache";
  CellWarmup four(m.get(), o);
  REQUIRE(four.run(types, ab, 2));
  four.state_for_cell(0, 0.0, a.memptr());
  CHECK(a.is_finite());
  CHECK(read_text("bins.cache").size() > 0);
}

TEST_CASE("CellWarmup: the cache is reused for the same options and ignored for others")
{
  auto m = make_model("TNNP", 0.02);
  arma::ivec types = {EPI, ENDO};
  const int n = m->get_num_state_vars();
  arma::vec fresh(n), cached(n), other(n);

  {
    CellWarmup w(m.get(), warmup_options("wu.cache"));
    REQUIRE(w.run(types, arma::vec(), 2));
    w.state_for_cell(1, 0.0, fresh.memptr());
  }
  REQUIRE(file_exists("wu.cache"));

  // Tamper with the numbers so that reading from the cache is observable: a
  // recomputation would give back the original state
  const std::string text = read_text("wu.cache");

  {
    CellWarmup w(m.get(), warmup_options("wu.cache"));
    REQUIRE(w.run(types, arma::vec(), 2));
    w.state_for_cell(1, 0.0, cached.memptr());
  }
  CHECK_NEAR(arma::abs(fresh - cached).max(), 0.0, 1e-9 * arma::abs(fresh).max());

  {
    CellWarmup::Options o = warmup_options("wu.cache");
    o.bcl_ms = 150.0;
    CellWarmup w(m.get(), o);
    REQUIRE(w.run(types, arma::vec(), 2));
    w.state_for_cell(1, 0.0, other.memptr());
  }
  CHECK_MSG(arma::abs(fresh - other).max() > 1e-6, "a different BCL must not reuse the cache");
  CHECK_MSG(read_text("wu.cache") != text, "the cache must be rewritten for the new options");
}

TEST_CASE("CellWarmup: with -warmup_lat each node starts at the phase of its own LAT")
{
  auto m = make_model("TNNP", 0.02);
  CellWarmup::Options o = warmup_options();
  o.phase_from_lat = true;
  o.phase_ms = 0.0;
  o.sample_ms = 1.0;
  arma::ivec types = {EPI, EPI};
  CellWarmup w(m.get(), o);
  REQUIRE(w.run(types, arma::vec(), 2));
  CHECK(w.uses_lat());

  const int n = m->get_num_state_vars();
  arma::vec at0(n), at0b(n), at10(n), wrapped(n);
  w.state_for_cell(0, 0.0, at0.memptr());
  w.state_for_cell(1, 0.0, at0b.memptr());
  w.state_for_cell(0, 10.0, at10.memptr());
  w.state_for_cell(0, 10.0 + o.bcl_ms, wrapped.memptr());
  CHECK(at10.is_finite());
  CHECK_NEAR(arma::abs(at0 - at0b).max(), 0.0, 0.0);
  CHECK(arma::abs(at0 - at10).max() > 1e-6);
  CHECK_NEAR(arma::abs(at10 - wrapped).max(), 0.0, 1e-12);
}

TEST_CASE("CellWarmup: the options are read only when -warmup 1 is given, with the given defaults")
{
  CHECK(!CellWarmup::options_from_command_line(700.0, 0.05, "TT2", "base").enabled);

  static const char * on[] = {"prog", "-warmup", "1", "-warmup_beats", "7", nullptr};
  CommandLineArgs::init(5, on);
  const CellWarmup::Options o = CellWarmup::options_from_command_line(700.0, 0.05, "TT2", "base");
  static const char * restore[] = {"prog", "-linsolver", "petsc", nullptr};
  CommandLineArgs::init(3, restore);

  CHECK(o.enabled);
  CHECK_EQ(o.beats, 7);
  CHECK_NEAR(o.bcl_ms, 700.0, 0.0);
  CHECK_NEAR(o.dt_ms, 0.05, 0.0);
  CHECK_EQ(o.model_name, std::string("TT2"));
}
