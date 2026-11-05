/*
 * Unit tests of src/cardiac, the parts that can run without a full 3D simulation.
 *
 * What is covered and why:
 *  - ConductivityModel: the three anisotropy types, regions and name conversion.
 *  - Stimuli: reading from TOML, and which nodes receive what at which time.
 *  - Monodomain element and global matrices: through a subclass that exposes them. The
 *    stiffness carries the conductivity tensor, and the global matrices must combine
 *    mass and stiffness with the theta weight.
 *  - Regazzoni2020 (0D closed loop, no 3D model): conservation of blood volume, a
 *    limit cycle with physiological pressures and valve behaviour, the timing rules
 *    (RR scaling, atrial anchoring, clamping), and the 3D coupling hook.
 *  - PiecewisePVLoop (port of pulse.cycle) against a lumped cavity whose constraints are
 *    solved exactly: the five phases in order and at the prescribed times, the preload
 *    ramp, the 3-element Windkessel (affine law and drain), the valve rules, the exact
 *    volume of the isovolumic phases, stroke volume and ejection fraction.
 *  - CirculationModel implementations (Regazzoni, Piecewise, Legacy) driven through
 *    the CavityMechanics interface by a lumped mock, so that the couplings run in
 *    milliseconds: the pressure inversion reproduces the 0D volume, the -t ceiling
 *    reduces the beats, the run plan is right, and save_state/restore_state continue a
 *    run exactly as if it had not stopped.
 *  - cell_coupling: profiles per cell model, the unknown model error, and TaCurve.
 *  - LandCoupling: nodal active stiffness from the cell states.
 *  - UVC data transfer (kdtree_uvc): identical meshes reproduce the field, the values
 *    stay within the source range, categorical fields keep their labels.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <memory>

#include "cardiac/conductivity.hpp"
#include "cardiac/stimulus.hpp"
#include "cardiac/monodomain.hpp"
#include "cardiac/regazzoni2020.hpp"
#include "cardiac/piecewise_pvloop.hpp"
#include "cardiac/circulation_model.hpp"
#include "cardiac/regazzoni_circulation.hpp"
#include "cardiac/piecewise_circulation.hpp"
#include "cardiac/legacy_circulation.hpp"
#include "cardiac/cell_coupling.hpp"
#include "odes/torord_land.hpp"
#include "cardiac/uvc/kdtree_uvc.hpp"

using namespace cardiax_test;

// ---- ConductivityModel ------------------------------------------------------

namespace
{
const arma::vec3 kF = {1, 0, 0}, kS = {0, 1, 0}, kN = {0, 0, 1};
}

TEST_CASE("ConductivityModel: isotropic uses sigma_l on the diagonal")
{
  ConductivityModel c;
  c.set_type(ISOTROPIC);
  c.set_default({2.0, 3.0, 4.0});
  const arma::mat33 s = c.at(0, 3, kF, kS, kN);
  CHECK_NEAR(s(0, 0), 2.0, 0.0);
  CHECK_NEAR(s(2, 2), 2.0, 0.0);
  CHECK_NEAR(s(0, 1), 0.0, 0.0);
}

TEST_CASE("ConductivityModel: transversely isotropic with fibres along an oblique direction")
{
  ConductivityModel c;
  c.set_type(TRANSVERSELY_ISOTROPIC);
  c.set_default({5.0, 1.0, 99.0});
  const arma::vec3 f = {std::sqrt(0.5), std::sqrt(0.5), 0.0};
  const arma::mat33 s = c.at(0, 3, f, kS, kN);
  CHECK_NEAR(arma::dot(f, s * f), 5.0, 1e-14);
  const arma::vec3 t = {-std::sqrt(0.5), std::sqrt(0.5), 0.0};
  CHECK_NEAR(arma::dot(t, s * t), 1.0, 1e-14);
  CHECK_NEAR(s(2, 2), 1.0, 1e-14);
}

TEST_CASE("ConductivityModel: orthotropic recovers each conductivity along its axis")
{
  ConductivityModel c;
  c.set_type(ORTHOTROPIC);
  c.set_default({5.0, 2.0, 1.0});
  const double a = 0.6;
  const arma::vec3 f = {std::cos(a), std::sin(a), 0}, s_ = {-std::sin(a), std::cos(a), 0};
  const arma::mat33 s = c.at(0, 3, f, s_, kN);
  CHECK_NEAR(arma::dot(f, s * f), 5.0, 1e-14);
  CHECK_NEAR(arma::dot(s_, s * s_), 2.0, 1e-14);
  CHECK_NEAR(arma::dot(kN, s * kN), 1.0, 1e-14);
  CHECK_NEAR(arma::dot(f, s * s_), 0.0, 1e-14);
}

TEST_CASE("ConductivityModel: a region overrides the default")
{
  ConductivityModel c;
  c.set_default({1.0, 1.0, 1.0});
  c.set_region(7, {9.0, 9.0, 9.0});
  CHECK_NEAR(c.values(7).sigma_l, 9.0, 0.0);
  CHECK_NEAR(c.values(3).sigma_l, 1.0, 0.0);
  CHECK_NEAR(c.at(7, 3, kF, kS, kN)(1, 1), 9.0, 0.0);
}

TEST_CASE("ConductivityModel: names convert both ways, and an unknown name throws")
{
  CHECK(ConductivityModel::from_string("isotropic") == ISOTROPIC);
  CHECK(ConductivityModel::from_string("transversely_isotropic") == TRANSVERSELY_ISOTROPIC);
  CHECK(ConductivityModel::from_string("orthotropic") == ORTHOTROPIC);
  CHECK_THROWS(ConductivityModel::from_string("nope"));
  CHECK_EQ(ConductivityModel::to_string(ORTHOTROPIC), std::string("Orthotropic"));
}

// ---- Stimuli ----------------------------------------------------------------

namespace
{

toml::table stim_cfg()
{
  return toml::parse(
    "[stimuli]\n"
    "tissue = [ { start = 1.0, duration = 2.0, value = -50.0, min = [0, 0, 0], max = [0.1, 1, 1] },\n"
    "           { start = 10.0, duration = 1.0, value = -20.0, min = [0.9, 0, 0], max = [1, 1, 1] } ]\n");
}

} // namespace

TEST_CASE("Stimuli: TOML entries select the nodes in the box during the interval")
{
  BoxSpec b;
  b.nx = 4; b.ny = 1; b.nz = 1;
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");

  Stimuli st;
  CHECK_EQ(st.read_toml(stim_cfg(), "stimuli.tissue"), 2);

  std::set<uint> nodes;
  double value = 0.0;
  bool apply = true;

  st.check(0.5, m, nodes, &value, &apply);
  CHECK(!apply);

  st.check(2.0, m, nodes, &value, &apply);
  CHECK(apply);
  CHECK_NEAR(value, -50.0, 0.0);
  CHECK_EQ(nodes.size(), std::size_t(4));          // the 4 nodes at x = 0
  CHECK_EQ(nodes.count(box_node(b, 0, 1, 1)), std::size_t(1));
  CHECK_EQ(nodes.count(box_node(b, 1, 0, 0)), std::size_t(0));

  nodes.clear();
  st.check(5.0, m, nodes, &value, &apply);
  CHECK(!apply);

  nodes.clear();
  st.check(10.5, m, nodes, &value, &apply);
  CHECK(apply);
  CHECK_NEAR(value, -20.0, 0.0);
  CHECK_EQ(nodes.count(box_node(b, 4, 0, 0)), std::size_t(1));
}

TEST_CASE("Stimuli: an absent section gives zero stimuli, and a malformed entry throws")
{
  Stimuli st;
  toml::table empty;
  CHECK_EQ(st.read_toml(empty, "stimuli.tissue"), 0);

  CHECK_THROWS(st.read_toml(toml::parse(
    "[stimuli]\ntissue = [ { start = 1.0, value = -50.0, min = [0, 0, 0], max = [1, 1, 1] } ]\n"),
    "stimuli.tissue"));
  CHECK_THROWS(st.read_toml(toml::parse(
    "[stimuli]\ntissue = [ { start = 1.0, duration = 1.0, value = 1.0, min = [0, 0], max = [1, 1, 1] } ]\n"),
    "stimuli.tissue"));
}

TEST_CASE("Stimuli: added directly, with the same behaviour")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 1;
  write_tet_box("m.xml", b);
  Mesh m;
  m.read_xml("m.xml");
  Stimuli st;
  st.add(0.0, 1.0, -7.0, -1, 2, -1, 2, -1, 2);
  std::set<uint> nodes;
  double v = 0.0;
  bool apply = false;
  st.check(0.5, m, nodes, &v, &apply);
  CHECK(apply);
  CHECK_EQ(nodes.size(), std::size_t(8));
  CHECK_NEAR(v, -7.0, 0.0);
}

// ---- Monodomain matrices ----------------------------------------------------

namespace
{

struct OpenMonodomain : public Monodomain
{
  using Monodomain::calc_elmat_stiff;
  using Monodomain::calc_elmat_mass;
  using Monodomain::fespace;
  using Monodomain::Ai;
  using Monodomain::Mi;
  using Monodomain::assemble_matrices;
  using Monodomain::mesh;
  using Monodomain::vm;
  using CardiacProblem::parameters;
};

} // namespace

TEST_CASE("Monodomain: the element stiffness carries the anisotropic conductivity")
{
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  write_tet_box("slab.xml", b);

  std::string mesh = "slab.xml", cell = "NP", ode = "ExplicitEuler";
  OpenMonodomain md;
  md.setup(mesh, cell, ode, 0.1, 1.0, 1.0, 1.0);
  md.set_parameters(toml::parse(
    "[physical]\nconductivity_type = \"orthotropic\"\n"
    "sigma_f = 3.0\nsigma_s = 2.0\nsigma_n = 1.0\nsurface_to_volume = 0.14\n"
    "[numerics]\ntheta = 0.5\nmaxnz = 40\n"));
  md.init();

  FiniteElement & fe = md.fespace.createFE(0);
  const int n = fe.get_ndof();
  arma::mat Ke(n, n), Me(n, n);

  // u = x on the element: the energy is sigma_xx * volume, so it picks sigma_f
  auto energy = [&](int e, int d) {
    md.calc_elmat_stiff(e, fe, Ke);
    std::vector<int> pn;
    md.mesh->get_element_pt_nums(e, pn);
    arma::vec u(n);
    for (int a = 0; a < n; a++) u(a) = md.mesh->get_point(pn[a])(d);
    return arma::dot(u, Ke * u);
  };
  md.calc_elmat_mass(0, fe, Me);
  const double vol = arma::accu(Me);
  CHECK_NEAR(energy(0, 0) / vol, 3.0, 1e-12);
  CHECK_NEAR(energy(0, 1) / vol, 2.0, 1e-12);
  CHECK_NEAR(energy(0, 2) / vol, 1.0, 1e-12);

  md.calc_elmat_stiff(0, fe, Ke);
  CHECK_NEAR(arma::abs(Ke - Ke.t()).max(), 0.0, 1e-13);
  CHECK_NEAR(arma::abs(Ke * arma::ones<arma::vec>(n)).max(), 0.0, 1e-13);
  delete &fe;
}

TEST_CASE("Monodomain: a constant Vm is a fixed point of the diffusion matrices")
{
  // A*1 and M*1 are both the lumped mass: diffusion alone does not change a uniform Vm
  BoxSpec b;
  b.nx = b.ny = b.nz = 2;
  write_tet_box("slab2.xml", b);

  std::string mesh = "slab2.xml", cell = "NP", ode = "ExplicitEuler";
  OpenMonodomain md;
  md.setup(mesh, cell, ode, 0.1, 1.0, 1.0, 1.0);
  md.set_parameters(toml::parse(
    "[physical]\nconductivity_type = \"isotropic\"\nsigma_f = 1.0\nsigma_s = 1.0\nsigma_n = 1.0\n"));
  md.init();

  const int n = md.mesh->get_n_points();
  petsc::Vector one, a, m;
  one.create(n);
  a.create(n);
  m.create(n);
  one.set(1.0);
  md.Ai.mult(one, a);
  md.Mi.mult(one, m);

  arma::vec va(n), vm_(n);
  a.get_data(va.memptr());
  m.get_data(vm_.memptr());
  // The stiffness is multiplied by dt/kappa ~ 1e8 here, which scales its roundoff as well
  CHECK_NEAR(arma::abs(va - vm_).max(), 0.0, 1e-6);
  CHECK_NEAR(arma::accu(va), 1.0, 1e-6);   // the unit cube
}

// ---- Regazzoni2020 ----------------------------------------------------------

namespace
{
const double kMmHg = 133.322387415;

struct BeatStats { double vtot_gap = 0, plv_max = 0, plv_min = 1e30, vlv_max = 0, vlv_min = 1e30; };

BeatStats run_0d(Regazzoni2020 & c, double t0, int steps, double dt = 1.0)
{
  BeatStats s;
  const double v0 = c.total_volume();
  for (int i = 0; i < steps; i++)
  {
    c.step(t0 + i * dt, dt);
    s.plv_max = std::max(s.plv_max, c.pressure_LV());
    s.plv_min = std::min(s.plv_min, c.pressure_LV());
    s.vlv_max = std::max(s.vlv_max, c[Regazzoni2020::V_LV]);
    s.vlv_min = std::min(s.vlv_min, c[Regazzoni2020::V_LV]);
    s.vtot_gap = std::max(s.vtot_gap, std::fabs(c.total_volume() - v0));
  }
  return s;
}

} // namespace

TEST_CASE("Regazzoni2020: blood volume is conserved over several beats")
{
  Regazzoni2020 c;
  const BeatStats s = run_0d(c, 0.0, 5000);
  CHECK_MSG(s.vtot_gap < 1e-8 * c.total_volume(), "drift " << s.vtot_gap << " mL");
}

TEST_CASE("Regazzoni2020: reaches a periodic limit cycle with physiological values")
{
  Regazzoni2020 c;
  const int beat = static_cast<int>(c.period());
  double t = 0.0;
  for (int b = 0; b < 30; b++, t += c.period())
    run_0d(c, t, beat);

  std::array<double, Regazzoni2020::NSTATE> y0;
  for (int i = 0; i < Regazzoni2020::NSTATE; i++) y0[i] = c[i];

  const BeatStats s = run_0d(c, t, beat);

  for (int i = 0; i < 4; i++)
    CHECK_MSG(std::fabs(c[i] - y0[i]) < 0.5, "state " << i << ": " << y0[i] << " -> " << c[i]);

  CHECK_MSG(s.plv_max / kMmHg > 70.0 && s.plv_max / kMmHg < 200.0, "peak LV " << s.plv_max / kMmHg);
  CHECK_MSG(s.plv_min / kMmHg < 25.0, "min LV " << s.plv_min / kMmHg);
  const double sv = s.vlv_max - s.vlv_min;
  CHECK_MSG(sv > 20.0 && sv < 150.0, "stroke volume " << sv);
  CHECK_GOLDEN("regazzoni0d.plv_max", s.plv_max, 1e-6);
  CHECK_GOLDEN("regazzoni0d.stroke_volume", sv, 1e-6);
}

TEST_CASE("Regazzoni2020: the valves are diodes, with a small reverse leak at most")
{
  Regazzoni2020 c;
  double t = 0.0;
  const int beat = static_cast<int>(c.period());
  for (int b = 0; b < 5; b++, t += c.period()) run_0d(c, t, beat);
  for (int i = 0; i < beat; i++)
  {
    c.step(t + i, 1.0);
    if (c.pressure_LA() > c.pressure_LV() + 1.0) CHECK(c.flow_MV() > 0.0);
    if (c.pressure_LV() > c.pressure_LA() + 1.0) CHECK(c.flow_MV() < 0.0);
    // Closed: resistance much larger than the forward one, so a much smaller flow
    if (c.pressure_LV() > c.pressure_LA() + 133.0)
      CHECK_MSG(std::fabs(c.flow_MV()) < 0.01, c.flow_MV());
  }
}

TEST_CASE("Regazzoni2020: the period regenerates the chamber timing")
{
  Regazzoni2020 c;
  const double diastole_1000 = c.diastolic_window();
  CHECK_NEAR(c.period(), 1000.0, 0.0);
  CHECK(!c.timing_was_clamped());

  c.set_period(600.0);
  CHECK_NEAR(c.period(), 600.0, 0.0);
  CHECK_MSG(c.diastolic_window() < diastole_1000, "the shorter period comes out of diastole");

  // Periods too short for the atria: the atrial systole is cut, with a flag
  Regazzoni2020 tight;
  tight.set_period(200.0);
  CHECK(tight.timing_was_clamped() || tight.diastolic_window() <= 0.0);
  CHECK(!tight.describe_timing().empty());
}

TEST_CASE("Regazzoni2020: activation time and EMD change the timing")
{
  Regazzoni2020 a, b;
  b.set_activation_time(100.0);
  CHECK_NEAR(b.timing().t_act, 100.0, 0.0);
  b.set_emd(10.0);
  CHECK_NEAR(b.timing().emd, 10.0, 0.0);

  // Delaying activation delays the pressure rise
  const auto rise_time = [](Regazzoni2020 & c) {
    double t_peak = 0.0, p_peak = 0.0;
    for (int b = 0; b < 3; b++)
      for (int i = 0; i < 1000; i++)
      {
        c.step(b * 1000.0 + i, 1.0);
        if (b == 2 && c.pressure_LV() > p_peak) { p_peak = c.pressure_LV(); t_peak = i; }
      }
    return t_peak;
  };
  CHECK(rise_time(b) > rise_time(a));
}

TEST_CASE("Regazzoni2020: rebalancing restores the total volume through the venous reservoirs")
{
  Regazzoni2020 c;
  const double vtot = c.total_volume();
  const double p_ven0 = c[Regazzoni2020::p_VEN_SYS];
  c.set_initial_volumes(c[Regazzoni2020::V_LV] - 30.0, c[Regazzoni2020::V_RV] - 40.0);
  CHECK_NEAR(c.total_volume(), vtot - 70.0, 1e-9);

  const double dv = c.rebalance_total_volume(vtot);
  CHECK_NEAR(dv, 70.0, 1e-9);
  CHECK_NEAR(c.total_volume(), vtot, 1e-9);
  CHECK(c[Regazzoni2020::p_VEN_SYS] > p_ven0);
}

TEST_CASE("Regazzoni2020: with 3D coupling, the pressures come from the callback")
{
  Regazzoni2020 c;
  double seen_vlv = -1.0, seen_t = -1.0;
  int calls = 0;
  c.set_BiV_coupling([&](double vlv, double, double t, double & plv, double & prv) {
    seen_vlv = vlv;
    seen_t = t;
    plv = 1234.0;
    prv = 567.0;
    calls++;
  });
  const double v_expected = c[Regazzoni2020::V_LV];
  c.step(0.0, 1.0);
  CHECK(calls >= 1);
  CHECK_NEAR(seen_vlv, v_expected, 1.0);
  CHECK(seen_t >= 0.0);
  CHECK_NEAR(c.pressure_LV(), 1234.0, 0.0);
  CHECK_NEAR(c.pressure_RV(), 567.0, 0.0);
}

TEST_CASE("Regazzoni2020: a history row has one value per column of the header")
{
  Regazzoni2020 c;
  c.step(0.0, 1.0);
  std::istringstream head(Regazzoni2020::history_header()), row(c.history_row(1.0));
  std::string tok;
  int n_head = 0, n_row = 0;
  while (head >> tok) if (tok != "#") n_head++;
  while (row >> tok) n_row++;
  CHECK_EQ(n_head, n_row);
}

// ---- Lumped cavity mock -----------------------------------------------------

namespace
{

//! Time-varying compliance in series: V = V0 + C (p - Pact(t)), with an active
//! pressure that peaks at t_act + 100 ms in each period. It behaves enough like a
//! ventricle for the phase logic and the Newton of the 0D inversion to work on it.
struct LumpedCavity : public CavityMechanics
{
  double V0_lv = 80.0, V0_rv = 90.0;
  double C = 0.002;                //!< mL/Pa
  double Pmax = 16000.0;           //!< Pa
  double period = 1000.0, t_act = 0.0, T_sys = 300.0;
  double vlv, vrv;
  int solves = 0;
  double last_t = 0.0;

  LumpedCavity() : vlv(V0_lv), vrv(V0_rv) {}

  double active(double t) const
  {
    double s = std::fmod(t - t_act, period);
    if (s < 0.0) s += period;
    return (s < T_sys) ? Pmax * std::pow(std::sin(M_PI * s / T_sys), 2) : 0.0;
  }

  void apply_pressures(double t, double plv, double prv) override
  {
    last_t = t;
    vlv = V0_lv + C * (plv - active(t));
    vrv = V0_rv + C * (prv - 0.3 * active(t));
    solves++;
  }
  double volume_LV() override { return vlv; }
  double volume_RV() override { return vrv; }
  double volume_total() override { return vlv + vrv; }
  bool verbose() const override { return false; }
};

// Runs the loop as CardiacElectromechanic does, from step `first` to `last`
void run_steps(CirculationModel & m, int first, int last)
{
  const RunPlan & plan = m.plan();
  for (int i = first; i < last; i++)
  {
    const StepContext s = {i, plan.nsteps, (i + 1) * plan.step_ms, plan.step_ms, false};
    m.announce(s);
    m.advance(s);
  }
}

} // namespace

// ---- PiecewisePVLoop --------------------------------------------------------

namespace
{

typedef PiecewisePVLoop::Constraint PwConstraint;

// Pressure that satisfies a constraint exactly on the lumped cavity
// V = V0 + C (p - a), which is linear in p
double lumped_pressure(const PwConstraint & c, double V0, double C, double a)
{
  switch (c.kind)
  {
    case PwConstraint::PRESSURE: return c.P;
    case PwConstraint::VOLUME:   return (c.V - V0) / C + a;
    case PwConstraint::AFFINE:   return (c.A + c.B * (V0 - C * a)) / (1.0 - c.B * C);
  }
  return std::nan("");
}

// One step of the controller on the lumped cavity, with the constraints solved exactly
// (what the monolithic Newton of pulse does)
void lumped_step(PiecewisePVLoop & pw, LumpedCavity & mech, double t, double dt)
{
  const PwConstraint cl = pw.constraint_LV(t, dt), cr = pw.constraint_RV(t, dt);
  const double a = mech.active(t);
  const double pl = lumped_pressure(cl, mech.V0_lv, mech.C, a);
  const double pr = lumped_pressure(cr, mech.V0_rv, mech.C, 0.3 * a);
  mech.apply_pressures(t, pl, pr);
  pw.commit(t, dt, mech.volume_LV(), pl, mech.volume_RV(), pr);
}

} // namespace

TEST_CASE("PiecewisePVLoop: phase names, defaults of complete_cycle.py, activation and run length")
{
  PiecewisePVLoop pw;
  CHECK_EQ(std::string(PiecewisePVLoop::phase_name(PiecewisePVLoop::EJECTION)).empty(), false);

  // Numbering of pulse.cycle.Phase
  CHECK_EQ(int(PiecewisePVLoop::PRELOAD), 0);
  CHECK_EQ(int(PiecewisePVLoop::IVC), 1);
  CHECK_EQ(int(PiecewisePVLoop::EJECTION), 2);
  CHECK_EQ(int(PiecewisePVLoop::IVR), 3);
  CHECK_EQ(int(PiecewisePVLoop::FILLING), 4);

  CHECK_NEAR(pw.period(), 800.0, 0.0);
  CHECK_NEAR(pw.lv().t_zero, 50.0, 0.0);
  CHECK_NEAR(pw.lv().windkessel.p_init, 9000.0, 0.0);
  CHECK_NEAR(pw.rv().windkessel.p_init, 2500.0, 0.0);
  // 1.1 mmHg*s/mL and 1.5 mL/mmHg in Pa*ms/mL and mL/Pa
  CHECK_NEAR(pw.lv().windkessel.R, 1.1 * 133.322387415 * 1000.0, 1e-9);
  CHECK_NEAR(pw.lv().windkessel.C, 1.5 / 133.322387415, 1e-15);

  pw.lv().t_end_diastole = 150.0;
  CHECK_NEAR(pw.activation_time(), 150.0, 0.0);
  CHECK_NEAR(pw.total_time(3), 2400.0, 0.0);
  CHECK_EQ(pw.beat_of(1.0), 0);
  CHECK_EQ(pw.beat_of(800.0), 0);   // the step that ends at the period is still beat 0
  CHECK_EQ(pw.beat_of(801.0), 1);
}

TEST_CASE("PiecewisePVLoop: the preload ramp, 0 -> preload_pressure then -> p_end_diastole")
{
  PiecewisePVLoop pw;
  pw.initialize(0.0, 80.0, 0.0, 90.0, 0.0);
  const PiecewisePVLoop::ChamberParams & p = pw.lv();
  auto P = [&](double t) {
    const PwConstraint c = pw.constraint_LV(t, 1.0);
    CHECK(c.kind == PwConstraint::PRESSURE);
    return c.P;
  };
  CHECK_NEAR(P(0.5 * p.t_zero), 0.5 * p.preload_pressure, 1e-9);
  CHECK_NEAR(P(p.t_zero), p.preload_pressure, 1e-9);
  CHECK_NEAR(P(p.t_end_diastole), p.p_end_diastole, 1e-9);
  CHECK_NEAR(P(0.5 * (p.t_zero + p.t_end_diastole)),
             0.5 * (p.preload_pressure + p.p_end_diastole), 1e-9);

  // In the next beats the ramp starts from preload_pressure, not from 0
  PiecewisePVLoop::ChamberState s = pw.state_LV();
  s.n_beats = 1;
  pw.set_states(s, pw.state_RV());
  CHECK_NEAR(P(p.period + p.t_zero), p.preload_pressure, 1e-9);
  CHECK_NEAR(P(p.period + p.t_end_diastole), p.p_end_diastole, 1e-9);
}

TEST_CASE("PiecewisePVLoop: the 3-element Windkessel, ejection law and drain")
{
  PiecewisePVLoop pw;
  const PiecewisePVLoop::Windkessel & w = pw.lv().windkessel;
  const double P_c = 10000.0, V_n = 100.0, dt = 1.0, V = 99.5;

  // The pressure imposed by the affine law is the one the Windkessel ends at
  double A, B;
  w.ejection_law(P_c, V_n, dt, A, B);
  CHECK(B < 0.0);
  double P_c1, Q;
  w.advance(P_c, V, V_n, dt, true, true, P_c1, Q);
  CHECK_NEAR(Q, 0.5, 1e-12);
  CHECK_NEAR(A + B * V, P_c1 + w.Rc * Q, 1e-9 * P_c);
  // Implicit Euler of C dP_c/dt = Q - P_c/R
  CHECK_NEAR(w.C * (P_c1 - P_c) / dt, Q - P_c1 / w.R, 1e-12);

  // Valve closed: P_c is frozen until the first ejection, then drains through R
  w.advance(P_c, V, V_n, dt, false, false, P_c1, Q);
  CHECK_NEAR(P_c1, P_c, 0.0);
  CHECK_NEAR(Q, 0.0, 0.0);
  w.advance(P_c, V, V_n, dt, false, true, P_c1, Q);
  CHECK_NEAR(P_c1, P_c / (1.0 + dt / (w.R * w.C)), 1e-9);
}

TEST_CASE("PiecewisePVLoop: goes through the five phases in order on a lumped cavity")
{
  LumpedCavity mech;
  PiecewisePVLoop pw;
  mech.period = pw.period();
  mech.t_act = pw.activation_time();
  pw.initialize(0.0, mech.volume_LV(), 0.0, mech.volume_RV(), 0.0);
  CHECK(pw.state_LV().phase == PiecewisePVLoop::PRELOAD);

  std::vector<int> order_lv, order_rv;
  double p_max = 0.0, t_ivc = -1.0, t_back = -1.0, v_gap_ivc = 0.0, v_gap_ivr = 0.0;

  const double dt = 1.0;
  const int nsteps = static_cast<int>((pw.period() + pw.lv().t_zero + 10.0) / dt);
  for (int i = 0; i < nsteps; i++)
  {
    const double t = (i + 1) * dt;
    const int before = pw.state_LV().phase;
    lumped_step(pw, mech, t, dt);
    const PiecewisePVLoop::ChamberState & s = pw.state_LV();

    if (order_lv.empty() || order_lv.back() != s.phase) order_lv.push_back(s.phase);
    const int ph_rv = pw.state_RV().phase;
    if (order_rv.empty() || order_rv.back() != ph_rv) order_rv.push_back(ph_rv);
    p_max = std::max(p_max, s.pressure_n);

    // The step was solved under `before`
    if (before == PiecewisePVLoop::IVC)
      v_gap_ivc = std::max(v_gap_ivc, std::fabs(s.volume_n - s.end_dia_vol));
    if (before == PiecewisePVLoop::IVR)
      v_gap_ivr = std::max(v_gap_ivr, std::fabs(s.volume_n - s.end_sys_vol));
    if (s.phase == PiecewisePVLoop::IVC && t_ivc < 0.0) t_ivc = t;
    if (before == PiecewisePVLoop::FILLING && s.phase == PiecewisePVLoop::PRELOAD) t_back = t;
  }

  const std::vector<int> expected = {0, 1, 2, 3, 4, 0};
  auto show = [](const std::vector<int> & v) {
    std::ostringstream os;
    for (int p : v) os << p << " ";
    return os.str();
  };
  CHECK_MSG(order_lv == expected, show(order_lv));
  CHECK_MSG(order_rv == expected, show(order_rv));

  // Phase changes at the prescribed times
  CHECK_NEAR(t_ivc, pw.lv().t_end_diastole, 0.0);
  CHECK_NEAR(t_back, pw.period() + pw.lv().t_zero, 0.0);
  CHECK_EQ(pw.state_LV().n_beats, 1);

  // Isovolumic phases keep the volume exactly (the constraint is solved, not penalised)
  CHECK_MSG(v_gap_ivc < 1e-9, "volume moved " << v_gap_ivc << " mL in IVC");
  CHECK_MSG(v_gap_ivr < 1e-9, "volume moved " << v_gap_ivr << " mL in IVR");
  CHECK(p_max > pw.lv().windkessel.p_init);

  const double sv = pw.stroke_volume_LV();
  CHECK_MSG(sv > 0.0, "stroke volume " << sv);
  CHECK(pw.ejection_fraction_LV() > 0.0);
  CHECK(pw.ejection_fraction_LV() < 100.0);
  CHECK_NEAR(pw.ejection_fraction_LV(), 100.0 * sv / pw.state_LV().end_dia_vol, 1e-12);
  CHECK(pw.stroke_volume_RV() > 0.0);
}

TEST_CASE("PiecewisePVLoop: the valve opens when P > P_c, and stays open min_ejection_duration")
{
  PiecewisePVLoop pw;
  pw.initialize(0.0, 100.0, 0.0, 100.0, 0.0);
  PiecewisePVLoop::ChamberState s = pw.state_LV();
  s.phase = PiecewisePVLoop::IVC;
  s.end_dia_vol = 100.0;
  pw.set_states(s, pw.state_RV());
  const double P_c = pw.state_LV().wdk_pressure_n;

  pw.commit(200.0, 1.0, 100.0, P_c - 1.0, 100.0, 0.0);
  CHECK(pw.state_LV().phase == PiecewisePVLoop::IVC);
  pw.commit(201.0, 1.0, 100.0, P_c + 1.0, 100.0, 0.0);
  CHECK(pw.state_LV().phase == PiecewisePVLoop::EJECTION);
  CHECK(pw.state_LV().has_ejected);

  // No outflow at all, but the ejection lasts min_ejection_duration before closing
  const double t_open = 201.0, d = pw.lv().min_ejection_duration;
  double t = t_open;
  while (pw.state_LV().phase == PiecewisePVLoop::EJECTION && t < t_open + 10.0 * d)
  {
    t += 1.0;
    pw.commit(t, 1.0, 100.0, 2.0 * P_c, 100.0, 0.0);
  }
  CHECK(pw.state_LV().phase == PiecewisePVLoop::IVR);
  CHECK_MSG(t > t_open + d && t <= t_open + d + 1.0, "closed at " << t);
  CHECK_NEAR(pw.state_LV().end_sys_vol, 100.0, 0.0);
}

TEST_CASE("PiecewisePVLoop: configure reads the XML block, and -rr wins over it")
{
  write_text("pw.xml",
    "<mesh celltype=\"tetrahedron\" dim=\"3\">"
    "<pvloop_piecewise period=\"700\" t_zero=\"40\" t_end_diastole=\"100\">"
    "<lv p_init=\"8000\" fill_rate=\"0.2\"/><rv p_fill=\"150\"/>"
    "</pvloop_piecewise></mesh>");
  PiecewisePVLoop a;
  a.configure("pw.xml");
  CHECK_NEAR(a.period(), 700.0, 0.0);
  CHECK_NEAR(a.lv().t_zero, 40.0, 0.0);
  CHECK_NEAR(a.rv().t_end_diastole, 100.0, 0.0);
  CHECK_NEAR(a.lv().windkessel.p_init, 8000.0, 0.0);
  CHECK_NEAR(a.lv().fill_rate, 0.2, 0.0);
  CHECK_NEAR(a.rv().p_fill, 150.0, 0.0);
  CHECK_NEAR(a.rv().windkessel.p_init, 2500.0, 0.0);   // default kept
  CHECK(!a.describe().empty());

  PiecewisePVLoop b;
  b.configure("pw.xml", 500.0);
  CHECK_NEAR(b.period(), 500.0, 0.0);
  CHECK_NEAR(b.rv().period, 500.0, 0.0);

  std::istringstream head(PiecewisePVLoop::history_header());
  std::string tok;
  int cols = 0;
  while (head >> tok && tok.front() != '(') if (tok != "#") cols++;
  a.initialize(0.0, 80.0, 0.0, 90.0, 0.0);
  std::istringstream row(a.history_row(0.0));
  int vals = 0;
  while (row >> tok) vals++;
  CHECK_EQ(cols, 11);
  CHECK_EQ(cols, vals);
}

// ---- Circulation models with the mock ---------------------------------------

namespace
{
const char * kEmptyMesh = "<mesh celltype=\"tetrahedron\" dim=\"3\"></mesh>";
}

TEST_CASE("RegazzoniCirculation: run plan from the time base, and the -t ceiling reduces beats")
{
  write_text("r.xml", kEmptyMesh);
  {
    TimeBase tb;
    tb.step_ms = 2.0;
    tb.beats = 3;
    tb.period_ms = 800.0;
    RegazzoniCirculation r;
    r.configure("r.xml", tb);
    CHECK_EQ(r.plan().beats, 3);
    CHECK_EQ(r.plan().nsteps, 3 * 400);
    CHECK_NEAR(r.plan().run_ms, 2400.0, 0.0);
    CHECK_NEAR(r.plan().period_ms, 800.0, 0.0);
  }
  {
    TimeBase tb;
    tb.ceiling_ms = 1700.0;
    tb.beats = 5;
    tb.period_ms = 800.0;
    RegazzoniCirculation r;
    r.configure("r.xml", tb);
    CHECK_EQ(r.plan().beats, 2);
  }
  {
    TimeBase tb;
    tb.ceiling_ms = 100.0;
    tb.period_ms = 800.0;
    RegazzoniCirculation r;
    r.configure("r.xml", tb);
    CHECK_MSG(r.plan().beats == 1, "never below one beat");
  }
  {
    // A coupling step larger than the period is refused
    TimeBase tb;
    tb.step_ms = 5000.0;
    tb.period_ms = 800.0;
    RegazzoniCirculation r;
    CHECK_THROWS(r.configure("r.xml", tb));
  }
}

TEST_CASE("RegazzoniCirculation: the pressure inversion reproduces the 0D volumes")
{
  write_text("r.xml", kEmptyMesh);
  TimeBase tb;
  tb.step_ms = 2.0;
  tb.beats = 2;
  RegazzoniCirculation r;
  r.configure("r.xml", tb);

  LumpedCavity mech;
  r.begin(mech);
  run_steps(r, 0, r.plan().nsteps);
  r.end();

  CHECK(mech.solves >= r.plan().nsteps);

  // The 0D state is what the 3D side must end up reproducing
  const arma::mat h = read_table("r.xml_circulation.dat");
  REQUIRE(h.n_rows == arma::uword(r.plan().nsteps));
  CHECK(h.is_finite());
  CHECK_NEAR(h(h.n_rows - 1, 0), r.plan().run_ms, 1e-9);

  // A beating ventricle: the volume moves, and comes back within the beat
  const arma::vec vlv = h.col(1);
  CHECK_MSG(vlv.max() - vlv.min() > 5.0, "LV volume range " << vlv.max() - vlv.min());
}

TEST_CASE("RegazzoniCirculation: conserves the stressed blood volume when coupled")
{
  write_text("r.xml", kEmptyMesh);
  TimeBase tb;
  tb.step_ms = 2.0;
  tb.beats = 1;
  RegazzoniCirculation r;
  r.configure("r.xml", tb);
  LumpedCavity mech;
  r.begin(mech);

  CirculationModel::State st0;
  r.save_state(st0);
  run_steps(r, 0, r.plan().nsteps);
  CirculationModel::State st1;
  r.save_state(st1);
  r.end();

  auto vtot = [](const CirculationModel::State & st) {
    Regazzoni2020 c;
    for (const auto & kv : st)
      if (kv.second.size() == Regazzoni2020::NSTATE)
      {
        for (int i = 0; i < Regazzoni2020::NSTATE; i++) c[i] = kv.second[i];
        return c.total_volume();
      }
    return std::nan("");
  };
  const double v0 = vtot(st0), v1 = vtot(st1);
  REQUIRE(std::isfinite(v0));
  CHECK_MSG(std::fabs(v1 - v0) < 1e-6 * v0, v0 << " -> " << v1);
}

TEST_CASE("PiecewiseCirculation: the plan is beats*period, and the constraints are met")
{
  write_text("p.xml", kEmptyMesh);
  TimeBase tb;
  tb.step_ms = 1.0;
  tb.beats = 1;
  tb.period_ms = 600.0;
  PiecewiseCirculation p;
  p.configure("p.xml", tb);

  const PiecewisePVLoop defaults;
  const RunPlan & plan = p.plan();
  CHECK_NEAR(plan.run_ms, 600.0, 0.0);
  CHECK_NEAR(plan.stim_offset_ms, defaults.lv().t_end_diastole, 0.0);
  CHECK_EQ(plan.nsteps, 600);

  LumpedCavity mech;
  mech.period = 600.0;
  mech.t_act = plan.stim_offset_ms;
  p.begin(mech);
  run_steps(p, 0, plan.nsteps);
  p.end();

  // Pressure steps take one solve; the others need the Newton on the pressures
  CHECK(mech.solves > plan.nsteps);

  const arma::mat h = read_table("p.xml_piecewise.dat");
  REQUIRE(h.n_rows == arma::uword(plan.nsteps));
  CHECK(h.is_finite());

  // Columns: t, then per cavity: phase, P, V, Q, P_c
  const PiecewisePVLoop::Windkessel wk[2] = { defaults.lv().windkessel, defaults.rv().windkessel };
  int seen[2][5] = {};
  for (int c = 0; c < 2; c++)
  {
    const int k = 1 + 5 * c;
    double v_locked = 0.0;
    for (arma::uword i = 0; i < h.n_rows; i++)
    {
      const int ph = static_cast<int>(h(i, k));
      REQUIRE(ph >= 0);
      REQUIRE(ph <= 4);
      seen[c][ph]++;
      const double P = h(i, k + 1), V = h(i, k + 2), Q = h(i, k + 3), Pc = h(i, k + 4);
      if (ph == PiecewisePVLoop::PRELOAD || ph == PiecewisePVLoop::EJECTION) v_locked = V;
      if (ph == PiecewisePVLoop::IVC || ph == PiecewisePVLoop::IVR)
        CHECK_MSG(std::fabs(V - v_locked) < 1e-3,
                  "cavity " << c << ", t = " << h(i, 0) << ": V moved " << V - v_locked << " mL");
      if (ph == PiecewisePVLoop::FILLING && i > 0)
        CHECK_NEAR(V - h(i - 1, k + 2), (c == 0 ? defaults.lv() : defaults.rv()).fill_rate, 1e-3);
      // Windkessel: the cavity pressure is P_c + R_c Q at the end of the step
      if (ph == PiecewisePVLoop::EJECTION)
        CHECK_NEAR(P, Pc + wk[c].Rc * Q, 1e-3 * (wk[c].Rc + 1.0 / wk[c].C));
    }
  }
  for (int c = 0; c < 2; c++)
    for (int ph = 0; ph < 5; ph++)
      CHECK_MSG(seen[c][ph] > 0, "cavity " << c << " never solved under phase " << ph);
}

namespace
{

template <typename Model>
void check_restart_is_exact(const std::string & mesh, const std::string & history,
                            double period, double t_act)
{
  write_text(mesh, kEmptyMesh);
  TimeBase tb;
  tb.step_ms = 1.0;
  tb.beats = 1;
  tb.period_ms = period;

  // Reference, straight through
  Model full;
  full.configure(mesh, tb);
  const int n = full.plan().nsteps;
  const int cut = n / 2;
  LumpedCavity mech_full;
  mech_full.period = period;
  mech_full.t_act = t_act > 0.0 ? t_act : full.plan().stim_offset_ms;
  full.begin(mech_full);

  run_steps(full, 0, cut);
  CirculationModel::State saved;
  full.save_state(saved);
  LumpedCavity mech_at_cut = mech_full;
  run_steps(full, cut, n);
  CirculationModel::State ref;
  full.save_state(ref);
  full.end();

  // Restarted from the saved state, in a new object
  Model again;
  again.configure(mesh, tb);
  LumpedCavity mech_again = mech_at_cut;
  again.begin(mech_again);
  again.restore_state(saved);
  mech_again.vlv = mech_at_cut.vlv;
  mech_again.vrv = mech_at_cut.vrv;
  run_steps(again, cut, n);
  CirculationModel::State res;
  again.save_state(res);
  again.end();

  REQUIRE(ref.size() == res.size());
  for (const auto & kv : ref)
  {
    REQUIRE(res.count(kv.first) == 1);
    const std::vector<double> & b = res.at(kv.first);
    REQUIRE(kv.second.size() == b.size());
    for (std::size_t i = 0; i < b.size(); i++)
      CHECK_MSG(rel_gap(kv.second[i], b[i]) < 1e-12 || std::fabs(kv.second[i] - b[i]) < 1e-12,
                kv.first << "[" << i << "]: " << kv.second[i] << " vs " << b[i]);
  }
}

} // namespace

TEST_CASE("RegazzoniCirculation: save and restore continue the run exactly")
{
  check_restart_is_exact<RegazzoniCirculation>("rs.xml", "rs.xml_circulation.dat", 400.0, 0.0);
}

TEST_CASE("PiecewiseCirculation: save and restore continue the run exactly")
{
  check_restart_is_exact<PiecewiseCirculation>("ps.xml", "ps.xml_piecewise.dat", 400.0, 0.0);
}

TEST_CASE("CirculationModel: restore_state with a missing entry throws")
{
  write_text("r.xml", kEmptyMesh);
  TimeBase tb;
  RegazzoniCirculation r;
  r.configure("r.xml", tb);
  LumpedCavity mech;
  r.begin(mech);
  CirculationModel::State empty;
  CHECK_THROWS(r.restore_state(empty));

  PiecewiseCirculation p;
  write_text("p.xml", kEmptyMesh);
  p.configure("p.xml", tb);
  p.begin(mech);
  CHECK_THROWS(p.restore_state(empty));
}

// ---- LegacyCirculation ------------------------------------------------------

namespace
{

std::string pvloop_mesh(double total_ms, int steps, bool biv = false)
{
  std::ostringstream os;
  os << "<mesh celltype=\"tetrahedron\" dim=\"3\">"
     << "<pvloop size=\"" << steps << "\" total_time=\"" << total_ms << "\" passive_time=\"100\""
     << " C_art=\"0.014\" C_ven=\"0.3\" R_ao=\"3850\" V_ven_zero=\"3300\" P_o=\"500\" B_LA=\"0.049\""
     << " V_art_zero=\"580\" p_art=\"10800\" p_ven=\"1600\" R_mv=\"1750\" R_ven=\"1400\" R_per=\"140000\""
     << " E_es_LA=\"60\" A_LA=\"58.67\" T_max=\"200\" tau=\"25\" biv=\"" << (biv ? 1 : 0) << "\">";
  for (int i = 0; i < steps; i++)
  {
    const double t = total_ms * i / (steps - 1);
    os << "<step time=\"" << t << "\" ";
    if (biv) os << "pressure_lv=\"" << 100.0 * i << "\" pressure_rv=\"" << 30.0 * i << "\"";
    else     os << "pressure=\"" << 100.0 * i << "\"";
    os << "/>";
  }
  os << "</pvloop></mesh>";
  return os.str();
}

} // namespace

TEST_CASE("LegacyCirculation: the coupling step is the sampling of the curve, and -t cuts it")
{
  write_text("l.xml", pvloop_mesh(100.0, 11));
  {
    LegacyCirculation l;
    TimeBase tb;
    tb.step_ms = 0.5;   // ignored by the legacy model
    l.configure("l.xml", tb);
    CHECK_EQ(l.plan().nsteps, 10);
    CHECK_NEAR(l.plan().step_ms, 10.0, 1e-12);
    CHECK_NEAR(l.plan().run_ms, 100.0, 0.0);
    CHECK_NEAR(l.plan().ep_dt_hint_ms, 10.0, 1e-12);
    CHECK(l.plan().store_stress);
  }
  {
    LegacyCirculation l;
    TimeBase tb;
    tb.ceiling_ms = 50.0;
    l.configure("l.xml", tb);
    CHECK_EQ(l.plan().nsteps, 5);
  }
}

TEST_CASE("LegacyCirculation: malformed curves are rejected")
{
  write_text("e.xml", "<mesh><pvloop total_time=\"100\"></pvloop></mesh>");
  write_text("s.xml", pvloop_mesh(1.0, 5));   // total_time in seconds
  LegacyCirculation a, b;
  TimeBase tb;
  CHECK_THROWS(a.configure("e.xml", tb));
  CHECK_THROWS(b.configure("s.xml", tb));
}

namespace
{
//! The legacy Newton matches volume_total() to the LV of the 0D model: only one cavity here
struct OneCavity : public LumpedCavity
{
  std::vector<double> plv, prv;
  void apply_pressures(double t, double a, double b) override
  {
    plv.push_back(a);
    prv.push_back(b);
    LumpedCavity::apply_pressures(t, a, b);
  }
  double volume_total() override { return vlv; }
};
}

TEST_CASE("LegacyCirculation: Newton on the volume, with the Windkessel volumes conserved")
{
  write_text("l.xml", pvloop_mesh(100.0, 11));
  LegacyCirculation l;
  TimeBase tb;
  l.configure("l.xml", tb);

  OneCavity mech;
  l.begin(mech);
  run_steps(l, 0, l.plan().nsteps);
  l.end();

  REQUIRE(mech.plv.size() > std::size_t(l.plan().nsteps));
  for (double p : mech.plv) CHECK(std::isfinite(p));
  CHECK_NEAR(mech.prv.back(), 0.0, 0.0);   // no RV curve without biv

  // t, p_lv, V_lv, 0, p_art, p_ven, p_la, V_art, V_ven, V_la, ...
  const arma::mat h = read_table("l.xml_pvloop.txt");
  REQUIRE(h.n_rows == 10);
  CHECK(h.is_finite());
  arma::vec total = h.col(2) + h.col(7) + h.col(8) + h.col(9);
  CHECK_MSG(total.max() - total.min() < 1.0, "blood volume drift " << total.max() - total.min() << " mL");
}

TEST_CASE("LegacyCirculation: with biv, the RV pressure comes from the table")
{
  write_text("b.xml", pvloop_mesh(100.0, 11, true));
  LegacyCirculation l;
  TimeBase tb;
  l.configure("b.xml", tb);
  OneCavity mech;
  l.begin(mech);
  run_steps(l, 0, 5);
  CHECK_NEAR(mech.prv.back(), 30.0 * 5, 1e-9);
  l.end();
}

TEST_CASE("LegacyCirculation: save and restore continue the run exactly")
{
  write_text("l.xml", pvloop_mesh(100.0, 21));
  TimeBase tb;

  LegacyCirculation full;
  full.configure("l.xml", tb);
  OneCavity m1;
  full.begin(m1);
  run_steps(full, 0, 10);
  CirculationModel::State saved;
  full.save_state(saved);
  OneCavity m_cut = m1;
  run_steps(full, 10, 20);
  CirculationModel::State ref;
  full.save_state(ref);
  full.end();

  LegacyCirculation again;
  again.configure("l.xml", tb);
  OneCavity m2 = m_cut;
  again.begin(m2);
  again.restore_state(saved);
  m2.vlv = m_cut.vlv;
  m2.vrv = m_cut.vrv;
  run_steps(again, 10, 20);
  CirculationModel::State res;
  again.save_state(res);
  again.end();

  for (const auto & kv : ref)
  {
    REQUIRE(res.count(kv.first) == 1);
    REQUIRE(res.at(kv.first).size() == kv.second.size());
    for (std::size_t i = 0; i < kv.second.size(); i++)
      CHECK_MSG(std::fabs(kv.second[i] - res.at(kv.first)[i]) <= 1e-10 * (1.0 + std::fabs(kv.second[i])),
                kv.first << "[" << i << "]: " << kv.second[i] << " vs " << res.at(kv.first)[i]);
  }
}

// ---- cell_coupling ----------------------------------------------------------

TEST_CASE("CellProfile: ToRORd-Land and Kerckhoffs mappings, and unknown models")
{
  const CellProfile t = CellProfile::for_model("ToRORdLand", 1.0);
  CHECK_EQ(t.ta_var, 49);
  CHECK_NEAR(t.ta_unit, 1000.0, 0.0);
  CHECK_EQ(t.vm_var, 0);
  CHECK_EQ(t.land_xs, 43);
  CHECK_EQ(t.land_xw, 44);
  CHECK_EQ(t.land_zetas, 47);
  CHECK_EQ(t.land_zetaw, 48);
  CHECK(t.uses_stretch);
  CHECK(t.stim_amplitude < 0.0);

  const CellProfile k = CellProfile::for_model("Kerkoff2003", 7.5);
  CHECK_EQ(k.ta_var, -1);
  CHECK_NEAR(k.ta_unit, 7.5, 0.0);
  CHECK_EQ(k.vm_var, -1);
  CHECK_NEAR(k.stim_amplitude, 0.0, 0.0);
  CHECK(!k.uses_stretch);

  CHECK_THROWS(CellProfile::for_model("TT2", 1.0));
}

TEST_CASE("TaCurve: interpolation, the ends of the table, and the period")
{
  write_text("ta.txt", "# t Ta\n0 0\n100 1000\n200 500   # comment\n\n300 0\n");
  TaCurve c;
  CHECK(!c.active());
  c.load("ta.txt", -1.0);
  CHECK(c.active());
  CHECK_NEAR(c(50.0), 500.0, 1e-12);
  CHECK_NEAR(c(150.0), 750.0, 1e-12);
  CHECK_NEAR(c(100.0), 1000.0, 1e-12);
  CHECK_NEAR(c(-5.0), 0.0, 0.0);
  CHECK_NEAR(c(1000.0), 0.0, 0.0);

  TaCurve p;
  p.load("ta.txt", 400.0);
  CHECK_NEAR(p(450.0), 500.0, 1e-12);
  CHECK_NEAR(p(350.0), 0.0, 0.0);
  CHECK_NEAR(p(-350.0), 500.0, 1e-12);
}

TEST_CASE("TaCurve: bad files are rejected")
{
  TaCurve c;
  CHECK_THROWS(c.load("no_such_file.txt", -1.0));
  write_text("short.txt", "0 0\n");
  CHECK_THROWS(c.load("short.txt", -1.0));
  write_text("back.txt", "0 0\n10 1\n5 2\n");
  CHECK_THROWS(c.load("back.txt", -1.0));
}

TEST_CASE("LandCoupling: nodal active stiffness from the crossbridge states")
{
  LandCoupling land;
  CHECK(!land.has_stiffness_states());
  const CellProfile t = CellProfile::for_model("ToRORdLand", 1.0);
  land.set_state_indices(t.land_xs, t.land_xw, t.land_zetas, t.land_zetaw);
  CHECK(land.has_stiffness_states());

  Cells cells(3, "ToRORdLand", "ExplicitEuler", 0.02, 10.0, 1.0);
  cells.set_solver_time_unit_ms(1.0);
  cells.init();
  arma::vec xs = {0.1, 0.2, 0.0}, xw = {0.05, 0.0, 0.0};
  cells.set_var(t.land_xs, xs);
  cells.set_var(t.land_xw, xw);

  arma::vec lambda = {1.0, 1.1, 1.0};
  arma::vec ka;
  land.active_stiffness(cells, lambda, ka);
  REQUIRE(ka.n_elem == 3);
  CHECK_NEAR(ka(0), land_active_stiffness(1.0, 0.1, 0.05), 1e-12);
  CHECK_NEAR(ka(1), land_active_stiffness(1.1, 0.2, 0.0), 1e-12);
  CHECK_NEAR(ka(2), 0.0, 0.0);
}

// ---- UVC transfer -----------------------------------------------------------

namespace
{

// Regular grid in (ab, tm, rt) with xyz = uvc, all in the LV (tv = 0)
UVCDataTransfer uvc_grid(int n, const std::string & label)
{
  UVCDataTransfer m;
  m.label = label;
  m.nen = 4;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
      {
        const double a = i / (n - 1.0), t = j / (n - 1.0), r = k / (n - 1.0);
        m.xyz.push_back(a);
        m.xyz.push_back(t);
        m.xyz.push_back(r);
        m.ab.push_back(a);
        m.tm.push_back(t);
        m.rt.push_back(r);
        m.tv.push_back(0.0);
      }
  m.n_points = n * n * n;
  return m;
}

} // namespace

TEST_CASE("UVC transfer: nodes that coincide copy the source value")
{
  UVCDataTransfer src = uvc_grid(5, "source"), tgt = uvc_grid(5, "target");
  UVCParameters p;
  p.quiet = true;
  p.k = 1;
  KdtreeUVC kd;
  REQUIRE(kdtree_uvc_build(src, tgt, p, kd));
  CHECK_EQ(kd.n_source, src.n_points);
  CHECK_EQ(kd.n_target, tgt.n_points);

  std::vector<double> f(src.n_points), out;
  for (int i = 0; i < src.n_points; i++) f[i] = src.xyz[3*i] + 2.0 * src.xyz[3*i+1];
  REQUIRE(kdtree_uvc_transfer(kd, f, false, out));
  REQUIRE(int(out.size()) == tgt.n_points);
  for (int i = 0; i < tgt.n_points; i++)
    CHECK_NEAR(out[i], f[i], 1e-12);
}

TEST_CASE("UVC transfer: smooth fields stay in range, weights are a partition of unity")
{
  UVCDataTransfer src = uvc_grid(6, "source"), tgt = uvc_grid(4, "target");
  UVCParameters p;
  p.quiet = true;
  KdtreeUVC kd;
  REQUIRE(kdtree_uvc_build(src, tgt, p, kd));

  for (int i = 0; i < kd.n_target; i++)
  {
    double sum = 0.0;
    for (int j = 0; j < kd.n_neighbors[i]; j++) sum += kd.weights[kd.k_max * i + j];
    CHECK_NEAR(sum, 1.0, 1e-10);
  }

  std::vector<double> f(src.n_points), out;
  for (int i = 0; i < src.n_points; i++) f[i] = 3.0 + src.xyz[3*i];
  REQUIRE(kdtree_uvc_transfer(kd, f, false, out));
  for (double v : out)
  {
    CHECK(v >= 3.0 - 1e-12);
    CHECK(v <= 4.0 + 1e-12);
  }
}

TEST_CASE("UVC transfer: categorical fields keep their labels")
{
  UVCDataTransfer src = uvc_grid(5, "source"), tgt = uvc_grid(4, "target");
  UVCParameters p;
  p.quiet = true;
  KdtreeUVC kd;
  REQUIRE(kdtree_uvc_build(src, tgt, p, kd));
  std::vector<double> f(src.n_points), out;
  for (int i = 0; i < src.n_points; i++) f[i] = src.xyz[3*i] > 0.5 ? 2.0 : 7.0;
  REQUIRE(kdtree_uvc_transfer(kd, f, true, out));
  for (int i = 0; i < kd.n_target; i++)
    CHECK(out[i] == 2.0 || out[i] == 7.0);
  CHECK_NEAR(out[0], 7.0, 0.0);
  CHECK_NEAR(out[kd.n_target - 1], 2.0, 0.0);
}

TEST_CASE("UVC transfer: non-finite source values are skipped, wrong sizes are refused")
{
  UVCDataTransfer src = uvc_grid(4, "source"), tgt = uvc_grid(3, "target");
  UVCParameters p;
  p.quiet = true;
  KdtreeUVC kd;
  REQUIRE(kdtree_uvc_build(src, tgt, p, kd));
  // Node (1,1,1) of the source grid is at 1/3 and coincides with no target node. A
  // coincident neighbour would take the whole Gaussian weight, leaving nothing to renormalize.
  std::vector<double> f(src.n_points, 1.0), out;
  f[1 + 4 + 16] = std::nan("");
  REQUIRE(kdtree_uvc_transfer(kd, f, false, out));
  for (double v : out) CHECK_NEAR(v, 1.0, 1e-12);

  std::vector<double> too_short(3, 1.0);
  CHECK(!kdtree_uvc_transfer(kd, too_short, false, out));

  UVCDataTransfer none;
  KdtreeUVC kd2;
  CHECK(!kdtree_uvc_build(none, tgt, p, kd2));
}

TEST_CASE("distinct_values_uvc: sorted unique values")
{
  const std::vector<double> d = distinct_values_uvc({3.0, 1.0, 3.0, 2.0, 1.0});
  REQUIRE(d.size() == 3);
  CHECK_NEAR(d[0], 1.0, 0.0);
  CHECK_NEAR(d[2], 3.0, 0.0);
}
