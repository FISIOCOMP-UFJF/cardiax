#include "regazzoni_circulation.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "util/command_line_args.h"

using std::cout;
using std::endl;

void RegazzoniCirculation::configure(const std::string & mshfile, const TimeBase & tb)
{
  base_ = mshfile;

  double hr  = CommandLineArgs::read("-hr",  -1.0,
      "frequencia cardiaca em Hz (default 1.0 = 60 bpm)", "circulacao");
  double bpm = CommandLineArgs::read("-bpm", -1.0,
      "frequencia cardiaca em bpm (alternativa a -hr)", "circulacao");
  if (bpm > 0.0) hr = bpm / 60.0;
  if (hr  > 0.0) circ_.set_heart_rate(hr);

  const double tact = CommandLineArgs::read("-tact", -1.0,
      "instante da ativacao eletrica ventricular, em s", "circulacao");
  if (tact >= 0.0) circ_.set_activation_time(tact);

  const double emd = CommandLineArgs::read("-emd", -1.0,
      "atraso eletromecanico estimulo -> inicio de Ta, em s", "circulacao");
  if (emd >= 0.0) circ_.set_emd(emd);

  cout << "\n" << circ_.describe_timing() << endl;

  dt_s_ = CommandLineArgs::read("-dtc", 1.0e-3,
      "passo de tempo da circulacao, em s", "circulacao");
  steps_per_beat_ = static_cast<int>(std::round(circ_.period() / dt_s_));
  if (dt_s_ <= 0.0 || steps_per_beat_ < 1)
    throw std::runtime_error("-dtc invalido: precisa ser > 0 e menor que o periodo cardiaco");

  const double period_ms = circ_.period() * 1000.0;

  plan_.beats             = limit_beats(tb, tb.beats, period_ms);
  plan_.period_ms         = period_ms;
  plan_.run_ms            = plan_.beats * period_ms;
  plan_.step_ms           = dt_s_ * 1000.0;
  plan_.nsteps            = plan_.beats * steps_per_beat_;
  plan_.stim_offset_ms    = 0.0;
  plan_.span_margin_steps = 2;
  plan_.default_out_every = 8;
}

void RegazzoniCirculation::do_begin()
{
  mmHg_to_p_ = mech_->pressure_unit_per_mmHg();

  // Inicializa o 0D a partir dos volumes 3D reais.
  // Sobrescrever V_LV/V_RV com os da malha MUDA o volume sanguineo estressado:
  // o estado inicial publicado tem V_LV = 118.5 e V_RV = 166.2 mL, e uma malha
  // tipica esta dezenas de mL abaixo disso. Sem compensacao o modelo inteiro
  // roda subenchido, e o efeito e invisivel no log.
  const double vtot0 = CommandLineArgs::read("-circ_vtot",
      circ_.total_volume(),
      "volume sanguineo estressado em mL (default: o do estado inicial)",
      "circulacao");

  const double v_lv0 = mech_->volume_LV();
  const double v_rv0 = mech_->volume_RV();
  circ_.set_initial_volumes(v_lv0, v_rv0);
  const double dv_realloc = circ_.rebalance_total_volume(vtot0);

  cout << " Initial cavity volumes from the 3D mesh: LV = " << v_lv0
       << " mL, RV = " << v_rv0 << " mL" << endl;
  cout << " Stressed blood volume held at " << vtot0 << " mL ("
       << dv_realloc << " mL reallocated to the venous reservoirs)" << endl;
  cout << " Pressure conversion: 1 mmHg = " << mmHg_to_p_
       << " model pressure units" << endl;

  // O 0D chama de volta o 3D atraves deste acoplamento
  circ_.set_BiV_coupling(
    [this](double Vlv, double Vrv, double t, double & plv, double & prv)
    {
      this->pressures_for_volumes(Vlv, Vrv, t, plv, prv);
    });

  open_history(base_ + "_circulation.dat", Regazzoni2020::history_header());

  cout << " Heart period = " << circ_.period() << " s, dt = " << dt_s_
       << " s -> " << steps_per_beat_ << " steps/beat, " << plan_.beats
       << " beat(s)" << endl;
  cout << "\n Legend: [CIRC] 0D circulation | [PRESS] pressure inversion |"
       << " [STATE] volumes/pressures | [FLOW] valve flows" << endl;
}

void RegazzoniCirculation::announce(const StepContext & s)
{
  const int beat = s.index / steps_per_beat_ + 1;

  if (s.index % steps_per_beat_ == 0)
    cout << "\n--- beat " << beat << "/" << plan_.beats << " ---" << endl;

  if (s.log)
    print_status(s, "beat " + std::to_string(beat) + "/" + std::to_string(plan_.beats));
}

void RegazzoniCirculation::advance(const StepContext & s)
{
  // circ_.step() chama pressures_for_volumes() internamente, que roda as
  // solucoes 3D necessarias para casar os volumes que o 0D pede.
  if (s.log)
    cout << "  [CIRC ] advancing 0D model, dt = " << dt_s_ << " s" << endl;

  circ_.step(t_s_, dt_s_);
  t_s_ += dt_s_;

  if (s.log)
  {
    cout << "  [PRESS] " << last_newton_its_ << " Newton it, |dV| = "
         << last_dV_ << " mL, " << last_solves_ << " 3D solve(s)" << endl;
    cout << "  [STATE] LV: V = " << circ_[Regazzoni2020::V_LV]
         << " mL, p = " << circ_.pressure_LV() << " mmHg"
         << "  |  RV: V = " << circ_[Regazzoni2020::V_RV]
         << " mL, p = " << circ_.pressure_RV() << " mmHg" << endl;
    cout << "  [FLOW ] MV = " << circ_.flow_MV() << ", AV = " << circ_.flow_AV()
         << ", TV = " << circ_.flow_TV() << ", PV = " << circ_.flow_PV()
         << " mL/s  | Vtot = " << circ_.total_volume() << " mL" << endl;
  }

  if (history_.is_open())
  {
    history_ << circ_.history_row(t_s_) << "\n";
    history_.flush();
  }
}

void RegazzoniCirculation::end()
{
  if (history_.is_open()) history_.close();

  cout << "\n Total 3D solves: " << solves_
       << " (" << double(solves_) / std::max(1, plan_.nsteps)
       << " per circulation step)" << endl;
  cout << " Final blood volume: " << circ_.total_volume() << " mL" << endl;
}

// ----------------------------------------------------------------------
//  Inversao do mapa p -> V
// ----------------------------------------------------------------------

void RegazzoniCirculation::update_pv_jacobian(double t)
{
  const double dp = 0.05 * mmHg_to_p_;   // 0.05 mmHg perturbation

  // baseline
  mech_->apply_pressures(t * 1000.0, p_lv_cur_, p_rv_cur_);
  solves_++;
  double v0_lv = mech_->volume_LV();
  double v0_rv = mech_->volume_RV();

  // perturb LV pressure
  mech_->apply_pressures(t * 1000.0, p_lv_cur_ + dp, p_rv_cur_);
  solves_++;
  Jpv_(0,0) = (mech_->volume_LV() - v0_lv) / dp;
  Jpv_(1,0) = (mech_->volume_RV() - v0_rv) / dp;

  // perturb RV pressure
  mech_->apply_pressures(t * 1000.0, p_lv_cur_, p_rv_cur_ + dp);
  solves_++;
  Jpv_(0,1) = (mech_->volume_LV() - v0_lv) / dp;
  Jpv_(1,1) = (mech_->volume_RV() - v0_rv) / dp;

  Jpv_valid_ = true;

  if (mech_->verbose())
    cout << "    [PRESS] dV/dp = [" << Jpv_(0,0) << " " << Jpv_(0,1) << "; "
         << Jpv_(1,0) << " " << Jpv_(1,1) << "] mL per pressure unit" << endl;

  if (std::fabs(arma::det(Jpv_)) < 1e-12)
  {
    cout << "   [coupling] WARNING: dV/dp is near singular; the two cavities"
         << " respond almost identically to both pressures." << endl;
  }
}

void RegazzoniCirculation::pressures_for_volumes(double Vlv_target,
                                                    double Vrv_target,
                                                    double t,
                                                    double & plv_mmHg,
                                                    double & prv_mmHg)
{
  // INVARIANT: this routine must NOT advance the cell model. It performs
  // several 3D solves per circulation step, and the 0D solver may call it
  // more than once per step; advancing ephy here would make the active
  // tension run ahead of the circulation clock by a factor equal to the
  // number of solves. `ta` is set once per step in solve_circulation() and
  // is held fixed throughout the inversion.
  const double tol_vol = 0.05;   // mL
  const int    max_it  = 8;

  int solves0 = solves_;

  if (!Jpv_valid_)
  {
    if (mech_->verbose())
      cout << "    [PRESS] rebuilding dV/dp Jacobian (3 solves)" << endl;
    update_pv_jacobian(t);
  }

  arma::vec2 F;
  int it = 0;

  // evaluate at the current pressures
  mech_->apply_pressures(t * 1000.0, p_lv_cur_, p_rv_cur_);
  solves_++;
  F(0) = mech_->volume_LV() - Vlv_target;
  F(1) = mech_->volume_RV() - Vrv_target;

  double f0 = arma::norm(F, 2);

  while (arma::norm(F, 2) > tol_vol && it < max_it)
  {
    arma::vec2 dp;
    bool ok = arma::solve(dp, Jpv_, -F);
    if (!ok)
    {
      cout << "   [coupling] Jacobian solve failed; rebuilding" << endl;
      update_pv_jacobian(t);
      ok = arma::solve(dp, Jpv_, -F);
      if (!ok) break;
    }

    // damped update: cavity pressures should not jump wildly in one step
    double maxstep = 20.0 * mmHg_to_p_;
    double n = arma::norm(dp, 2);
    if (n > maxstep) dp *= maxstep / n;

    p_lv_cur_ += dp(0);
    p_rv_cur_ += dp(1);

    // Physiological bracket. If the 0D model asks for a volume the 3D model
    // cannot reach (too stiff, or a target outside its operating range),
    // Newton would otherwise run away to absurd pressures and take the whole
    // simulation with it.
    const double p_min = -20.0 * mmHg_to_p_;
    const double p_max = 400.0 * mmHg_to_p_;
    bool clipped = false;
    if (p_lv_cur_ < p_min) { p_lv_cur_ = p_min; clipped = true; }
    if (p_lv_cur_ > p_max) { p_lv_cur_ = p_max; clipped = true; }
    if (p_rv_cur_ < p_min) { p_rv_cur_ = p_min; clipped = true; }
    if (p_rv_cur_ > p_max) { p_rv_cur_ = p_max; clipped = true; }
    if (clipped)
    {
      cout << "   [coupling] WARNING t=" << t << ": pressure hit the"
           << " physiological bracket. The requested volumes (LV "
           << Vlv_target << " mL, RV " << Vrv_target
           << " mL) may be unreachable for this 3D model." << endl;
      Jpv_valid_ = false;
      break;
    }

    mech_->apply_pressures(t * 1000.0, p_lv_cur_, p_rv_cur_);
    solves_++;
    F(0) = mech_->volume_LV() - Vlv_target;
    F(1) = mech_->volume_RV() - Vrv_target;

    it++;
  }

  // If Newton stalled, the Jacobian is stale: force a rebuild next time.
  if (arma::norm(F, 2) > tol_vol)
  {
    Jpv_valid_ = false;
    cout << "   [coupling] t=" << t << " residual |dV| = "
         << arma::norm(F, 2) << " mL after " << it
         << " iterations (target " << tol_vol << "); Jacobian will be rebuilt"
         << endl;
  }
  else if (it > 3 || arma::norm(F,2) > 0.5 * f0)
  {
    Jpv_valid_ = false;   // slow convergence -> refresh
  }

  last_newton_its_ = it;
  last_dV_         = arma::norm(F, 2);
  last_solves_     = solves_ - solves0;

  plv_mmHg = p_lv_cur_ / mmHg_to_p_;
  prv_mmHg = p_rv_cur_ / mmHg_to_p_;
}