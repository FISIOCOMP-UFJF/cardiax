#include "regazzoni_circulation.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "util/command_line_args.h"

using std::cout;
using std::endl;

// Limites do Newton, em Pa (equivalem a 0.05, 20 e 400 mmHg)
namespace
{
constexpr double kJacobianPerturbation = 6.666119371;
constexpr double kMaxNewtonStep        = 2666.447748;
constexpr double kPressureMin          = -2666.447748;
constexpr double kPressureMax          = 53328.95497;
}

void RegazzoniCirculation::configure(const std::string & mshfile, const TimeBase & tb)
{
  base_ = mshfile;

  if (tb.period_ms > 0.0) circ_.set_period(tb.period_ms);

  const double tact = CommandLineArgs::read("-tact", -1.0,
      "instante da ativacao ventricular no ciclo, ms (default 0)", "circulacao");
  if (tact >= 0.0) circ_.set_activation_time(tact);

  const double emd = CommandLineArgs::read("-emd", -1.0,
      "atraso eletromecanico estimulo -> inicio de Ta, ms (default 40)", "circulacao");
  if (emd >= 0.0) circ_.set_emd(emd);

  cout << "\n" << circ_.describe_timing() << endl;

  dt_ms_ = tb.step_ms;
  steps_per_beat_ = static_cast<int>(std::round(circ_.period() / dt_ms_));
  if (dt_ms_ <= 0.0 || steps_per_beat_ < 1)
    throw std::runtime_error("-dt_mech invalido: precisa ser > 0 e menor que o periodo cardiaco");

  plan_.beats             = limit_beats(tb, tb.beats, circ_.period());
  plan_.period_ms         = circ_.period();
  plan_.run_ms            = plan_.beats * circ_.period();
  plan_.step_ms           = dt_ms_;
  plan_.nsteps            = plan_.beats * steps_per_beat_;
  plan_.stim_offset_ms    = 0.0;
  plan_.span_margin_steps = 2;
  plan_.default_out_every = 8;
}

void RegazzoniCirculation::do_begin()
{
  // Sobrescrever V_LV/V_RV com os da malha MUDA o volume sanguineo estressado
  // (o estado publicado tem 118.5/166.2 mL; uma malha tipica esta dezenas de mL
  // abaixo). Sem compensar, o modelo roda subenchido e o log nao mostra.
  const double vtot0 = CommandLineArgs::read("-circ_vtot", circ_.total_volume(),
      "volume sanguineo estressado, mL (default: o do estado inicial)", "circulacao");

  const double v_lv0 = mech_->volume_LV();
  const double v_rv0 = mech_->volume_RV();
  circ_.set_initial_volumes(v_lv0, v_rv0);
  const double dv = circ_.rebalance_total_volume(vtot0);

  cout << " Initial cavity volumes (3D mesh): LV = " << v_lv0 << " mL, RV = "
       << v_rv0 << " mL\n Stressed blood volume held at " << vtot0 << " mL ("
       << dv << " mL reallocated to the venous reservoirs)" << endl;

  // O 0D chama de volta o 3D atraves deste acoplamento
  circ_.set_BiV_coupling(
    [this](double Vlv, double Vrv, double t, double & plv, double & prv)
    { this->pressures_for_volumes(Vlv, Vrv, t, plv, prv); });

  open_history(base_ + "_circulation.dat", Regazzoni2020::history_header());

  cout << " Heart period = " << circ_.period() << " ms, dt = " << dt_ms_ << " ms -> "
       << steps_per_beat_ << " steps/beat, " << plan_.beats << " beat(s)\n"
       << "\n Legend: [CIRC] 0D | [PRESS] pressure inversion | [STATE] volumes and"
       << " pressures | [FLOW] valve flows" << endl;
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
  // circ_.step() chama pressures_for_volumes(), que faz as solucoes 3D.
  circ_.step(t_ms_, dt_ms_);
  t_ms_ += dt_ms_;

  if (s.log)
  {
    cout << "  [PRESS] " << last_newton_its_ << " Newton it, |dV| = " << last_dV_
         << " mL, " << last_solves_ << " 3D solve(s)\n"
         << "  [STATE] LV: V = " << circ_[Regazzoni2020::V_LV] << " mL, p = "
         << circ_.pressure_LV() << " Pa  |  RV: V = " << circ_[Regazzoni2020::V_RV]
         << " mL, p = " << circ_.pressure_RV() << " Pa\n"
         << "  [FLOW ] MV = " << circ_.flow_MV() << ", AV = " << circ_.flow_AV()
         << ", TV = " << circ_.flow_TV() << ", PV = " << circ_.flow_PV()
         << " mL/ms | Vtot = " << circ_.total_volume() << " mL" << endl;
  }

  if (history_.is_open())
  {
    history_ << circ_.history_row(t_ms_) << "\n";
    history_.flush();
  }
}

void RegazzoniCirculation::end()
{
  if (history_.is_open()) history_.close();
  cout << "\n Total 3D solves: " << solves_ << " ("
       << double(solves_) / std::max(1, plan_.nsteps) << " per step)\n"
       << " Final blood volume: " << circ_.total_volume() << " mL" << endl;
}

// ----------------------------------------------------------------------
//  Inversao do mapa p -> V
// ----------------------------------------------------------------------

void RegazzoniCirculation::update_pv_jacobian(double t)
{
  const double dp = kJacobianPerturbation;

  mech_->apply_pressures(t, p_lv_cur_, p_rv_cur_);
  solves_++;
  const double v0_lv = mech_->volume_LV();
  const double v0_rv = mech_->volume_RV();

  mech_->apply_pressures(t, p_lv_cur_ + dp, p_rv_cur_);
  solves_++;
  Jpv_(0,0) = (mech_->volume_LV() - v0_lv) / dp;
  Jpv_(1,0) = (mech_->volume_RV() - v0_rv) / dp;

  mech_->apply_pressures(t, p_lv_cur_, p_rv_cur_ + dp);
  solves_++;
  Jpv_(0,1) = (mech_->volume_LV() - v0_lv) / dp;
  Jpv_(1,1) = (mech_->volume_RV() - v0_rv) / dp;

  Jpv_valid_ = true;

  if (mech_->verbose())
    cout << "    [PRESS] dV/dp = [" << Jpv_(0,0) << " " << Jpv_(0,1) << "; "
         << Jpv_(1,0) << " " << Jpv_(1,1) << "] mL/Pa" << endl;

  if (std::fabs(arma::det(Jpv_)) < 1e-12)
    cout << "   [coupling] WARNING: dV/dp is near singular; the two cavities"
         << " respond almost identically to both pressures." << endl;
}

void RegazzoniCirculation::pressures_for_volumes(double Vlv_target, double Vrv_target,
                                                 double t, double & plv, double & prv)
{
  const double tol_vol = 0.05;   // mL
  const int    max_it  = 8;
  const int    solves0 = solves_;

  if (!Jpv_valid_)
  {
    if (mech_->verbose()) cout << "    [PRESS] rebuilding dV/dp Jacobian (3 solves)" << endl;
    update_pv_jacobian(t);
  }

  arma::vec2 F;
  int it = 0;

  mech_->apply_pressures(t, p_lv_cur_, p_rv_cur_);
  solves_++;
  F(0) = mech_->volume_LV() - Vlv_target;
  F(1) = mech_->volume_RV() - Vrv_target;
  const double f0 = arma::norm(F, 2);

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

    // Passo amortecido: a pressao da cavidade nao pode saltar de um passo a outro
    const double n = arma::norm(dp, 2);
    if (n > kMaxNewtonStep) dp *= kMaxNewtonStep / n;

    p_lv_cur_ += dp(0);
    p_rv_cur_ += dp(1);

    // Faixa fisiologica: se o 0D pede um volume que o 3D nao alcanca, o Newton
    // fugiria para pressoes absurdas e levaria a simulacao junto.
    bool clipped = false;
    for (double * p : { &p_lv_cur_, &p_rv_cur_ })
    {
      if (*p < kPressureMin) { *p = kPressureMin; clipped = true; }
      if (*p > kPressureMax) { *p = kPressureMax; clipped = true; }
    }
    if (clipped)
    {
      cout << "   [coupling] WARNING t=" << t << " ms: pressure hit the physiological"
           << " bracket. The requested volumes (LV " << Vlv_target << " mL, RV "
           << Vrv_target << " mL) may be unreachable for this 3D model." << endl;
      Jpv_valid_ = false;
      break;
    }

    mech_->apply_pressures(t, p_lv_cur_, p_rv_cur_);
    solves_++;
    F(0) = mech_->volume_LV() - Vlv_target;
    F(1) = mech_->volume_RV() - Vrv_target;
    it++;
  }

  // Newton estagnou, ou convergiu devagar: o Jacobiano esta velho
  if (arma::norm(F, 2) > tol_vol)
  {
    Jpv_valid_ = false;
    cout << "   [coupling] t=" << t << " ms: residual |dV| = " << arma::norm(F, 2)
         << " mL after " << it << " iterations (target " << tol_vol
         << "); Jacobian will be rebuilt" << endl;
  }
  else if (it > 3 || arma::norm(F, 2) > 0.5 * f0)
    Jpv_valid_ = false;

  last_newton_its_ = it;
  last_dV_         = arma::norm(F, 2);
  last_solves_     = solves_ - solves0;

  plv = p_lv_cur_;
  prv = p_rv_cur_;
}