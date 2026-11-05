#include "piecewise_circulation.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "util/command_line_args.h"

using std::cout;
using std::endl;

namespace
{
typedef PiecewisePVLoop::ChamberState ChamberState;

constexpr double kJacobianPerturbation = 6.666119371;   // Pa (0.05 mmHg)
constexpr double kMaxNewtonStep        = 5000.0;        // Pa per iteration

// ChamberState <-> vector, field by field in declaration order
std::vector<double> pack(const ChamberState & s)
{
  return { double(s.phase), double(s.n_beats), s.last_phase_change, s.volume_n,
           s.pressure_n, s.dvol_n, s.wdk_pressure_n, s.outflow_n,
           s.has_ejected ? 1.0 : 0.0, s.end_dia_vol, s.end_sys_vol,
           double(s.solved_under) };
}

ChamberState unpack(const std::vector<double> & v)
{
  if (v.size() != 12)
    throw std::runtime_error("Checkpoint: estado de cavidade do -pw com tamanho errado");
  ChamberState s;
  s.phase             = static_cast<PiecewisePVLoop::Phase>(static_cast<int>(v[0]));
  s.n_beats           = static_cast<int>(v[1]);
  s.last_phase_change = v[2];
  s.volume_n          = v[3];
  s.pressure_n        = v[4];
  s.dvol_n            = v[5];
  s.wdk_pressure_n    = v[6];
  s.outflow_n         = v[7];
  s.has_ejected       = (v[8] != 0.0);
  s.end_dia_vol       = v[9];
  s.end_sys_vol       = v[10];
  s.solved_under      = static_cast<PiecewisePVLoop::Phase>(static_cast<int>(v[11]));
  return s;
}
} // namespace

void PiecewiseCirculation::configure(const std::string & mshfile, const TimeBase & tb)
{
  base_ = mshfile;

  pw_.configure(mshfile, tb.period_ms);
  cout << "\n" << pw_.describe() << endl;

  tol_V_  = CommandLineArgs::read("-pwtol", tol_V_,
      "tolerancia das restricoes de volume/Windkessel, mL", "pvloop por fases");
  max_it_ = CommandLineArgs::read("-pwmaxit", max_it_,
      "iteracoes maximas do Newton nas pressoes por passo", "pvloop por fases");

  // Run length is beats*period from t = 0
  plan_.beats             = limit_beats(tb, tb.beats, pw_.period());
  plan_.period_ms         = pw_.period();
  plan_.run_ms            = pw_.total_time(plan_.beats);
  plan_.step_ms           = tb.step_ms;
  plan_.nsteps            = static_cast<int>(std::round(plan_.run_ms / tb.step_ms));
  // Activation starts at the end of diastole
  plan_.stim_offset_ms    = pw_.activation_time();
  plan_.span_margin_steps = 2;
  plan_.default_out_every = 8;
}

void PiecewiseCirculation::do_begin()
{
  // Start from the unloaded mesh volume at zero pressure
  const double v_lv0 = mech_->volume_LV();
  const double v_rv0 = mech_->volume_RV();
  pw_.initialize(0.0, v_lv0, 0.0, v_rv0, 0.0);
  p_cur_[0] = p_cur_[1] = 0.0;
  Jpv_.zeros();

  cout << " Volumes iniciais da malha: VE = " << v_lv0 << " mL, VD = " << v_rv0 << " mL" << endl;
  if (v_rv0 > v_lv0)
    cout << " *** NOTA: o volume da cavidade DIREITA excede o da esquerda nesta malha."
         << " Confirme os marcadores (30 = VE, 20 = VD)." << endl;

  open_history(base_ + "_piecewise.dat", PiecewisePVLoop::history_header());

  cout << " Rodada: " << plan_.run_ms << " ms = " << plan_.nsteps << " passos de "
       << plan_.step_ms << " ms (" << plan_.beats << " batimento(s) de "
       << plan_.period_ms << " ms)\n"
       << " Restricoes impostas por Newton nas pressoes: |residuo| <= " << tol_V_
       << " mL, ate " << max_it_ << " iteracoes\n"
       << "\n Legenda: [PW] fases | [PRESS] Newton nas pressoes" << endl;
}

void PiecewiseCirculation::announce(const StepContext & s)
{
  const int beat = pw_.beat_of(s.t_ms);
  if (beat != beat_reported_ && beat >= 0)
  {
    cout << "\n--- batimento " << beat + 1 << "/" << plan_.beats << " ---" << endl;
    beat_reported_ = beat;
  }
  if (s.log)
    print_status(s, std::string("VE ") + PiecewisePVLoop::phase_name(pw_.state_LV().phase)
                  + " | VD " + PiecewisePVLoop::phase_name(pw_.state_RV().phase));
}

// ----------------------------------------------------------------------
//  Constraints -> pressures
// ----------------------------------------------------------------------

void PiecewiseCirculation::mech_solve(double t, const double p[2], double V[2])
{
  mech_->apply_pressures(t, p[0], p[1]);
  solves_++;
  V[0] = mech_->volume_LV();
  V[1] = mech_->volume_RV();
}

void PiecewiseCirculation::rebuild_jacobian(double t, const double p[2], const double V[2],
                                            const bool free[2])
{
  for (int j = 0; j < 2; j++)
  {
    if (!free[j]) continue;
    double pp[2] = { p[0], p[1] }, Vp[2];
    pp[j] += kJacobianPerturbation;
    mech_solve(t, pp, Vp);
    Jpv_(0, j) = (Vp[0] - V[0]) / kJacobianPerturbation;
    Jpv_(1, j) = (Vp[1] - V[1]) / kJacobianPerturbation;
    col_valid_[j] = true;
  }
  if (mech_->verbose())
    cout << "    [PRESS] dV/dp = [" << Jpv_(0,0) << " " << Jpv_(0,1) << "; "
         << Jpv_(1,0) << " " << Jpv_(1,1) << "] mL/Pa" << endl;
}

void PiecewiseCirculation::solve_constrained(double t, const Constraint c[2],
                                             double p[2], double V[2])
{
  // Initial guess: last converged pressure; prescribed cavities take P
  bool free[2];
  int n_free = 0;
  for (int i = 0; i < 2; i++)
  {
    free[i] = (c[i].kind != Constraint::PRESSURE);
    p[i] = free[i] ? p_cur_[i] : c[i].P;
    n_free += free[i];
  }

  mech_solve(t, p, V);
  last_its_ = 0;
  last_res_ = 0.0;
  if (n_free == 0) return;

  // Residual in mL. AFFINE: p = A + B V  <=>  V = (p - A)/B
  auto residual = [&](const double pp[2], const double VV[2], double F[2])
  {
    double r = 0.0;
    for (int i = 0; i < 2; i++)
    {
      F[i] = 0.0;
      if (c[i].kind == Constraint::VOLUME) F[i] = VV[i] - c[i].V;
      if (c[i].kind == Constraint::AFFINE) F[i] = VV[i] - (pp[i] - c[i].A) / c[i].B;
      r = std::max(r, std::fabs(F[i]));
    }
    return r;
  };
  // dF_i/dp_i of the constraint term (0 for VOLUME)
  double g[2];
  for (int i = 0; i < 2; i++)
    g[i] = (c[i].kind == Constraint::AFFINE) ? -1.0 / c[i].B : 0.0;

  double F[2];
  double res = residual(p, V, F);
  int it = 0;
  bool rebuilt_here = false;

  while (res > tol_V_ && it < max_it_)
  {
    if ((free[0] && !col_valid_[0]) || (free[1] && !col_valid_[1]))
    {
      rebuild_jacobian(t, p, V, free);
      rebuilt_here = true;
    }

    // Newton system restricted to the free cavities
    arma::vec2 dp(arma::fill::zeros);
    if (n_free == 2)
    {
      arma::mat22 J = Jpv_;
      J(0, 0) += g[0];
      J(1, 1) += g[1];
      arma::vec2 rhs = { -F[0], -F[1] };
      if (!arma::solve(dp, J, rhs, arma::solve_opts::no_approx))
        dp.zeros();
    }
    else
    {
      const int i = free[0] ? 0 : 1;
      const double Jii = Jpv_(i, i) + g[i];
      if (Jii != 0.0) dp(i) = -F[i] / Jii;
    }

    if (!dp.is_finite() || arma::norm(dp, 2) == 0.0)
    {
      if (rebuilt_here) break;
      col_valid_[0] = col_valid_[1] = false;
      continue;
    }
    const double n = arma::norm(dp, 2);
    if (n > kMaxNewtonStep) dp *= kMaxNewtonStep / n;

    double p_new[2] = { p[0] + dp(0), p[1] + dp(1) }, V_new[2];
    mech_solve(t, p_new, V_new);
    it++;

    // Broyden update of dV/dp (only free columns change: dp is 0 elsewhere)
    arma::vec2 dV = { V_new[0] - V[0], V_new[1] - V[1] };
    Jpv_ += ((dV - Jpv_ * dp) * dp.t()) / arma::dot(dp, dp);

    p[0] = p_new[0]; p[1] = p_new[1];
    V[0] = V_new[0]; V[1] = V_new[1];
    const double res_old = res;
    res = residual(p, V, F);

    // Poor convergence: rebuild the Jacobian on the next iteration
    if (res > 0.5 * res_old && !rebuilt_here) col_valid_[0] = col_valid_[1] = false;
  }

  if (res > tol_V_)
  {
    failures_++;
    col_valid_[0] = col_valid_[1] = false;
    cout << "   [PRESS] *** t=" << t << " ms: restricao nao satisfeita, |residuo| = "
         << res << " mL apos " << it << " iteracoes (tol " << tol_V_ << ")" << endl;
  }
  last_its_ = it;
  last_res_ = res;
}

void PiecewiseCirculation::advance(const StepContext & s)
{
  // 1. constraints from the state at the start of the step
  const Constraint c[2] = { pw_.constraint_LV(s.t_ms, s.dt_ms),
                            pw_.constraint_RV(s.t_ms, s.dt_ms) };

  // 2. impose them with the Newton on the pressures
  double p[2], V[2];
  const int solves0 = solves_;
  solve_constrained(s.t_ms, c, p, V);
  p_cur_[0] = p[0];
  p_cur_[1] = p[1];

  // 3. Windkessel and phase transitions on the converged state
  pw_.commit(s.t_ms, s.dt_ms, V[0], p[0], V[1], p[1]);

  if (s.log)
  {
    const ChamberState & l = pw_.state_LV();
    const ChamberState & r = pw_.state_RV();
    cout << "  [PRESS] " << last_its_ << " it, |res| = " << last_res_ << " mL, "
         << solves_ - solves0 << " solucao(oes) 3D\n"
         << "  [PW   ] VE: P = " << l.pressure_n << " Pa, V = " << l.volume_n << " mL, P_c = "
         << l.wdk_pressure_n << " Pa, Q = " << l.outflow_n << " mL/ms | VD: P = "
         << r.pressure_n << " Pa, V = " << r.volume_n << " mL, P_c = " << r.wdk_pressure_n
         << " Pa" << endl;
  }

  if (history_.is_open())
  {
    history_ << pw_.history_row(s.t_ms) << "\n";
    history_.flush();
  }
}

void PiecewiseCirculation::end()
{
  if (history_.is_open()) history_.close();

  cout << "\n Total de solucoes 3D: " << solves_ << " ("
       << (plan_.nsteps > 0 ? double(solves_) / plan_.nsteps : 0.0) << " por passo)" << endl;
  if (failures_ > 0)
    cout << " *** ATENCAO: " << failures_ << " passo(s) terminaram sem satisfazer a"
         << " restricao (-pwtol " << tol_V_ << " mL). Aumente -pwmaxit." << endl;

  const ChamberState & l = pw_.state_LV();
  const ChamberState & r = pw_.state_RV();
  cout << " VE: V_ED = " << l.end_dia_vol << " mL, V_ES = " << l.end_sys_vol << " mL, VS = "
       << pw_.stroke_volume_LV() << " mL, FE = " << pw_.ejection_fraction_LV() << " %\n"
       << " VD: V_ED = " << r.end_dia_vol << " mL, V_ES = " << r.end_sys_vol << " mL, VS = "
       << pw_.stroke_volume_RV() << " mL, FE = " << pw_.ejection_fraction_RV() << " %" << endl;
}

void PiecewiseCirculation::save_state(State & st) const
{
  st["lv"]       = pack(pw_.state_LV());
  st["rv"]       = pack(pw_.state_RV());
  st["pressure"] = { p_cur_[0], p_cur_[1] };
  st["Jpv"]      = { Jpv_(0,0), Jpv_(1,0), Jpv_(0,1), Jpv_(1,1),
                     col_valid_[0] ? 1.0 : 0.0, col_valid_[1] ? 1.0 : 0.0 };
  st["counters"] = { double(solves_), double(beat_reported_), double(failures_) };
}

void PiecewiseCirculation::restore_state(const State & st)
{
  pw_.set_states(unpack(state_entry(st, "lv")), unpack(state_entry(st, "rv")));

  const std::vector<double> & p = state_entry(st, "pressure");
  p_cur_[0] = p.at(0);
  p_cur_[1] = p.at(1);

  const std::vector<double> & J = state_entry(st, "Jpv");
  Jpv_(0,0) = J.at(0); Jpv_(1,0) = J.at(1); Jpv_(0,1) = J.at(2); Jpv_(1,1) = J.at(3);
  col_valid_[0] = (J.at(4) != 0.0);
  col_valid_[1] = (J.at(5) != 0.0);

  const std::vector<double> & c = state_entry(st, "counters");
  solves_        = static_cast<int>(c.at(0));
  beat_reported_ = static_cast<int>(c.at(1));
  failures_      = static_cast<int>(c.at(2));
}
