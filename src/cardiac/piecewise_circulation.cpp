#include "piecewise_circulation.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using std::cout;
using std::endl;

namespace
{
typedef PiecewisePVLoop::ChamberState ChamberState;

// ChamberState <-> vetor, campo a campo, na ordem da declaracao
std::vector<double> pack(const ChamberState & s)
{
  return { double(s.phase), s.P, s.P_art, s.P_prev, s.V, s.V_prev, s.Vdot,
           s.V_lock, s.V_DS, s.V_ED, s.V_ES, s.V_ref, s.has_ejected ? 1.0 : 0.0,
           s.C_cav, s.t_phase, double(s.n_clipped) };
}

ChamberState unpack(const std::vector<double> & v)
{
  if (v.size() != 16)
    throw std::runtime_error("Checkpoint: estado de cavidade do -pw com tamanho errado");
  ChamberState s;
  s.phase  = static_cast<PiecewisePVLoop::Phase>(static_cast<int>(v[0]));
  s.P      = v[1];  s.P_art  = v[2];  s.P_prev = v[3];
  s.V      = v[4];  s.V_prev = v[5];  s.Vdot   = v[6];
  s.V_lock = v[7];  s.V_DS   = v[8];  s.V_ED   = v[9];
  s.V_ES   = v[10]; s.V_ref  = v[11];
  s.has_ejected = (v[12] != 0.0);
  s.C_cav  = v[13]; s.t_phase = v[14];
  s.n_clipped = static_cast<int>(v[15]);
  return s;
}
} // namespace

void PiecewiseCirculation::configure(const std::string & mshfile, const TimeBase & tb)
{
  base_ = mshfile;

  pw_.configure(mshfile, tb.period_ms);
  cout << "\n" << pw_.describe() << endl;

  // A rodada e a fase 0 (t_init) seguida de `beats` periodos RR
  plan_.beats             = limit_beats(tb, tb.beats, pw_.timing().RR, pw_.timing().t_init);
  plan_.period_ms         = pw_.timing().RR;
  plan_.run_ms            = pw_.total_time(plan_.beats);
  plan_.step_ms           = tb.step_ms;
  plan_.nsteps            = static_cast<int>(std::round(plan_.run_ms / tb.step_ms));
  // A ativacao ocorre ao fim da fase 1 (no artigo, "Phase 2: Isovolumetric
  // contraction. Electrical activation ensues"). Sem este deslocamento ela
  // cairia em t = 0, no meio da rampa de inflacao passiva.
  plan_.stim_offset_ms    = pw_.stimulus_offset();
  plan_.span_margin_steps = 2;
  plan_.default_out_every = 8;
}

void PiecewiseCirculation::do_begin()
{
  const double v_lv0 = mech_->volume_LV();
  const double v_rv0 = mech_->volume_RV();
  pw_.init(v_lv0, v_rv0);

  cout << " Volumes iniciais da malha: VE = " << v_lv0 << " mL, VD = " << v_rv0 << " mL" << endl;
  if (v_rv0 > v_lv0)
    cout << " *** NOTA: o volume da cavidade DIREITA excede o da esquerda nesta malha."
         << " Confirme os marcadores (30 = VE, 20 = VD) antes de calibrar a ejecao." << endl;

  open_history(base_ + "_piecewise.dat", PiecewisePVLoop::history_header());

  cout << " Rodada: " << plan_.run_ms << " ms = " << plan_.nsteps << " passos de "
       << plan_.step_ms << " ms (" << plan_.beats << " batimento(s) de "
       << plan_.period_ms << " ms apos a fase 0)\n"
       << "\n Legenda: [PW] controlador de fases | [MECH] mecanica 3D" << endl;
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

void PiecewiseCirculation::advance(const StepContext & s)
{
  // 1. pressoes do passo, a partir de V e dV/dt do passo anterior
  pw_.advance_pressure(s.t_ms, s.dt_ms);

  if (s.log)
    cout << "  [PW   ] p_VE = " << pw_.pressure_LV() << " Pa, p_VD = " << pw_.pressure_RV()
         << " Pa  (P_art VE = " << pw_.state_LV().P_art << " Pa)" << endl;

  // 2. UMA solucao mecanica
  mech_->apply_pressures(s.t_ms, pw_.pressure_LV(), pw_.pressure_RV());
  solves_++;

  // 3. volumes de volta ao controlador (decidem as transicoes de fase)
  const double v_lv = mech_->volume_LV();
  const double v_rv = mech_->volume_RV();
  pw_.update_volumes(s.t_ms, s.dt_ms, v_lv, v_rv);

  if (s.log)
    cout << "  [MECH ] V_VE = " << v_lv << " mL (dV/dt = " << pw_.state_LV().Vdot
         << " mL/ms) | V_VD = " << v_rv << " mL (dV/dt = " << pw_.state_RV().Vdot
         << " mL/ms)" << endl;

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

  const PiecewisePVLoop::ChamberState & slv = pw_.state_LV();
  const PiecewisePVLoop::ChamberState & srv = pw_.state_RV();

  cout << " VE: V_ED = " << slv.V_ED << " mL, V_ES = " << slv.V_ES << " mL, VS = "
       << pw_.stroke_volume_LV() << " mL, FE = " << pw_.ejection_fraction_LV() << " %\n"
       << " VD: V_ED = " << srv.V_ED << " mL, V_ES = " << srv.V_ES << " mL, VS = "
       << pw_.stroke_volume_RV() << " mL, FE = " << pw_.ejection_fraction_RV() << " %" << endl;

  // O teto de |dP| e uma rede de seguranca: se atuou muito, o ganho das
  // penalidades esta alto demais para este passo.
  if (slv.n_clipped > 0 || srv.n_clipped > 0)
    cout << " *** ATENCAO: o teto -pwdpmax atuou " << slv.n_clipped << " vez(es) no VE e "
         << srv.n_clipped << " no VD. Reduza -pwcpiso/-pwcviso ou -dt_mech antes de"
         << " interpretar o laco." << endl;

  if (slv.V_ES >= slv.V_ED)
    cout << " *** ATENCAO: V_ES >= V_ED no VE: nao houve ejecao liquida. Verifique"
         << " -tascale, o limiar -pwpej e se a ativacao caiu no fim da fase 1." << endl;
}
void PiecewiseCirculation::save_state(State & st) const
{
  st["lv"]       = pack(pw_.state_LV());
  st["rv"]       = pack(pw_.state_RV());
  st["counters"] = { double(solves_), double(beat_reported_) };
}

void PiecewiseCirculation::restore_state(const State & st)
{
  pw_.set_states(unpack(state_entry(st, "lv")), unpack(state_entry(st, "rv")));

  const std::vector<double> & c = state_entry(st, "counters");
  solves_        = static_cast<int>(c.at(0));
  beat_reported_ = static_cast<int>(c.at(1));
}
