#include "piecewise_circulation.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "util/command_line_args.h"

using std::cout;
using std::endl;

void PiecewiseCirculation::configure(const std::string & mshfile, const TimeBase & tb)
{
  base_ = mshfile;

  pw_.configure(mshfile);
  cout << "\n" << pw_.describe() << endl;

  double dt_ms = CommandLineArgs::read("-dtp", 1.0,
      "passo do controlador do PV loop por fases, em ms", "pvloop por fases");
  if (dt_ms <= 0.0) dt_ms = 1.0;

  // A rodada comeca com a fase 0 (t_init) e segue com `beats` periodos RR
  plan_.beats             = limit_beats(tb, tb.beats, pw_.timing().RR, pw_.timing().t_init);
  plan_.period_ms         = pw_.timing().RR;
  plan_.run_ms            = pw_.total_time(plan_.beats);
  plan_.step_ms           = dt_ms;
  plan_.nsteps            = static_cast<int>(std::round(plan_.run_ms / dt_ms));
  // A ativacao ocorre ao fim da fase 1: no artigo, "Phase 2: Isovolumetric
  // contraction. Electrical activation ensues". Sem este deslocamento ela
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

  cout << " Volumes iniciais da malha: VE = " << v_lv0
       << " mL, VD = " << v_rv0 << " mL" << endl;
  cout << " Conversao de pressao: 1 kPa = " << mech_->pressure_unit_per_kPa()
       << " unidades do modelo" << endl;

  if (v_rv0 > v_lv0)
    cout << " *** NOTA: o volume da cavidade DIREITA excede o da esquerda"
         << " nesta malha. Confirme os marcadores (30 = VE, 20 = VD) antes"
         << " de calibrar os limiares de ejecao." << endl;

  open_history(base_ + "_piecewise.dat", PiecewisePVLoop::history_header());

  cout << " Rodada: " << plan_.run_ms << " ms = " << plan_.nsteps
       << " passos de " << plan_.step_ms << " ms (" << plan_.beats
       << " batimento(s) de " << plan_.period_ms << " ms apos a fase 0)" << endl;
  cout << "\n Legenda: [PW] controlador de fases | [MECH] mecanica 3D" << endl;
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
  // ---- controlador de fases: p_VE/p_VD a partir de V e dV/dt do passo anterior
  pw_.advance_pressure(s.t_ms, s.dt_ms);

  const double p_lv_kPa = pw_.pressure_LV();
  const double p_rv_kPa = pw_.pressure_RV();

  if (s.log)
    cout << "  [PW   ] p_VE = " << p_lv_kPa << " kPa, p_VD = " << p_rv_kPa
         << " kPa  (P_art VE = " << pw_.state_LV().P_art << " kPa)" << endl;

  // ---- UMA solucao mecanica ------------------------------------------
  // A conversao kPa -> unidade de pressao do modelo e a MESMA que a tensao
  // ativa recebe, de proposito: um fator de 1000 entre a pressao da cavidade
  // e a tensao ativa seria quase invisivel nos resultados e fatal na fisiologia.
  const double kPa_to_p = mech_->pressure_unit_per_kPa();
  mech_->apply_pressures(s.t_ms, p_lv_kPa * kPa_to_p, p_rv_kPa * kPa_to_p);
  solves_++;

  // ---- volumes de volta ao controlador -------------------------------
  const double v_lv = mech_->volume_LV();
  const double v_rv = mech_->volume_RV();
  pw_.update_volumes(s.t_ms, s.dt_ms, v_lv, v_rv);

  if (s.log)
    cout << "  [MECH ] V_VE = " << v_lv << " mL (dV/dt = "
         << pw_.state_LV().Vdot << " mL/ms)"
         << " | V_VD = " << v_rv << " mL (dV/dt = "
         << pw_.state_RV().Vdot << " mL/ms)" << endl;

  if (history_.is_open())
  {
    history_ << pw_.history_row(s.t_ms) << "\n";
    history_.flush();
  }
}

void PiecewiseCirculation::end()
{
  if (history_.is_open()) history_.close();

  cout << "\n Total de solucoes 3D: " << solves_
       << " (" << (plan_.nsteps > 0 ? double(solves_) / plan_.nsteps : 0.0)
       << " por passo)" << endl;

  const PiecewisePVLoop::ChamberState & slv = pw_.state_LV();
  const PiecewisePVLoop::ChamberState & srv = pw_.state_RV();

  cout << " VE: V_ED = " << slv.V_ED << " mL, V_ES = " << slv.V_ES
       << " mL, VS = " << pw_.stroke_volume_LV() << " mL, FE = "
       << pw_.ejection_fraction_LV() << " %" << endl;
  cout << " VD: V_ED = " << srv.V_ED << " mL, V_ES = " << srv.V_ES
       << " mL, VS = " << pw_.stroke_volume_RV() << " mL, FE = "
       << pw_.ejection_fraction_RV() << " %" << endl;

  // O teto de |dP| e uma rede de seguranca. Se atuou muito, o ganho das
  // penalidades esta alto demais para este dt, e a fase isovolumetrica esta
  // sendo suavizada de um jeito que a calibracao nao ve.
  if (slv.n_clipped > 0 || srv.n_clipped > 0)
    cout << " *** ATENCAO: o teto -pwdpmax atuou " << slv.n_clipped
         << " vez(es) no VE e " << srv.n_clipped << " no VD. Reduza"
         << " -pwcpiso/-pwcviso ou -dtp antes de interpretar o laco." << endl;

  if (slv.V_ES >= slv.V_ED)
    cout << " *** ATENCAO: V_ES >= V_ED no VE: nao houve ejecao liquida."
         << " Verifique -tascale, o limiar -pwpej e se a ativacao caiu no"
         << " fim da fase 1." << endl;
}