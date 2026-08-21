#include "regazzoni2020.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>

Regazzoni2020::Regazzoni2020()
  : clamped(false),
    has_coupling(false),
    p_la(0), p_lv(0), p_ra(0), p_rv(0),
    q_mv(0), q_av(0), q_tv(0), q_pv(0)
{
  tim = default_timing();
  par = default_parameters();
  y   = default_initial_state();

  // Os TC/TR/tC de default_parameters() sao placeholders; a fonte da verdade
  // e a especificacao tim. Derivar aqui garante que qualquer objeto recem
  // construido ja esteja consistente, mesmo que ninguem chame set_heart_rate.
  apply_timing();
}

Regazzoni2020::Timing Regazzoni2020::default_timing()
{
  Timing t;
  t.t_act        = 0.00;   // solve_circulation() estimula em t = 0, T, 2T...
  t.emd          = 0.04;   // ToRORd-Land: estimulo -> inicio da subida de Ta
  t.TC_atria     = 0.17;   // sistole atrial, ~independente da frequencia
  t.TR_atria     = 0.17;
  t.TC_vent_ref  = 0.25;   // valores publicados, no periodo de referencia
  t.TR_vent_ref  = 0.40;
  t.RR_ref       = 1.00;   // os tempos publicados sao consistentes com 1.0 s
  t.sys_exponent = 0.50;   // Bazett
  t.atrial_max_fill = 0.90;
  return t;
}

Regazzoni2020::Parameters Regazzoni2020::default_parameters()
{
  Parameters p;

  // 1.0 Hz. O valor anterior (1.25) era incompativel com os tC publicados:
  // com RR = 0.8 s, tC_LA = 0.90 cai em 0.90 mod 0.8 = 0.10, exatamente o
  // tC do ventriculo, e o atrio passava a contrair dentro da sistole.
  p.HR = 1.0;

  //          EA      EB      TC    TR    tC    V0
  // TC/TR/tC abaixo sao apenas placeholders: apply_timing() os sobrescreve.
  p.LA = {  0.07,   0.18,  0.17, 0.17, 0.87,  4.0 };
  p.LV = {  4.482,  0.170, 0.25, 0.40, 0.04, 42.0 };
  p.RA = {  0.06,   0.07,  0.17, 0.17, 0.87,  4.0 };
  p.RV = {  0.200,  0.029, 0.25, 0.40, 0.04, 16.0 };

  const Valve v = { 0.0075, 75006.2 };
  p.MV = v;  p.AV = v;  p.TV = v;  p.PV = v;

  // systemic
  p.SYS.R_AR  = 0.733;   p.SYS.C_AR  = 1.372;   p.SYS.L_AR  = 5.0e-3;
  p.SYS.R_VEN = 0.320;   p.SYS.C_VEN = 11.363;  p.SYS.L_VEN = 5.0e-4;
  // pulmonary
  p.PUL.R_AR  = 0.046;   p.PUL.C_AR  = 20.0;    p.PUL.L_AR  = 5.0e-4;
  p.PUL.R_VEN = 0.0015;  p.PUL.C_VEN = 16.0;    p.PUL.L_VEN = 5.0e-4;

  return p;
}

std::array<double,Regazzoni2020::NSTATE> Regazzoni2020::default_initial_state()
{
  std::array<double,NSTATE> s;
  s[V_LA]      =  87.183;   // mL
  s[V_LV]      = 118.520;
  s[V_RA]      =  86.833;
  s[V_RV]      = 166.177;
  s[p_AR_SYS]  =  87.675;   // mmHg
  s[p_VEN_SYS] =  35.898;
  s[p_AR_PUL]  =  19.545;
  s[p_VEN_PUL] =  15.004;
  s[Q_AR_SYS]  =  71.104;   // mL/s
  s[Q_VEN_SYS] =  94.039;
  s[Q_AR_PUL]  =  94.084;
  s[Q_VEN_PUL] = 473.279;
  return s;
}

// ===========================================================================
//  Tempos das camaras derivados da frequencia
// ===========================================================================

double Regazzoni2020::wrap(double x, double p)
{
  double r = std::fmod(x, p);
  if (r < 0.0) r += p;
  return r;
}

void Regazzoni2020::apply_timing()
{
  const double RR = period();
  clamped = false;

  // ---- ventriculos ---------------------------------------------------
  // Escalonamento sublinear da sistole. Com sys_exponent = 0.5 (Bazett) a
  // duracao sistolica varia com a raiz do periodo, que e o que a regressao
  // de Weissler para QS2 produz aproximadamente entre 40 e 120 bpm.
  const double s = std::pow(RR / tim.RR_ref, tim.sys_exponent);

  double TC_v = tim.TC_vent_ref * s;
  double TR_v = tim.TR_vent_ref * s;

  // A ativacao mecanica do ventriculo comeca um EMD depois do estimulo. No
  // caminho acoplado a elastancia ventricular nao e usada, mas manter os dois
  // caminhos na mesma fase e o que permite usar o modo 0D puro como
  // diagnostico comparavel.
  const double t_vent = wrap(tim.t_act + tim.emd, RR);

  par.LV.TC = par.RV.TC = TC_v;
  par.LV.TR = par.RV.TR = TR_v;
  par.LV.tC = par.RV.tC = t_vent;

  // ---- atrios --------------------------------------------------------
  // A duracao da sistole atrial e tratada como independente da frequencia.
  // O que a frequencia muda e a janela diastolica disponivel; se ela ficar
  // pequena demais, a sistole atrial e encurtada em vez de invadir a sistole
  // ventricular.
  const double diast = RR - (TC_v + TR_v);

  double TC_a = tim.TC_atria;
  double TR_a = tim.TR_atria;

  if (diast <= 0.0)
  {
    // Sistole ventricular maior que o proprio ciclo: HR alto demais para a
    // duracao sistolica configurada. Deixa uma sistole atrial minima e
    // sinaliza; quem chamou decide o que fazer.
    TC_a = 0.05 * RR;
    TR_a = 0.05 * RR;
    clamped = true;
  }
  else if (TC_a > tim.atrial_max_fill * diast)
  {
    const double f = (tim.atrial_max_fill * diast) / TC_a;
    TC_a *= f;
    TR_a *= f;
    clamped = true;
  }

  // ANCORAGEM PELO FIM: a contracao atrial termina exatamente quando a
  // pressao ventricular comeca a subir (t_act + emd). Esta e a linha que
  // torna a sincronia AV automatica em qualquer frequencia -- tC atrial NAO
  // e uma fracao fixa do ciclo.
  const double t_atria = wrap(tim.t_act + tim.emd - TC_a, RR);

  par.LA.TC = par.RA.TC = TC_a;
  par.LA.TR = par.RA.TR = TR_a;
  par.LA.tC = par.RA.tC = t_atria;
}

void Regazzoni2020::set_heart_rate(double hr_hz)
{
  if (hr_hz <= 0.0) return;
  par.HR = hr_hz;
  apply_timing();
}

void Regazzoni2020::set_activation_time(double t_act)
{
  tim.t_act = t_act;
  apply_timing();
}

void Regazzoni2020::set_emd(double emd)
{
  tim.emd = emd;
  apply_timing();
}

double Regazzoni2020::diastolic_window() const
{
  return period() - (par.LV.TC + par.LV.TR);
}

std::string Regazzoni2020::describe_timing() const
{
  const double RR = period();
  std::ostringstream os;
  os << std::fixed << std::setprecision(3);

  os << " Cardiac timing: HR = " << par.HR * 60.0 << " bpm, RR = " << RR
     << " s\n";
  os << "   ventricular activation (stimulus) at t = " << tim.t_act
     << " s, EMD = " << tim.emd << " s\n";
  os << "   atria      : contract [" << par.LA.tC << ", "
     << wrap(par.LA.tC + par.LA.TC, RR) << "), relax "
     << par.LA.TR << " s"
     << (par.LA.tC + par.LA.TC > RR ? "   (wraps)" : "") << "\n";
  os << "   ventricles : contract [" << par.LV.tC << ", "
     << wrap(par.LV.tC + par.LV.TC, RR) << "), relax "
     << par.LV.TR << " s   (0D only; ignored when coupled)\n";
  os << "   diastolic window = " << diastolic_window() << " s ("
     << 100.0 * diastolic_window() / RR << " % of RR)";

  if (clamped)
    os << "\n   *** WARNING: the atrial systole had to be shortened to fit"
          " the cycle.\n"
          "       The check uses the 0D ventricular window (TC+TR). When the"
          " 3D model\n"
          "       supplies p_LV/p_RV that window is only an estimate: what"
          " matters is\n"
          "       whether Ta has actually relaxed before t = "
       << std::fixed << std::setprecision(3) << par.LA.tC
       << " s. Raise sys_exponent,\n"
          "       shorten TC_vent_ref/TR_vent_ref, or lower the heart rate.";

  return os.str();
}

// ===========================================================================
//  Volume
// ===========================================================================

double Regazzoni2020::rebalance_total_volume(double vtot)
{
  const double dv = vtot - total_volume();
  const double C  = par.SYS.C_VEN + par.PUL.C_VEN;
  if (C <= 0.0) return 0.0;

  // Deslocamento comum de pressao nos dois reservatorios venosos: eles sao os
  // vasos de capacitancia e absorvem a diferenca sem alterar o gradiente
  // arterio-venoso que define o ponto de operacao.
  const double dp = dv / C;
  y[p_VEN_SYS] += dp;
  y[p_VEN_PUL] += dp;

  return dv;
}

// ===========================================================================

double Regazzoni2020::activation(const Chamber & c, double t) const
{
  // Blanco & Feijoo time-varying elastance:
  //   contraction  0 <= (t-tC) mod RR < TC  -> 0.5(1 - cos(pi/TC (t-tC)))
  //   relaxation   0 <= (t-tR) mod RR < TR  -> 0.5(1 + cos(pi/TR (t-tR)))
  //   otherwise    rest -> 0
  const double RR = period();
  const double tR = c.tC + c.TC;

  double e = 0.0;

  double s1 = std::fmod(t - c.tC, RR);
  if (s1 < 0.0) s1 += RR;
  if (s1 < c.TC)
    e += 0.5 * (1.0 - std::cos(M_PI / c.TC * s1));

  double s2 = std::fmod(t - tR, RR);
  if (s2 < 0.0) s2 += RR;
  if (s2 < c.TR)
    e += 0.5 * (1.0 + std::cos(M_PI / c.TR * s2));

  return std::min(std::max(e, 0.0), 1.0);
}

double Regazzoni2020::chamber_pressure(const Chamber & c, double V, double t) const
{
  double E = c.EA * activation(c, t) + c.EB;
  return E * (V - c.V0);
}

double Regazzoni2020::valve_resistance(const Valve & v, double p1, double p2) const
{
  // Smooth transition between Rmin (open, p1 > p2) and Rmax (closed).
  // Interpolated in log10 so the resistance spans orders of magnitude
  // without introducing a discontinuity that would break the ODE solver.
  double x  = p2 - p1;
  double sh = std::atan(M_PI / 2.0 * x * 200.0) / M_PI + 0.5;  // smooth step
  double e  = std::log10(v.Rmin)
              + (std::log10(v.Rmax) - std::log10(v.Rmin)) * sh;
  return std::pow(10.0, e);
}

void Regazzoni2020::rhs(double t, const std::array<double,NSTATE> & yy,
                        std::array<double,NSTATE> & dy)
{
  // ---- chamber pressures -------------------------------------------
  p_la = chamber_pressure(par.LA, yy[V_LA], t);
  p_ra = chamber_pressure(par.RA, yy[V_RA], t);

  if (has_coupling)
  {
    // ventricular pressures come from the 3D model
    coupling(yy[V_LV], yy[V_RV], t, p_lv, p_rv);
  }
  else
  {
    p_lv = chamber_pressure(par.LV, yy[V_LV], t);
    p_rv = chamber_pressure(par.RV, yy[V_RV], t);
  }

  // ---- valve flows --------------------------------------------------
  q_mv = flux_through_valve(par.MV, p_la, p_lv);              // LA -> LV
  q_av = flux_through_valve(par.AV, p_lv, yy[p_AR_SYS]);      // LV -> aorta
  q_tv = flux_through_valve(par.TV, p_ra, p_rv);              // RA -> RV
  q_pv = flux_through_valve(par.PV, p_rv, yy[p_AR_PUL]);      // RV -> pulm

  // ---- chamber volumes ----------------------------------------------
  dy[V_LA] = yy[Q_VEN_PUL] - q_mv;
  dy[V_LV] = q_mv - q_av;
  dy[V_RA] = yy[Q_VEN_SYS] - q_tv;
  dy[V_RV] = q_tv - q_pv;

  // ---- compartment pressures (compliance) ---------------------------
  dy[p_AR_SYS]  = (q_av         - yy[Q_AR_SYS] ) / par.SYS.C_AR;
  dy[p_VEN_SYS] = (yy[Q_AR_SYS] - yy[Q_VEN_SYS]) / par.SYS.C_VEN;
  dy[p_AR_PUL]  = (q_pv         - yy[Q_AR_PUL] ) / par.PUL.C_AR;
  dy[p_VEN_PUL] = (yy[Q_AR_PUL] - yy[Q_VEN_PUL]) / par.PUL.C_VEN;

  // ---- compartment flows (inertance + resistance) -------------------
  dy[Q_AR_SYS]  = -(par.SYS.R_AR  * yy[Q_AR_SYS]
                    + yy[p_VEN_SYS] - yy[p_AR_SYS] ) / par.SYS.L_AR;
  dy[Q_VEN_SYS] = -(par.SYS.R_VEN * yy[Q_VEN_SYS]
                    + p_ra          - yy[p_VEN_SYS]) / par.SYS.L_VEN;
  dy[Q_AR_PUL]  = -(par.PUL.R_AR  * yy[Q_AR_PUL]
                    + yy[p_VEN_PUL] - yy[p_AR_PUL] ) / par.PUL.L_AR;
  dy[Q_VEN_PUL] = -(par.PUL.R_VEN * yy[Q_VEN_PUL]
                    + p_la          - yy[p_VEN_PUL]) / par.PUL.L_VEN;
}

void Regazzoni2020::step(double t, double dt)
{
  std::array<double,NSTATE> dy;
  rhs(t, y, dy);
  for (int i = 0; i < NSTATE; i++)
    y[i] += dt * dy[i];
}

double Regazzoni2020::total_volume() const
{
  return y[V_LA] + y[V_LV] + y[V_RA] + y[V_RV]
       + par.SYS.C_AR  * y[p_AR_SYS]
       + par.SYS.C_VEN * y[p_VEN_SYS]
       + par.PUL.C_AR  * y[p_AR_PUL]
       + par.PUL.C_VEN * y[p_VEN_PUL];
}

std::string Regazzoni2020::history_header()
{
  return "# t V_LV p_LV V_RV p_RV V_LA p_LA V_RA p_RA "
         "p_AR_SYS p_VEN_SYS p_AR_PUL p_VEN_PUL "
         "Q_MV Q_AV Q_TV Q_PV V_tot";
}

std::string Regazzoni2020::history_row(double t) const
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(6);
  os << t << " "
     << y[V_LV] << " " << p_lv << " "
     << y[V_RV] << " " << p_rv << " "
     << y[V_LA] << " " << p_la << " "
     << y[V_RA] << " " << p_ra << " "
     << y[p_AR_SYS] << " " << y[p_VEN_SYS] << " "
     << y[p_AR_PUL] << " " << y[p_VEN_PUL] << " "
     << q_mv << " " << q_av << " " << q_tv << " " << q_pv << " "
     << total_volume();
  return os.str();
}
