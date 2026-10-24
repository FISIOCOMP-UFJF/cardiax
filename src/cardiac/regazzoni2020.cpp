#include "regazzoni2020.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>

Regazzoni2020::Regazzoni2020()
  : clamped(false), has_coupling(false),
    p_la(0), p_lv(0), p_ra(0), p_rv(0),
    q_mv(0), q_av(0), q_tv(0), q_pv(0)
{
  tim = default_timing();
  par = default_parameters();
  y   = default_initial_state();
  apply_timing();      // os TC/TR/tC de default_parameters() sao placeholders
}

Regazzoni2020::Timing Regazzoni2020::default_timing()
{
  Timing t;
  t.t_act           = 0.0;     // estimulo em t = 0, RR, 2RR...
  t.emd             = 40.0;    // ToRORd-Land: estimulo -> subida de Ta
  t.TC_atria        = 170.0;
  t.TR_atria        = 170.0;
  t.TC_vent_ref     = 250.0;   // publicados, validos em RR_ref
  t.TR_vent_ref     = 400.0;
  t.RR_ref          = 1000.0;
  t.sys_exponent    = 0.5;     // Bazett
  t.atrial_max_fill = 0.9;
  return t;
}

// Valores de Regazzoni et al. (mmHg, s, mL) ja convertidos para Pa, ms, mL.
Regazzoni2020::Parameters Regazzoni2020::default_parameters()
{
  Parameters p;
  p.RR = 1000.0;

  //          EA            EB            TC     TR     tC     V0
  // TC/TR/tC sao placeholders: apply_timing() os sobrescreve.
  p.LA = {   9.3325671190505,  23.9980297347,  170.0, 170.0, 870.0,  4.0 };
  p.LV = { 597.55094039403,    22.66480586055,  250.0, 400.0,  40.0, 42.0 };
  p.RA = {   7.9993432449,   9.3325671190505, 170.0, 170.0, 870.0,  4.0 };
  p.RV = {  26.664477483,    3.866349235035, 250.0, 400.0,  40.0, 16.0 };

  const Valve v = { 999.9179056125, 10000005654.927 };
  p.MV = v;  p.AV = v;  p.TV = v;  p.PV = v;

  //                     R_AR         C_AR          L_AR        R_VEN        C_VEN         L_VEN
  p.SYS = { 97725.309975195,  0.0102908448206024, 666611.937075, 42663.1639728, 0.0852294968633419, 66661.1937075 };
  p.PUL = {  6132.82982109, 0.150012315169131,   66661.1937075,  199.9835811225, 0.120009852135305,  66661.1937075 };

  return p;
}

std::array<double,Regazzoni2020::NSTATE> Regazzoni2020::default_initial_state()
{
  std::array<double,NSTATE> s;
  s[V_LA]      =  87.183;      // mL
  s[V_LV]      = 118.520;
  s[V_RA]      =  86.833;
  s[V_RV]      = 166.177;
  s[p_AR_SYS]  = 11689.0403166101;  // Pa
  s[p_VEN_SYS] =  4786.00706342367;
  s[p_AR_PUL]  =  2605.78606202618;
  s[p_VEN_PUL] =  2000.36910077466;
  s[Q_AR_SYS]  = 0.071104;     // mL/ms
  s[Q_VEN_SYS] = 0.094039;
  s[Q_AR_PUL]  = 0.094084;
  s[Q_VEN_PUL] = 0.473279;
  return s;
}

// ===========================================================================
//  Tempos das camaras derivados do periodo
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

  // Ventriculos: sistole sublinear com o periodo
  const double s = std::pow(RR / tim.RR_ref, tim.sys_exponent);
  const double TC_v = tim.TC_vent_ref * s;
  const double TR_v = tim.TR_vent_ref * s;
  const double t_vent = wrap(tim.t_act + tim.emd, RR);

  par.LV.TC = par.RV.TC = TC_v;
  par.LV.TR = par.RV.TR = TR_v;
  par.LV.tC = par.RV.tC = t_vent;

  // Atrios: duracao independente do periodo, mas encurtada se a diastole nao
  // comportar (HR alto), em vez de invadir a sistole ventricular
  const double diast = RR - (TC_v + TR_v);
  double TC_a = tim.TC_atria;
  double TR_a = tim.TR_atria;

  if (diast <= 0.0)
  {
    TC_a = TR_a = 0.05 * RR;
    clamped = true;
  }
  else if (TC_a > tim.atrial_max_fill * diast)
  {
    const double f = (tim.atrial_max_fill * diast) / TC_a;
    TC_a *= f;
    TR_a *= f;
    clamped = true;
  }

  // Ancoragem pelo FIM: a contracao atrial termina quando a pressao
  // ventricular comeca a subir
  const double t_atria = wrap(tim.t_act + tim.emd - TC_a, RR);

  par.LA.TC = par.RA.TC = TC_a;
  par.LA.TR = par.RA.TR = TR_a;
  par.LA.tC = par.RA.tC = t_atria;
}

void Regazzoni2020::set_period(double rr_ms)
{
  if (rr_ms <= 0.0) return;
  par.RR = rr_ms;
  apply_timing();
}

void Regazzoni2020::set_activation_time(double t_act_ms) { tim.t_act = t_act_ms; apply_timing(); }
void Regazzoni2020::set_emd(double emd_ms)               { tim.emd   = emd_ms;   apply_timing(); }

double Regazzoni2020::diastolic_window() const
{
  return period() - (par.LV.TC + par.LV.TR);
}

std::string Regazzoni2020::describe_timing() const
{
  const double RR = period();
  std::ostringstream os;
  os << std::fixed << std::setprecision(1);

  os << " Cardiac timing: RR = " << RR << " ms (" << 60000.0 / RR << " bpm)\n";
  os << "   ventricular activation at t = " << tim.t_act << " ms, EMD = "
     << tim.emd << " ms\n";
  os << "   atria      : contract [" << par.LA.tC << ", "
     << wrap(par.LA.tC + par.LA.TC, RR) << ") ms, relax " << par.LA.TR << " ms"
     << (par.LA.tC + par.LA.TC > RR ? "   (wraps)" : "") << "\n";
  os << "   ventricles : contract [" << par.LV.tC << ", "
     << wrap(par.LV.tC + par.LV.TC, RR) << ") ms, relax " << par.LV.TR
     << " ms   (0D only; ignored when coupled)\n";
  os << "   diastolic window = " << diastolic_window() << " ms ("
     << 100.0 * diastolic_window() / RR << " % of RR)";

  if (clamped)
    os << "\n   *** WARNING: the atrial systole had to be shortened to fit the"
          " cycle. When the 3D model supplies p_LV/p_RV the 0D window (TC+TR) is"
          " only an estimate: what matters is whether Ta has relaxed before t = "
       << par.LA.tC << " ms. Lower sys_exponent/TC_vent_ref/TR_vent_ref or"
          " increase -rr.";
  return os.str();
}

double Regazzoni2020::rebalance_total_volume(double vtot)
{
  const double dv = vtot - total_volume();
  const double C  = par.SYS.C_VEN + par.PUL.C_VEN;
  if (C <= 0.0) return 0.0;

  const double dp = dv / C;          // [Pa]: mL / (mL/Pa)
  y[p_VEN_SYS] += dp;
  y[p_VEN_PUL] += dp;
  return dv;
}

// ===========================================================================

double Regazzoni2020::activation(const Chamber & c, double t) const
{
  // Blanco & Feijoo: contracao 0.5(1-cos), relaxamento 0.5(1+cos), repouso 0
  const double RR = period();
  const double tR = c.tC + c.TC;

  double e = 0.0;
  double s1 = std::fmod(t - c.tC, RR);
  if (s1 < 0.0) s1 += RR;
  if (s1 < c.TC) e += 0.5 * (1.0 - std::cos(M_PI / c.TC * s1));

  double s2 = std::fmod(t - tR, RR);
  if (s2 < 0.0) s2 += RR;
  if (s2 < c.TR) e += 0.5 * (1.0 + std::cos(M_PI / c.TR * s2));

  return std::min(std::max(e, 0.0), 1.0);
}

double Regazzoni2020::chamber_pressure(const Chamber & c, double V, double t) const
{
  return (c.EA * activation(c, t) + c.EB) * (V - c.V0);
}

double Regazzoni2020::valve_resistance(const Valve & v, double p1, double p2) const
{
  // Transicao suave entre Rmin (aberta, p1 > p2) e Rmax (fechada), em log10.
  // 1.50012315169131 /Pa = 200 /mmHg, a nitidez do degrau do modelo publicado.
  const double sharp = 1.50012315169131;
  const double x  = p2 - p1;
  const double sh = std::atan(M_PI / 2.0 * x * sharp) / M_PI + 0.5;
  const double e  = std::log10(v.Rmin) + (std::log10(v.Rmax) - std::log10(v.Rmin)) * sh;
  return std::pow(10.0, e);
}

void Regazzoni2020::rhs(double t, const std::array<double,NSTATE> & yy,
                        std::array<double,NSTATE> & dy)
{
  p_la = chamber_pressure(par.LA, yy[V_LA], t);
  p_ra = chamber_pressure(par.RA, yy[V_RA], t);

  if (has_coupling)
    coupling(yy[V_LV], yy[V_RV], t, p_lv, p_rv);   // pressoes vem do 3D
  else
  {
    p_lv = chamber_pressure(par.LV, yy[V_LV], t);
    p_rv = chamber_pressure(par.RV, yy[V_RV], t);
  }

  q_mv = flux_through_valve(par.MV, p_la, p_lv);              // LA -> LV
  q_av = flux_through_valve(par.AV, p_lv, yy[p_AR_SYS]);      // LV -> aorta
  q_tv = flux_through_valve(par.TV, p_ra, p_rv);              // RA -> RV
  q_pv = flux_through_valve(par.PV, p_rv, yy[p_AR_PUL]);      // RV -> pulmonar

  dy[V_LA] = yy[Q_VEN_PUL] - q_mv;
  dy[V_LV] = q_mv - q_av;
  dy[V_RA] = yy[Q_VEN_SYS] - q_tv;
  dy[V_RV] = q_tv - q_pv;

  dy[p_AR_SYS]  = (q_av          - yy[Q_AR_SYS] ) / par.SYS.C_AR;
  dy[p_VEN_SYS] = (yy[Q_AR_SYS]  - yy[Q_VEN_SYS]) / par.SYS.C_VEN;
  dy[p_AR_PUL]  = (q_pv          - yy[Q_AR_PUL] ) / par.PUL.C_AR;
  dy[p_VEN_PUL] = (yy[Q_AR_PUL]  - yy[Q_VEN_PUL]) / par.PUL.C_VEN;

  dy[Q_AR_SYS]  = -(par.SYS.R_AR  * yy[Q_AR_SYS]  + yy[p_VEN_SYS] - yy[p_AR_SYS] ) / par.SYS.L_AR;
  dy[Q_VEN_SYS] = -(par.SYS.R_VEN * yy[Q_VEN_SYS] + p_ra          - yy[p_VEN_SYS]) / par.SYS.L_VEN;
  dy[Q_AR_PUL]  = -(par.PUL.R_AR  * yy[Q_AR_PUL]  + yy[p_VEN_PUL] - yy[p_AR_PUL] ) / par.PUL.L_AR;
  dy[Q_VEN_PUL] = -(par.PUL.R_VEN * yy[Q_VEN_PUL] + p_la          - yy[p_VEN_PUL]) / par.PUL.L_VEN;
}

void Regazzoni2020::step(double t, double dt)
{
  std::array<double,NSTATE> dy;
  rhs(t, y, dy);
  for (int i = 0; i < NSTATE; i++) y[i] += dt * dy[i];
}

double Regazzoni2020::total_volume() const
{
  return y[V_LA] + y[V_LV] + y[V_RA] + y[V_RV]
       + par.SYS.C_AR  * y[p_AR_SYS] + par.SYS.C_VEN * y[p_VEN_SYS]
       + par.PUL.C_AR  * y[p_AR_PUL] + par.PUL.C_VEN * y[p_VEN_PUL];
}

std::string Regazzoni2020::history_header()
{
  return "# t[ms] V_LV[mL] p_LV[Pa] V_RV p_RV V_LA p_LA V_RA p_RA "
         "p_AR_SYS p_VEN_SYS p_AR_PUL p_VEN_PUL "
         "Q_MV[mL/ms] Q_AV Q_TV Q_PV V_tot[mL]";
}

std::string Regazzoni2020::history_row(double t) const
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(6);
  os << t << " "
     << y[V_LV] << " " << p_lv << " " << y[V_RV] << " " << p_rv << " "
     << y[V_LA] << " " << p_la << " " << y[V_RA] << " " << p_ra << " "
     << y[p_AR_SYS] << " " << y[p_VEN_SYS] << " "
     << y[p_AR_PUL] << " " << y[p_VEN_PUL] << " "
     << q_mv << " " << q_av << " " << q_tv << " " << q_pv << " "
     << total_volume();
  return os.str();
}