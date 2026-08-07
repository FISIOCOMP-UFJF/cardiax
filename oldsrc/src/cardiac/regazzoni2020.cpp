#include "regazzoni2020.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>

Regazzoni2020::Regazzoni2020()
  : has_coupling(false),
    p_la(0), p_lv(0), p_ra(0), p_rv(0),
    q_mv(0), q_av(0), q_tv(0), q_pv(0)
{
  par = default_parameters();
  y   = default_initial_state();
}

Regazzoni2020::Parameters Regazzoni2020::default_parameters()
{
  Parameters p;
  p.HR = 1.25;                     // Hz  -> 0.8 s period, 75 bpm

  //          EA      EB      TC    TR    tC    V0
  p.LA = {  0.07,   0.18,  0.17, 0.17, 0.90,  4.0 };
  p.LV = {  4.482,  0.170, 0.25, 0.40, 0.10, 42.0 };
  p.RA = {  0.06,   0.07,  0.17, 0.17, 0.90,  4.0 };
  p.RV = {  0.200,  0.029, 0.25, 0.40, 0.10, 16.0 };

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
