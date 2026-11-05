#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "piecewise_pvloop.hpp"
#include "util/command_line_args.h"
#include "util/pugixml.hpp"

// ======================================================================
//  Local helpers
// ======================================================================

namespace
{
  constexpr double kMmHg = 133.322387415;   // Pa

  // mmHg*s/mL -> Pa*ms/mL  and  mL/mmHg -> mL/Pa
  constexpr double res(double mmHg_s_per_mL) { return mmHg_s_per_mL * kMmHg * 1000.0; }
  constexpr double comp(double mL_per_mmHg)  { return mL_per_mmHg / kMmHg; }

  double xml_attr(const pugi::xml_node & node, const char * name, double def)
  {
    if (!node) return def;
    pugi::xml_attribute a = node.attribute(name);
    return a ? a.as_double() : def;
  }
}

// ======================================================================
//  Windkessel
// ======================================================================

// Implicit Euler with Q = -(V - V_n)/dt, D = 1 + dt/(R C):
//   P_c = (P_c,n + dt Q/C)/D,  P = P_c + Rc Q  =>  P = A + B V
void PiecewisePVLoop::Windkessel::ejection_law(double P_c, double V_n, double dt,
                                               double & A, double & B) const
{
  const double D = 1.0 + dt / (R * C);
  A = P_c / D + V_n / (C * D) + Rc * V_n / dt;
  B = -(1.0 / (C * D) + Rc / dt);
}

void PiecewisePVLoop::Windkessel::advance(double P_c, double V, double V_n, double dt,
                                          bool ejecting, bool has_ejected,
                                          double & P_c_new, double & Q) const
{
  // Valve closed: P_c is held before the first ejection, then drains through R
  const double D = 1.0 + dt / (R * C);
  if (ejecting)
  {
    Q = -(V - V_n) / dt;
    P_c_new = (P_c + dt * Q / C) / D;
  }
  else
  {
    Q = 0.0;
    P_c_new = has_ejected ? P_c / D : P_c;
  }
}

// ======================================================================
//  Construction and defaults (BiV)
// ======================================================================

PiecewisePVLoop::PiecewisePVLoop()
{
  par_lv.t_zero                = 50.0;
  par_lv.preload_pressure      = 500.0;
  par_lv.t_end_diastole        = 120.0;
  par_lv.p_end_diastole        = 1000.0;
  par_lv.p_fill                = 500.0;
  par_lv.period                = 800.0;
  par_lv.windkessel.p_init     = 9000.0;
  par_lv.windkessel.C          = comp(1.5);
  par_lv.windkessel.R          = res(1.1);
  par_lv.windkessel.Rc         = res(0.03);
  par_lv.fill_rate             = 0.104;
  par_lv.min_ejection_duration = 10.0;
  par_lv.dvol_eps              = 1e-4;

  par_rv = par_lv;
  par_rv.preload_pressure      = 200.0;
  par_rv.p_end_diastole        = 400.0;
  par_rv.p_fill                = 200.0;
  par_rv.windkessel.p_init     = 2500.0;
  par_rv.windkessel.C          = comp(4.0);
  par_rv.windkessel.R          = res(0.5);
  par_rv.windkessel.Rc         = res(0.01);
  par_rv.fill_rate             = 0.057;

  initialize(0.0, 0.0, 0.0, 0.0, 0.0);
}

const char * PiecewisePVLoop::phase_name(Phase p)
{
  switch (p)
  {
    case PRELOAD:  return "PRELOAD";
    case IVC:      return "ISOVOLUMIC_CONTRACTION";
    case EJECTION: return "EJECTION";
    case IVR:      return "ISOVOLUMIC_RELAXATION";
    case FILLING:  return "FILLING";
  }
  return "?";
}

// ======================================================================
//  Configuration: XML first, then the command line
// ======================================================================

void PiecewisePVLoop::configure(const std::string & mshfile, double rr_ms)
{
  // ---- 1. <pvloop_piecewise> block of the mesh, if present --------------
  // Times are attributes of the block itself; the rest go in <lv>/<rv>.
  pugi::xml_document doc;
  pugi::xml_node root;
  if (doc.load_file(mshfile.c_str()))
    root = doc.child("mesh").child("pvloop_piecewise");

  std::cout << (root ? " <pvloop_piecewise> encontrado na malha"
                     : " <pvloop_piecewise> ausente na malha; usando os defaults")
            << std::endl;

  ChamberParams * pp[2] = { &par_lv, &par_rv };
  const char * child[2] = { "lv", "rv" };
  for (int c = 0; c < 2; c++)
  {
    ChamberParams & p = *pp[c];
    pugi::xml_node n = root ? root.child(child[c]) : pugi::xml_node();

    p.period         = xml_attr(root, "period",         p.period);
    p.t_zero         = xml_attr(root, "t_zero",         p.t_zero);
    p.t_end_diastole = xml_attr(root, "t_end_diastole", p.t_end_diastole);
    p.min_ejection_duration = xml_attr(root, "min_ejection_duration", p.min_ejection_duration);
    p.dvol_eps       = xml_attr(root, "dvol_eps",       p.dvol_eps);

    p.preload_pressure = xml_attr(n, "preload_pressure", p.preload_pressure);
    p.p_end_diastole   = xml_attr(n, "p_end_diastole",   p.p_end_diastole);
    p.p_fill           = xml_attr(n, "p_fill",           p.p_fill);
    p.fill_rate        = xml_attr(n, "fill_rate",        p.fill_rate);
    p.windkessel.p_init = xml_attr(n, "p_init",          p.windkessel.p_init);
    p.windkessel.C     = xml_attr(n, "compliance",       p.windkessel.C);
    p.windkessel.R     = xml_attr(n, "resistance",       p.windkessel.R);
    p.windkessel.Rc    = xml_attr(n, "characteristic_impedance", p.windkessel.Rc);
  }

  // ---- 2. command line ----------------------------------------------
  const char * G = "pvloop por fases";
  if (rr_ms > 0.0) par_lv.period = rr_ms;

  const double t_zero = CommandLineArgs::read("-pwinit", par_lv.t_zero,
      "t_zero: fim da rampa 0 -> preload_pressure (1o batimento) e inicio do PRELOAD"
      " dos batimentos seguintes, ms", G);
  const double t_ed_dur = CommandLineArgs::read("-pwedp", par_lv.t_end_diastole - par_lv.t_zero,
      "duracao da rampa preload_pressure -> p_end_diastole: t_end_diastole = t_zero + isto, ms", G);
  const double min_ej = CommandLineArgs::read("-pwminphase", par_lv.min_ejection_duration,
      "min_ejection_duration: ejecao minima antes de testar o fechamento, ms", G);
  const double dveps = CommandLineArgs::read("-pwdveps", par_lv.dvol_eps,
      "dvol_eps: fim da ejecao quando dV > -dvol_eps, mL", G);
  for (ChamberParams * p : pp)
  {
    p->period = par_lv.period;
    p->t_zero = t_zero;
    p->t_end_diastole = t_zero + t_ed_dur;
    p->min_ejection_duration = min_ej;
    p->dvol_eps = dveps;
  }

  struct Flags { const char *pre, *ped, *pfill, *pinit, *C, *R, *Rc, *rate, *tag; };
  const Flags f[2] = {
    { "-pwpds",    "-pwedplv", "-pwpla", "-pwpej",    "-pwcart",    "-pwrart",    "-pwrc",    "-pwfillrate",    "VE" },
    { "-pwpds_rv", "-pwedprv", "-pwpra", "-pwpej_rv", "-pwcart_rv", "-pwrart_rv", "-pwrc_rv", "-pwfillrate_rv", "VD" } };
  for (int c = 0; c < 2; c++)
  {
    ChamberParams & p = *pp[c];
    const std::string t = std::string(" do ") + f[c].tag;
    p.preload_pressure = CommandLineArgs::read(f[c].pre, p.preload_pressure,
        (std::string("preload_pressure") + t + ", Pa").c_str(), G);
    p.p_end_diastole = CommandLineArgs::read(f[c].ped, p.p_end_diastole,
        (std::string("p_end_diastole") + t + ", Pa").c_str(), G);
    p.p_fill = CommandLineArgs::read(f[c].pfill, p.p_fill,
        (std::string("p_fill") + t + ": IVR -> FILLING quando P < p_fill, Pa").c_str(), G);
    p.windkessel.p_init = CommandLineArgs::read(f[c].pinit, p.windkessel.p_init,
        (std::string("Windkessel p_init") + t + " (P_c inicial), Pa").c_str(), G);
    p.windkessel.C = CommandLineArgs::read(f[c].C, p.windkessel.C,
        (std::string("Windkessel compliance") + t + ", mL/Pa").c_str(), G);
    p.windkessel.R = CommandLineArgs::read(f[c].R, p.windkessel.R,
        (std::string("Windkessel resistance") + t + ", Pa*ms/mL").c_str(), G);
    p.windkessel.Rc = CommandLineArgs::read(f[c].Rc, p.windkessel.Rc,
        (std::string("Windkessel characteristic_impedance") + t + ", Pa*ms/mL").c_str(), G);
    p.fill_rate = CommandLineArgs::read(f[c].rate, p.fill_rate,
        (std::string("fill rate (prescribed inflow)") + t + ", mL/ms").c_str(), G);
  }

  if (par_lv.t_end_diastole >= par_lv.period)
    std::cout << " *** ATENCAO: t_end_diastole (" << par_lv.t_end_diastole
              << " ms) >= periodo (" << par_lv.period << " ms)." << std::endl;
}

// ======================================================================
//  Cycle control
// ======================================================================

void PiecewisePVLoop::initialize(double t0, double V_lv, double P_lv, double V_rv, double P_rv)
{
  ChamberState * ss[2] = { &st_lv, &st_rv };
  const ChamberParams * pp[2] = { &par_lv, &par_rv };
  const double V[2] = { V_lv, V_rv }, P[2] = { P_lv, P_rv };
  for (int c = 0; c < 2; c++)
  {
    ChamberState & s = *ss[c];
    s.phase = PRELOAD;
    s.n_beats = 0;
    s.last_phase_change = t0;
    s.volume_n = V[c];
    s.pressure_n = P[c];
    s.dvol_n = 0.0;
    s.wdk_pressure_n = pp[c]->windkessel.p_init;
    s.outflow_n = 0.0;
    s.has_ejected = false;
    s.end_dia_vol = V[c];
    s.end_sys_vol = 0.0;
    s.solved_under = PRELOAD;
  }
}

double PiecewisePVLoop::preload_pressure(const ChamberState & st, const ChamberParams & p,
                                         double t)
{
  // First ramp 0 -> preload_pressure compares absolute t, so only the first
  // beat starts from 0; later beats ramp preload_pressure -> p_end_diastole.
  if (t <= p.t_zero)
    return p.preload_pressure * t / p.t_zero;
  if (p.t_end_diastole > p.t_zero)
    return (p.p_end_diastole - p.preload_pressure)
         * (t - st.n_beats * p.period - p.t_zero) / (p.t_end_diastole - p.t_zero)
         + p.preload_pressure;
  return st.pressure_n;
}

PiecewisePVLoop::Constraint PiecewisePVLoop::constraint(const ChamberState & st,
                                                        const ChamberParams & p,
                                                        double t, double dt)
{
  Constraint c = { Constraint::PRESSURE, 0.0, 0.0, 0.0, 0.0 };
  switch (st.phase)
  {
    case PRELOAD:
      c.kind = Constraint::PRESSURE;
      c.P = preload_pressure(st, p, t);
      break;
    case IVC:
      c.kind = Constraint::VOLUME;
      c.V = st.end_dia_vol;
      break;
    case EJECTION:
      c.kind = Constraint::AFFINE;
      p.windkessel.ejection_law(st.wdk_pressure_n, st.volume_n, dt, c.A, c.B);
      break;
    case IVR:
      c.kind = Constraint::VOLUME;
      c.V = st.end_sys_vol;
      break;
    case FILLING:
      c.kind = Constraint::VOLUME;
      c.V = st.volume_n + p.fill_rate * dt;
      break;
  }
  return c;
}

bool PiecewisePVLoop::advance_phase(ChamberState & st, const ChamberParams & p, double t,
                                    const char * tag)
{
  const Phase old = st.phase;
  const double since = t - st.last_phase_change;

  if (st.phase == PRELOAD
      && t >= st.n_beats * p.period + p.t_end_diastole && t >= p.t_zero)
  {
    st.phase = IVC; st.last_phase_change = t;
  }
  else if (st.phase == IVC)
  {
    if (st.pressure_n > st.wdk_pressure_n)
    {
      st.phase = EJECTION; st.last_phase_change = t;
      st.has_ejected = true;
    }
  }
  else if (st.phase == EJECTION && since > p.min_ejection_duration)
  {
    const bool outflow_stopped = st.dvol_n > -p.dvol_eps;
    const bool pressure_below_arterial = st.pressure_n < st.wdk_pressure_n;
    if (outflow_stopped || pressure_below_arterial)
    {
      st.phase = IVR; st.last_phase_change = t;
    }
  }
  else if (st.phase == IVR && st.pressure_n < p.p_fill)
  {
    st.phase = FILLING; st.last_phase_change = t;
  }
  else if (st.phase == FILLING)
  {
    if (p.period > 0.0 && t >= (st.n_beats + 1) * p.period + p.t_zero)
    {
      st.phase = PRELOAD; st.last_phase_change = t;
      st.n_beats += 1;
    }
  }

  const bool changed = (st.phase != old);
  if (changed)
    std::cout << "  [PW " << tag << "] " << phase_name(old) << " -> " << phase_name(st.phase)
              << " at t=" << t << " ms (V=" << st.volume_n << " mL, P=" << st.pressure_n
              << " Pa, P_c=" << st.wdk_pressure_n << " Pa)" << std::endl;
  return changed;
}

void PiecewisePVLoop::commit_one(ChamberState & s, const ChamberParams & p, double t, double dt,
                                 double V, double P, const char * tag, bool & changed)
{
  const bool ejecting = (s.phase == EJECTION);
  double P_c, Q;
  p.windkessel.advance(s.wdk_pressure_n, V, s.volume_n, dt, ejecting, s.has_ejected, P_c, Q);
  s.wdk_pressure_n = P_c;
  s.outflow_n = Q;
  s.dvol_n = V - s.volume_n;
  s.volume_n = V;
  s.pressure_n = P;
  s.solved_under = s.phase;
  if (s.phase == PRELOAD)
    s.end_dia_vol = V;

  if (advance_phase(s, p, t, tag))
  {
    changed = true;
    if (s.phase == IVC)
      s.end_dia_vol = V;
    else if (s.phase == IVR)
      s.end_sys_vol = V;
  }
}

bool PiecewisePVLoop::commit(double t, double dt, double V_lv, double P_lv,
                             double V_rv, double P_rv)
{
  bool changed = false;
  commit_one(st_lv, par_lv, t, dt, V_lv, P_lv, "LV", changed);
  commit_one(st_rv, par_rv, t, dt, V_rv, P_rv, "RV", changed);
  return changed;
}

// ======================================================================
//  Queries and output
// ======================================================================

int PiecewisePVLoop::beat_of(double t) const
{
  // A step ending exactly at k*period still belongs to beat k-1
  return static_cast<int>(std::floor((t - 1e-9) / period()));
}

double PiecewisePVLoop::ejection_fraction_LV() const
{
  if (st_lv.end_dia_vol <= 0.0) return 0.0;
  return 100.0 * stroke_volume_LV() / st_lv.end_dia_vol;
}

double PiecewisePVLoop::ejection_fraction_RV() const
{
  if (st_rv.end_dia_vol <= 0.0) return 0.0;
  return 100.0 * stroke_volume_RV() / st_rv.end_dia_vol;
}

std::string PiecewisePVLoop::describe() const
{
  std::ostringstream os;
  os << " PV loop por fases, quase-estatico\n";
  os << "   unidades: mL, Pa, ms\n";
  os << "   periodo = " << period() << " ms; t_zero = " << par_lv.t_zero
     << " ms; t_end_diastole = " << par_lv.t_end_diastole << " ms (inicio da ativacao)\n";
  os << "   min_ejection_duration = " << par_lv.min_ejection_duration
     << " ms; dvol_eps = " << par_lv.dvol_eps << " mL\n";

  const ChamberParams * pp[2] = { &par_lv, &par_rv };
  const char * tag[2] = { "VE", "VD" };
  for (int c = 0; c < 2; c++)
  {
    const ChamberParams & p = *pp[c];
    os << "   " << tag[c] << ": preload_pressure = " << p.preload_pressure
       << ", p_end_diastole = " << p.p_end_diastole << ", p_fill = " << p.p_fill << " Pa\n"
       << "       Windkessel: p_init = " << p.windkessel.p_init << " Pa, C = "
       << p.windkessel.C << " mL/Pa, R = " << p.windkessel.R << " Pa*ms/mL, R_c = "
       << p.windkessel.Rc << " Pa*ms/mL\n"
       << "       enchimento: vazao prescrita " << p.fill_rate << " mL/ms\n";
  }
  return os.str();
}

std::string PiecewisePVLoop::history_header()
{
  return "# t[ms] phase_LV P_LV[Pa] V_LV[mL] Q_LV[mL/ms] Pc_LV[Pa]"
         " phase_RV P_RV[Pa] V_RV[mL] Q_RV[mL/ms] Pc_RV[Pa]"
         "   (fases: 0 PRELOAD 1 IVC 2 EJECTION 3 IVR 4 FILLING;"
         " a fase e a do passo resolvido)";
}

std::string PiecewisePVLoop::history_row(double t) const
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(10) << t;
  for (const ChamberState * s : { &st_lv, &st_rv })
    os << " " << static_cast<int>(s->solved_under) << " " << s->pressure_n << " "
       << s->volume_n << " " << s->outflow_n << " " << s->wdk_pressure_n;
  return os.str();
}
