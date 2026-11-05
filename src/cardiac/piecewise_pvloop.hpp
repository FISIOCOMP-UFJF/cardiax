#ifndef PIECEWISE_PVLOOP_HPP
#define PIECEWISE_PVLOOP_HPP

#include <string>

/*!
 *  Five-phase cardiac cycle controller (Wang et al. 2025, eLife 14:RP106555,
 *  "Mechanical model").
 *
 *  Each cavity has one pressure unknown. At each step the current phase sets
 *  one constraint on it, built from the state at the START of the step:
 *
 *    0 PRELOAD    pressure: linear ramp (0 -> preload_pressure up to t_zero in
 *                 the first beat; then preload_pressure -> p_end_diastole
 *                 between t_zero and t_end_diastole of each beat)
 *    1 IVC        volume:  V = V_ED
 *    2 EJECTION   pressure affine in V: 3-element Windkessel, implicit Euler,
 *                     Q   = -(V - V_n)/dt
 *                     P_c = (P_c,n + dt Q/C) / (1 + dt/(R C))
 *                     P   = P_c + R_c Q          ->  P = A + B V
 *    3 IVR        volume:  V = V_ES
 *    4 FILLING    volume:  V = V_n + rate*dt (prescribed inflow)
 *
 *  Transitions, evaluated on the converged state at the END of the step:
 *    PRELOAD -> IVC      t >= n_beats*period + t_end_diastole
 *    IVC -> EJECTION     P > P_c (arterial valve opens)
 *    EJECTION -> IVR     after min_ejection_duration since opening:
 *                        dV > -dvol_eps (outflow stopped) or P < P_c
 *    IVR -> FILLING      P < p_fill
 *    FILLING -> PRELOAD  t >= (n_beats+1)*period + t_zero
 *  With the valve closed, P_c is held until the cavity first ejects and then
 *  drains through the peripheral resistance (implicit Euler).
 *
 *  This class is only the phase policy; it knows nothing of the mesh or the
 *  solver. The caller imposes the constraint (see PiecewiseCirculation) and
 *  calls constraint() and commit().
 *
 *  Units: time ms | pressure Pa | volume mL | resistance Pa*ms/mL |
 *         compliance mL/Pa | flow mL/ms.
 */
class PiecewisePVLoop
{
public:

  enum Phase
  {
    PRELOAD  = 0,
    IVC      = 1,   //!< isovolumic contraction
    EJECTION = 2,
    IVR      = 3,   //!< isovolumic relaxation
    FILLING  = 4
  };

  static const char * phase_name(Phase p);

  //! 3-element Windkessel. Stateless: P_c lives in ChamberState.
  struct Windkessel
  {
    double p_init;   //!< initial P_c [Pa]
    double C;        //!< compliance [mL/Pa]
    double R;        //!< peripheral resistance [Pa*ms/mL]
    double Rc;       //!< characteristic impedance [Pa*ms/mL] (0 = 2-element)

    //! Coefficients of P = A + B V imposed during ejection
    void ejection_law(double P_c, double V_n, double dt, double & A, double & B) const;

    //! Advances P_c over a converged step; returns (new P_c, Q)
    void advance(double P_c, double V, double V_n, double dt, bool ejecting,
                 bool has_ejected, double & P_c_new, double & Q) const;
  };

  //! Cycle parameters of one cavity
  struct ChamberParams
  {
    double t_zero;                 //!< [ms]
    double preload_pressure;       //!< [Pa]
    double t_end_diastole;         //!< [ms], from the start of the beat
    double p_end_diastole;         //!< [Pa]
    double p_fill;                 //!< [Pa]
    double period;                 //!< [ms]
    Windkessel windkessel;
    double fill_rate;              //!< prescribed inflow [mL/ms]
    double min_ejection_duration;  //!< [ms]
    double dvol_eps;               //!< [mL]
  };

  //! Cycle state of one cavity
  struct ChamberState
  {
    Phase  phase;
    int    n_beats;
    double last_phase_change;   //!< [ms]
    double volume_n;            //!< [mL]
    double pressure_n;          //!< [Pa]
    double dvol_n;              //!< [mL]
    double wdk_pressure_n;      //!< P_c [Pa]
    double outflow_n;           //!< Q [mL/ms]
    bool   has_ejected;
    double end_dia_vol;         //!< [mL]
    double end_sys_vol;         //!< [mL]
    Phase  solved_under;        //!< phase the last step was solved under (history)
  };

  //! Constraint of one step on the pressure of a cavity
  struct Constraint
  {
    enum Kind { PRESSURE, VOLUME, AFFINE } kind;
    double P;      //!< PRESSURE: p = P
    double V;      //!< VOLUME:   V(p) = V
    double A, B;   //!< AFFINE:   p = A + B V(p)   [Pa, Pa/mL]
  };

  PiecewisePVLoop();

  // ---- configuration -------------------------------------------------

  //! <pvloop_piecewise> block of the mesh (if present), then the -pw* flags.
  //! rr_ms > 0 (global -rr flag) overrides the period.
  void configure(const std::string & mshfile, double rr_ms = -1.0);

  ChamberParams & lv() { return par_lv; }
  ChamberParams & rv() { return par_rv; }
  const ChamberParams & lv() const { return par_lv; }
  const ChamberParams & rv() const { return par_rv; }

  // ---- life cycle -----------------------------------------------------

  //! Current volume and pressure; every cavity starts in PRELOAD
  void initialize(double t0, double V_lv, double P_lv, double V_rv, double P_rv);

  //! Constraints of the step ending at t, from the state at its start
  Constraint constraint_LV(double t, double dt) const { return constraint(st_lv, par_lv, t, dt); }
  Constraint constraint_RV(double t, double dt) const { return constraint(st_rv, par_rv, t, dt); }

  //! Converged step at t: advances the Windkessel and picks the next phase.
  //! Returns true if any cavity changed phase.
  bool commit(double t, double dt, double V_lv, double P_lv, double V_rv, double P_rv);

  // ---- queries --------------------------------------------------------

  const ChamberState & state_LV() const { return st_lv; }
  const ChamberState & state_RV() const { return st_rv; }
  void set_states(const ChamberState & lv, const ChamberState & rv) { st_lv = lv; st_rv = rv; }

  double period() const { return par_lv.period; }

  //! Activation onset of the first beat [ms]: end of diastole
  double activation_time() const { return par_lv.t_end_diastole; }

  //! Run length: num_beats * period [ms]
  double total_time(int num_beats) const { return num_beats * period(); }

  //! Beat of the step ending at t (0-based)
  int beat_of(double t) const;

  //! SV/EF from the latest V_ED/V_ES
  double stroke_volume_LV() const { return st_lv.end_dia_vol - st_lv.end_sys_vol; }
  double stroke_volume_RV() const { return st_rv.end_dia_vol - st_rv.end_sys_vol; }
  double ejection_fraction_LV() const;
  double ejection_fraction_RV() const;

  std::string describe() const;

  //! History: one row per step; the phase written is the one the step was
  //! SOLVED under.
  static std::string history_header();
  std::string history_row(double t) const;

private:
  ChamberParams par_lv, par_rv;
  ChamberState  st_lv,  st_rv;

  static Constraint constraint(const ChamberState & st, const ChamberParams & p,
                               double t, double dt);
  static double preload_pressure(const ChamberState & st, const ChamberParams & p, double t);
  static bool advance_phase(ChamberState & st, const ChamberParams & p, double t,
                            const char * tag);
  static void commit_one(ChamberState & st, const ChamberParams & p, double t, double dt,
                         double V, double P, const char * tag, bool & changed);
};

#endif
