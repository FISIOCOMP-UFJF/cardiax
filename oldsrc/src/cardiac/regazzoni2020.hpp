#ifndef REGAZZONI2020_HPP
#define REGAZZONI2020_HPP

#include <array>
#include <string>
#include <vector>
#include <functional>
#include <cmath>

/*!
 *  Closed-loop lumped-parameter (0D) circulation model.
 *
 *  Regazzoni et al. (2022), "A cardiac electromechanical model coupled with a
 *  lumped-parameter model for closed-loop blood circulation".
 *
 *  Four chambers (LA, LV, RA, RV), four valves (MV, AV, TV, PV) and two
 *  circuits (systemic and pulmonary), each split into arterial and venous
 *  RLC compartments.
 *
 *  ---------------------------------------------------------------------
 *  UNITS: mL, mmHg, s  (consistent with the original model)
 *      resistance  mmHg*s/mL
 *      compliance  mL/mmHg
 *      inertance   mmHg*s^2/mL
 *  The coupling with the 3D mechanics converts pressures at the interface.
 *  ---------------------------------------------------------------------
 *
 *  3D-0D coupling
 *  --------------
 *  By default every chamber uses a time-varying elastance
 *      p(t) = E(t) (V - V0).
 *  When a coupling callback is registered with set_BiV_coupling(), the LV and
 *  RV pressures come instead from the 3D finite element model: the 0D model
 *  hands it the current cavity volumes and receives back the pressures that
 *  reproduce them. The remaining chambers stay 0D.
 */
class Regazzoni2020
{
public:
  //! Number of ODE state variables
  static const int NSTATE = 12;

  //! State indices
  enum StateIdx {
    V_LA = 0, V_LV, V_RA, V_RV,
    p_AR_SYS, p_VEN_SYS, p_AR_PUL, p_VEN_PUL,
    Q_AR_SYS, Q_VEN_SYS, Q_AR_PUL, Q_VEN_PUL
  };

  //! Time-varying elastance parameters of a 0D chamber
  struct Chamber
  {
    double EA;   //!< active (max) elastance   [mmHg/mL]
    double EB;   //!< passive (min) elastance  [mmHg/mL]
    double TC;   //!< contraction duration     [s]
    double TR;   //!< relaxation duration      [s]
    double tC;   //!< contraction onset        [s]
    double V0;   //!< unloaded volume          [mL]
  };

  //! Diode-like valve: low resistance when open, high when closed
  struct Valve
  {
    double Rmin; //!< [mmHg*s/mL]
    double Rmax; //!< [mmHg*s/mL]
  };

  //! Arterial or venous RLC compartment
  struct Compartment
  {
    double R_AR, C_AR, L_AR;
    double R_VEN, C_VEN, L_VEN;
  };

  struct Parameters
  {
    double HR;                    //!< heart rate [Hz]
    Chamber LA, LV, RA, RV;
    Valve MV, AV, TV, PV;
    Compartment SYS, PUL;
  };

  //! Signature of the 3D coupling: given target cavity volumes and time,
  //! return the LV and RV pressures in mmHg.
  typedef std::function<void(double V_lv, double V_rv, double t,
                             double & p_lv, double & p_rv)> BiVCoupling;

  Regazzoni2020();

  //! Default parameter set from Regazzoni et al.
  static Parameters default_parameters();

  //! Default (already near periodic) initial state
  static std::array<double,NSTATE> default_initial_state();

  void set_parameters(const Parameters & p) { par = p; }
  Parameters & parameters() { return par; }

  void set_state(const std::array<double,NSTATE> & s) { y = s; }
  const std::array<double,NSTATE> & state() const { return y; }
  double & operator[](int i) { return y[i]; }

  //! Register the 3D mechanics coupling for LV and RV.
  void set_BiV_coupling(BiVCoupling f) { coupling = f; has_coupling = true; }
  void clear_BiV_coupling() { has_coupling = false; }

  //! Override only the initial LV/RV volumes (e.g. from the 3D mesh)
  void set_initial_volumes(double v_lv, double v_rv)
  { y[V_LV] = v_lv; y[V_RV] = v_rv; }

  //! Right-hand side of the ODE system. Fills dy and refreshes the derived
  //! quantities (pressures and flows) accessible through the getters below.
  void rhs(double t, const std::array<double,NSTATE> & yy,
           std::array<double,NSTATE> & dy);

  //! Advance one step. Forward Euler mirrors the reference implementation;
  //! it is conditionally stable, so dt must be small (1e-3 s is typical).
  void step(double t, double dt);

  //! Chamber period
  double period() const { return 1.0 / par.HR; }

  // ---- derived quantities, valid after the last rhs() evaluation ----
  double pressure_LA() const { return p_la; }
  double pressure_LV() const { return p_lv; }
  double pressure_RA() const { return p_ra; }
  double pressure_RV() const { return p_rv; }
  double flow_MV() const { return q_mv; }
  double flow_AV() const { return q_av; }
  double flow_TV() const { return q_tv; }
  double flow_PV() const { return q_pv; }

  //! Total blood volume currently held in the model
  double total_volume() const;

  //! Header and one row for a csv-style history file
  static std::string history_header();
  std::string history_row(double t) const;

private:
  Parameters par;
  std::array<double,NSTATE> y;

  BiVCoupling coupling;
  bool has_coupling;

  // derived (static) variables refreshed by rhs()
  double p_la, p_lv, p_ra, p_rv;
  double q_mv, q_av, q_tv, q_pv;

  //! Blanco time-varying elastance, normalised activation in [0,1]
  double activation(const Chamber & c, double t) const;

  //! p = E(t) (V - V0)
  double chamber_pressure(const Chamber & c, double V, double t) const;

  //! Smooth diode resistance; ~Rmin when p1 > p2 (open), ~Rmax otherwise
  double valve_resistance(const Valve & v, double p1, double p2) const;

  double flux_through_valve(const Valve & v, double p1, double p2) const
  { return (p1 - p2) / valve_resistance(v, p1, p2); }
};

#endif
