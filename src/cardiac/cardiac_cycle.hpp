#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal.hpp"
#include "regazzoni2020.hpp"
#include <fstream>

class CardiacElectromechanic
{
public:
  CardiacElectromechanic(const std::string &epmodel);

  ~CardiacElectromechanic();

  void config(const string &basename);

  // MonodomainDeformation &ref() { return ephy; }
  Eikonal &ref() { return ephy; }
  
  void solve();

  //! Solve with the closed-loop Regazzoni 0D circulation driving both
  //! cavities. The 0D model supplies target volumes, the 3D model returns
  //! the pressures that reproduce them.
  void solve_circulation(int num_beats = 1, double dt_circ = 1.0e-3);

  //! Pressure conversion: model pressure units per mmHg.
  //! Default 133.322 (mesh/material in Pa). Use 0.133322 for kPa.
  void set_pressure_unit_per_mmHg(double f) { mmHg_to_p = f; }

  //! Electrophysiology time units per second: 1 if the EP problem is set up
  //! in seconds, 1000 if in milliseconds. Left at 0 it is auto-detected from
  //! the EP final time.
  void set_ephy_time_per_second(double f) { ephy_time_per_second = f; }

  //! Where the active tension comes from.
  //! var_index < 0  : monitored value 0 (phenomenological models such as
  //!                  Kerckhoffs, which expose active_stress as monitored)
  //! var_index >= 0 : state variable var_index (ionic models such as
  //!                  ToRORdLand, whose Ta is state 49)
  //! scale multiplies the raw value.
  void set_active_tension_source(int var_index, double scale)
  { ta_var_index = var_index; ta_scale = scale; }

  //! State variable holding the membrane potential, for the per-node output.
  //! Ionic models (ToRORdLand) keep Vm in state 0; phenomenological models
  //! have no membrane potential at all and must pass -1, which disables the
  //! vm output instead of writing some unrelated state variable into it.
  void set_vm_source(int var_index) { vm_var_index = var_index; }

  //! Stimulus applied at each node once its LAT is reached.
  //! Amplitude in the cell-model current unit, duration in the solver time
  //! unit (2 ms = 0.002 when the solver runs in seconds). Amplitude 0
  //! disables it, which is correct for LAT-driven models like Kerckhoffs.
  void set_lat_stimulus(double amplitude, double duration)
  { stim_amplitude = amplitude; stim_duration = duration; }

public:
  void saveActiveStressToFile(const std::string &filename)
  {
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
      std::cerr << "Error: could not open file " << filename << std::endl;
    }
    for (size_t i = 0; i < activeStressCurve.size(); ++i)
    {
      outFile << timePoints[i] << " " << activeStressCurve[i] << "\n";
    }
    outFile.close();
    std::cout << "Active tension data saved to: " << filename << std::endl;
  }

private:
  double T_ref;
  double dt_mech;
  int n_cycles = 1;
  double Ta0 = 0.0;
  double P0;
  double dTa0;
  double c_art, c_ven, Rmv, Rven, Rao, Rper;
  double V_art_zero, V_ven_zero;
  double E_es_LA, A_LA, B_LA, Tmax, tau;
  double P_o, part, pven, stroke_volume;
  string filename;
  // MonodomainDeformation ephy;
  Eikonal ephy; 
  UpdatedLagrangian elas;
  std::vector<arma::mat33 *> vec_stress;
  std::vector<arma::mat33 *> vec_fib;
  std::vector<arma::mat33 *> vec_fib0;
  TimerSection timer;
  std::vector<double> curr_time;
  std::vector<double> recorded_time;
  std::vector<double> p_lv, p_rv, Ta_list, volume, p_art, p_ven, p_LA;
  std::vector<double> volume_rv;   //!< RV cavity volume history
  
  arma::vec ta;

  int vm_var_index = -1;           //!< state variable holding Vm (<0 = none)
  arma::vec vm_node;               //!< membrane potential per node [mV]

  //! Write the per-node membrane potential and active tension for the given
  //! output frame, at the same rate as the displacement field.
  void save_node_fields(int frame);

  void Solve_System(double tt, double pressure, double pressure2);

  // ---- closed-loop circulation coupling ----------------------------
  Regazzoni2020 circ;              //!< 0D closed-loop model
  double mmHg_to_p = 133.322387415;//!< model pressure unit per mmHg (Pa)
  double kPa_to_p = 1000.0;        //!< model pressure units per kPa (Pa)
  double ephy_time_per_second = 0.0;//!< EP time units per second (0 = detect)
  int    ta_var_index = -1;        //!< <0 = monitored 0; else state variable
  double ta_scale = 0.0;           //!< 0 = use T_ref
  double stim_amplitude = 0.0;     //!< LAT stimulus amplitude (0 = disabled)
  double stim_duration = 0.002;    //!< LAT stimulus duration (solver units)

  bool   quiet_solve = false;      //!< suppress the per-solve pressure line
  int    last_newton_its = 0;      //!< pressure-inversion stats, per step
  double last_dV = 0.0;
  int    last_solves = 0;
  int    circ_log_every = 25;      //!< detailed log every N circulation steps

  //! Read the active tension from the cell model into `ta`
  void update_active_tension();

  //! Fill `ka` with the nodal active stiffness of the Land submodel, used by
  //! the mechanics to stabilise the velocity feedback. Clears `ka` (which
  //! disables the stabilisation) whenever it does not apply.
  void update_active_stiffness();

  arma::vec ka;                    //!< active stiffness per node, model units

  //! State-variable indices of the Land distortion submodel, used only by the
  //! negative-Ta diagnostic. Set for ToRORdLand; left at -1 for any model
  //! whose Ta does not decompose this way, which disables the diagnostic.
  //!     Ta = Lfac*(Tref/dr)*[ (ZETAS+1)*XS + ZETAW*XW ]
  void set_land_state_indices(int xs, int xw, int zetas, int zetaw)
  { land_xs_index = xs; land_xw_index = xw;
    land_zetas_index = zetas; land_zetaw_index = zetaw; }

  int land_xs_index    = -1;
  int land_xw_index    = -1;
  int land_zetas_index = -1;
  int land_zetaw_index = -1;

  //! Scratch buffers for the diagnostic, kept as members so the per-step
  //! report does not reallocate four nodal vectors every call.
  arma::vec land_xs, land_xw, land_zetas, land_zetaw;

  //! Report WHICH term of the Land active tension went negative, and at what
  //! stretch rate. Called from update_active_tension() when ta.min() < 0.
  void report_negative_ta();

  // ---- mechano-electric feedback: fibre stretch --------------------
  // After each mechanical solve, lambda_f = sqrt(I4f) is evaluated per
  // element, averaged onto the nodes, and handed to the cell model, where it
  // drives h(lambda) and the length-dependent ca50 of the Land submodel.
  bool   lambda_coupling = false;  //!< -lam 1 enables it
  bool   lambda_rate_on  = false;  //!< -lamrate 1 also feeds d(lambda)/dt
  double lambda_min_clip = 0.7;    //!< guard rails on the value handed over
  double lambda_max_clip = 1.3;

  //! Guard rail on d(lambda)/dt, in the CELL MODEL's own time unit (1/ms for
  //! ToRORd). The Land distortion equation has
  //!     ZETAS_ss = (A/cds) * lambda_rate = 249 * lambda_rate,
  //! so lambda_rate = -0.004 /ms (4 lengths/s, the model's Vmax) already
  //! drives ZETAS to -1, where (ZETAS+1) changes sign and the active tension
  //! goes NEGATIVE. A negative Ta is a compressive fibre stress and inverts
  //! elements. The default 0.003 keeps ZETAS above -0.75 with the stock
  //! parameters. Set to 0 (-lamratemax 0) to disable the clamp entirely and
  //! recover the previous, unguarded behaviour.
  double lambda_rate_max_clip = 0.003;

  //! Number of INITIAL calls to update_fiber_stretch() during which the rate
  //! is forced to zero (-lamratedelay N, default 0 = previous behaviour).
  //! The passive inflation takes lambda from 1 to ~1.1 over the first few
  //! steps; differencing that transient produces d(lambda)/dt of order
  //! 0.01 /ms, several times the threshold at which the active tension
  //! changes sign, before the mechanics has settled into anything physical.
  //! lambda itself (and therefore ca50 and h(lambda)) stays coupled
  //! throughout -- only the distortion terms ZETAS/ZETAW are held off.
  //! -lamstab 1 (default) enables the Regazzoni-Quarteroni stabilisation of
  //! the velocity-dependent active tension. Set to 0 to recover the plain
  //! segregated scheme, which is unstable whenever the active stiffness
  //! exceeds the passive one -- useful only for reproducing older runs.
  bool stabilize_active = true;

  int  lambda_rate_delay = 0;
  int  lambda_step_count = 0;      //!< calls to update_fiber_stretch() so far

  arma::vec lambda_node;           //!< lambda_f per node, after clipping
  arma::vec lambda_raw;            //!< lambda_f per node, BEFORE clipping
  arma::vec lambda_prev;           //!< previous lambda_raw, for the rate
  arma::vec lambda_rate_node;      //!< d(lambda_f)/dt, cell-model time units
  bool   lambda_prev_valid = false;

  //! Evaluate lambda_f from the current mechanical state and push it into
  //! the cell model. dt_solver is the interval since the previous call, in
  //! the EP SOLVER's time unit (s or ms, per -tunit); it is converted to the
  //! cell model's own unit internally and is only used for the rate.
  void update_fiber_stretch(double dt_solver);
  double p_lv_cur = 0.0;           //!< current LV pressure [model units]
  double p_rv_cur = 0.0;           //!< current RV pressure [model units]
  arma::mat22 Jpv;                 //!< dV/dp Jacobian [mL / model pressure]
  bool Jpv_valid = false;
  int circ_solves = 0;             //!< 3D solves performed (cost counter)
  int circ_out_every = 8;          //!< save the mesh every N circulation steps
  int circ_n_frames = 1;           //!< frames allocated for the writer
  std::ofstream circ_file;

  //! Invert the 3D model: find the cavity pressures reproducing the target
  //! volumes. Returns pressures in mmHg through plv/prv.
  void pressures_for_volumes(double Vlv_target, double Vrv_target, double t,
                             double & plv_mmHg, double & prv_mmHg);

  //! Rebuild the 2x2 dV/dp Jacobian by finite differences (2 extra solves).
  void update_pv_jacobian(double t);


  std::vector<double> activeStressCurve;
  std::vector<double> timePoints;
  bool curveCalculated = false;
  bool has_eikonal = false; 
};

#endif