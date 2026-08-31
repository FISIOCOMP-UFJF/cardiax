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
      std::cerr << "Erro: Nao foi possivel abrir o arquivo" << filename << std::endl;
    }
    for (size_t i = 0; i < activeStressCurve.size(); ++i)
    {
      outFile << timePoints[i] << " " << activeStressCurve[i] << "\n";
    }
    outFile.close();
    std::cout << "Dados da tensão ativa salvos em: " << filename << std::endl;
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