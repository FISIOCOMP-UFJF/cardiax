#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"

/*!
 * Class for coupled electromechanics problem:
 *  - Electrics: Monodomain or Bidomain models
 *  - Mechanics: TotalLagrangian or UpdatedLagrangian formulations
 *               for non-linear elasticity
 */

class Electromechanic
{

public:

  //! Constructor
  Electromechanic(const std::string & epmodel);

  //! Destructor
  ~Electromechanic();

  //! Configure both problems
  void config(const string & basename);

  //! The CardiacProblem solver
  MonodomainDeformation & ref() { return ephy; }
  //BidomainDeformation & ref() { return ephy; }

  //! Solve monodomain and elasticity problems
  void solve();

 private:

  //! Time step for mechanical problem (miliseconds)
  double dt_mech;

  int n_cycles = 1;

  double Ta0 = 0.0;
  double P0;

  double dTa0;

  double C_art;
  double R_per;
  double P_o;
  double part;
  double stroke_volume;

  string filename;

  //! Cardiac problem
  //CardiacProblem * ep;
  MonodomainDeformation ephy;
  //BidomainDeformation ephy;

  //! Mechanical problem
  //TotalLagrangian elas;
  UpdatedLagrangian elas; //! TODO: Updated Lagrangian is not working!

  //! Vector of stress tensors (used to load cardiac EP problem)
  std::vector<arma::mat33*> vec_stress;

  //! Vector with the current fibers directions (f,s,n)
  std::vector<arma::mat33*> vec_fib;

  //! Vector with the reference fibers directions (f0,s0,n0)
  std::vector<arma::mat33*> vec_fib0;

  //! Timer for sections
  TimerSection timer;

  //std::vector< std::pair<double, double> > p_Ta_list;
  std::vector<double> curr_time;
  std::vector<double> recorded_time;
  std::vector<double> p_lv;
  std::vector<double> Ta_list;

  std::vector<double> volume;
  std::vector<double> p_art;
  std::vector<double> p_ven;
  std::vector<double> p_LA;

  void Solve_System(double active_stress, double pressure);

  double Isovolumetric_PressureUpdate(int i, double volume_constraint);

  double PressureUpdatePhase3(int i, double V_ED);

  double fullLumpedModel(int i, double volume_constraint);

};

#endif

