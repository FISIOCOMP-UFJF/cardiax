#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include <fstream>
#include <cmath>
#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal_reactive.hpp"

class Electromechanic
{
public:
  Electromechanic(const std::string &epmodel) {};
  ~Electromechanic();
  CardiacProblem &ref() { return *ephy; }  
  void config(const string &basename);
  void solve();
  void set_checkpoint_interval(double interval) {checkpoint_interval = interval;};

private:
  
  CardiacProblem* ephy; // Can be eikonal, Monodomain or MonodomainDeformation
  UpdatedLagrangian elas;
  
  std::string filename;
  std::string ep_type;
  std::string cell_model; 
  std::vector<arma::mat33*> vec_stress;
  std::vector<arma::mat33*> vec_fib;
  std::vector<arma::mat33*> vec_fib0;
  arma::vec ta;

  TimerSection timer;
  
  int checkpoint_interval; 
  int steps_per_mech = 1;   // EP steps between mechanics solves
  int steps_per_out  = 1;   // EP steps between output frames
  double t_final;
  double dt_mech;
  double dt_ephy;
  double dt_out = 1.0;
};

#endif
