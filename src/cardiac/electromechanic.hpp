#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal_reactive.hpp"
#include <fstream>

class Electromechanic
{
public:
  Electromechanic(const std::string &epmodel);
  ~Electromechanic();
  void config(const string &basename);
  CardiacProblem &ref() { return *ephy; }  
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
  
  TimerSection timer;
  
  arma::vec ta;
  
  int checkpoint_interval; 
  double dt_mech;
};

#endif
