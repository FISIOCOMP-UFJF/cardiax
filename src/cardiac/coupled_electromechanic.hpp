#ifndef MONOMECHANIC_H
#define MONOMECHANIC_H

#include "bidomain_deformation.hpp"
#include "monodomain_deformation.hpp"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"
#include "eikonal.hpp"
#include <fstream>

class Electromechanic
{
public:
  Electromechanic(const std::string &epmodel);

  ~Electromechanic();

  void config(const string &basename);

  CardiacProblem &ref() { return *ephy; }
  
  void solve();
 
private:
  double dt_mech;
  string filename;
  CardiacProblem* ephy; // Can be eikonal, monodomain or monodomainderformation
  std::string ep_type;
  UpdatedLagrangian elas;
  std::vector<arma::mat33 *> vec_stress;
  std::vector<arma::mat33 *> vec_fib;
  std::vector<arma::mat33 *> vec_fib0;
  TimerSection timer;
  string cell_model; 

  arma::vec ta;

};

#endif
