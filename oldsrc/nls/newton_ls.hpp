#ifndef NEWTON_LS_HPP
#define	NEWTON_LS_HPP

#include "nonlinear_problem.hpp"
#include "newton.hpp"

class NewtonLineSearch : public Newton
{
 public:
  NewtonLineSearch(NonlinearProblem * p);
  ~NewtonLineSearch();
  int solve();
  void line_search(double & eta0, double & eta, double & rtu0, double & rtu);
  

 private:
   int num_search;
   int max_search;
   double search;
};

#endif	/* NEWTON_LS_HPP */

