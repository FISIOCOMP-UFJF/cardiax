#ifndef SIMPLEODE_H
#define SIMPLEODE_H

#include "cellmodel.hpp"
#include "ode_solver.hpp"

/** A simple ODE with exact solution used to test ODESolvers

    Problem:
      du/dt = - 2 u
      u[0] = 5

    Solution
      u(t) = 5 e^{-at}

*/

class SimpleODE : public CellModel
{
 public:

  /// Default constructor
  SimpleODE();
  
  /// Set initial conditions for SimpleODE
  virtual void init(double * values) const;

  /// Compute RHS equations for SimpleODE
  virtual void equation(const double t, const double * sv, double * values);

  /// New member function to compute exact solution and save to file
  void dump_exact();

};

#endif
