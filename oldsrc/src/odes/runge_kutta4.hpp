#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

#include "ode_solver.hpp"
#include "cellmodel.hpp"

/** Concrete ODESolver 
    Derived class that implements the 4th order Runge-Kutta method
*/
class RungeKutta4 : public ODESolver
{
 public:

  // Constructor
  RungeKutta4 (CellModel *p);
  
  // Interface
  virtual void advance(double *y, double & t, double & dt);  

};

#endif
