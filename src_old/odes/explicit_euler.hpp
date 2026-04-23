#ifndef FORWARDEULER_H
#define FORWARDEULER_H

#include <armadillo>
#include "ode_solver.hpp"
#include "cellmodel.hpp"

/** 
    Concrete ODESolver derived class that implements forward Euler method
*/
class ExplicitEuler : public ODESolver
{
public:

  // Constructor
  ExplicitEuler (CellModel *model);
 
  // Interface
  virtual void advance(double *y, double & t, double & dt);

private:

  // Values of derivatives for computations
  arma::vec dydt;

};

#endif
