#ifndef IMPLICIT_ODE_SOLVER_H
#define IMPLICIT_ODE_SOLVER_H

#include "ode_solver.hpp"
#include "cellmodel.hpp"

class ImplicitODESolver : public ODESolver
{
public:
  
  // Constructor
  ImplicitODESolver(CellModel * ode);

  // Destructor
  virtual ~ImplicitODESolver() {};

  // Get the absolute tolerance of the Newton method
  double get_absolute_tolerance() const { return atol; }

  // Get the relative tolerance of the Newton method
  double get_relative_tolerance() const { return rtol; }

  // Change absolute tolerance
  void set_abolute_tolerance(double tol) { atol = tol; }

  // Change relative tolerance
  void set_relative_tolerance(double tol) { rtol = tol; }

  // Function to advance one time step
  virtual void advance(double *y, double & t, double & dt) = 0;

protected:

  // Some integers for Newton method
  int maxits;

  // Size of the Jacobian matrix
  int jsize;

  // Do we need to update the jacobian matrix?
  bool update_jacobian;

  // Absolute tolerance for the Newton method
  double atol;

  // Relative tolerance for the Newton method
  double rtol;

  // Jacobian matrix
  arma::mat jac;

  // Solution of the linear system (increment to x)
  arma::vec dy;

  // Residual vector
  arma::vec r;
 
};
 
#endif
