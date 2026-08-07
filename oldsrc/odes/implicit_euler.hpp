#ifndef IMPLICIT_EULER_H
#define IMPLICIT_EULER_H

#include "implicit_ode_solver.hpp"
#include "cellmodel.hpp"

class ImplicitEuler : public ImplicitODESolver
{
public:
  
  /// Constructor
  ImplicitEuler(CellModel * ode);

  /// Destructor
  ~ImplicitEuler() {};

  /// Function to advance one time step
  void advance(double *y, double & t, double & dt);

private:

  /// Matrix for the linear system
  arma::mat A;
  
  /// Identity matrix
  arma::mat I;

  /// Right-hand side of the linear system
  arma::vec b;

  /// Old y values (=y_n)
  arma::vec yo;

  /// Current Newton iterate
  arma::vec yn;
  arma::vec yc;

  /// Auxiliary vector for computing f
  arma::vec f;

};
 
#endif
