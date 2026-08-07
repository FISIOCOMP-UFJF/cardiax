#ifndef KERKOFF2003_H
#define KERKOFF2003_H

#include<cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/** A Two-Current Model for the Dynamics of Cardiac Membrane
    Colleen C. Mitchell and David G. Schaeffer
*/
class Kerkoff2003 : public CellModel
{
 public:
  

  //! The Kerckhoffs equations are written in SECONDS (tr = 0.075 s,
  //! tmax = 0.483 s), unlike the ionic models which use ms.
  double native_time_unit_ms() const override { return 1000.0; }

  //! State variable 1 holds the activation time
  int lat_var_index() const override { return 1; }
  /// Default constructor
  Kerkoff2003();
  
  /// Set initial conditions for MS model
  virtual void init(double * values) const;
  
  /// Compute RHS equations for the MS model
  virtual void equation(const double time, const double * sv, double * values);
 
  private:
    double active_stress; 
};

#endif