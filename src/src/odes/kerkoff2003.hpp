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
