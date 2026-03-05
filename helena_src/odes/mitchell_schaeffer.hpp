#ifndef MITCHELLSCHAEFFER_H
#define MITCHELLSCHAEFFER_H

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/** A Two-Current Model for the Dynamics of Cardiac Membrane
    Colleen C. Mitchell and David G. Schaeffer
*/
class MitchellSchaeffer : public CellModel
{
 public:
  
  /// Default constructor
  MitchellSchaeffer();
  
  /// Set initial conditions for MS model
  virtual void init(double * values) const;
  
  /// Compute RHS equations for the MS model
  virtual void equation(const double time, const double * sv, double * values);

private:

  double rateswitchGK(double v, double h);
  double rateswitchGK2(double v);

};

#endif
