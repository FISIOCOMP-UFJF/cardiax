#ifndef FITZHUGHNAGUMO_H
#define FITZHUGHNAGUMO_H

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/** FitzHugh-Nagumo model.
    For more details check the CellMl website.
*/
class FitzHughNagumo : public CellModel
{
 public:
  
  /// Default constructor
  FitzHughNagumo(); 
  
  /// Set initial conditions for FHN model
  virtual void init(double * values) const;
  
  /// Compute RHS equations for the FHN model
  virtual void equation(const double time, const double * sv, double * values);

};

#endif
