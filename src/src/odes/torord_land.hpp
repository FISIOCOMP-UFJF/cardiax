#ifndef TORORD_LAND_H_
#define TORORD_LAND_H_

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/*!
 *  Torord Land
*/
class TorordLand : public CellModel
{
 public:
  
  //! Default constructor
  TorordLand();

  //! Set initial conditions for model
  virtual void init(double * values) const;

  //! Compute RHS equations for the model
  virtual void equation(const double t, const double * sv, double * values);

private:
  double active;

};

#endif
