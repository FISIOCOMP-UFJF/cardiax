#ifndef TORORD_LAND_H_
#define TORORD_LAND_H_

#include <cmath>

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

//! Rush-Larsen para portas de Hodgkin-Huxley:  dy/dt = (y_inf - y)/tau
//! Solucao exata do passo congelando y_inf e tau:
//!     y(t+dt) = y_inf - (y_inf - y)*exp(-dt/tau)
static inline double rush_larsen(double y, double y_inf, double tau, double dt)
{
  return y_inf - (y_inf - y) * std::exp(-dt / tau);
}

//! Rush-Larsen na forma generica  dy/dt = a*y + b  (a < 0).
//! Cai para Euler explicito quando |a| e pequeno demais (evita 0/0).
static inline double rush_larsen_ab(double y, double a, double b, double dt)
{
  const double TOL = 1.0e-8;
  if (std::fabs(a) < TOL)
    return y + dt * (a * y + b);
  return std::exp(a * dt) * (y + b / a) - b / a;
}

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
