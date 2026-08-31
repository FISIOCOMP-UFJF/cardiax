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

//! Active stiffness Ka of the Land submodel, in kPa.
//!
//! Regazzoni & Quarteroni, CMAME 373 (2021) 113506, Eq. (50). Defining the
//! active stiffness as Ka = d(dTa/dt)/d(dlambda/dt) and applying it to
//!     Ta = h(lambda)*(Tref/dr)*[ (ZETAS+1)*XS + ZETAW*XW ]
//! gives the same expression with (ZETAS+1) -> As and ZETAW -> Aw:
//!     Ka = h(lambda)*(Tref/dr)*[ As*XS + Aw*XW ].
//! In this model As = Aw = A, so the bracket collapses to A*(XS + XW).
//!
//! This is the coefficient of the stabilisation term that the mechanics adds
//! to the active tension. It is NOT a new model parameter: it is fully
//! determined by the Land equations already implemented in equation().
//!
//! The constants mirror the "intact" branch (mode = 0) at the top of the Land
//! section in torord_land.cpp. They live here, with the cell model, rather
//! than leaking Land parameters into the mechanics module.
static inline double land_active_stiffness(double lambda, double XS, double XW)
{
  const double dr         = 0.25;
  const double wfrac      = 0.5;
  const double TOT_A      = 25.0;
  const double beta_0     = 2.3;
  const double Tref       = 120.0;   // mode 0 ("intact"); mode 1 uses 40.5
  const double lambda_min = 0.87;
  const double lambda_max = 1.2;

  const double A = (0.25 * TOT_A) / ((1.0 - dr) * wfrac + dr) * (dr / 0.25);

  const double lambda0 = std::fmin(lambda_max, lambda);
  const double Lfac = std::fmax(0.0, 1.0 + beta_0 * (lambda0
                        + std::fmin(lambda_min, lambda0) - (1.0 + lambda_min)));

  // equation() reads XS and XW through fmaxf(0, .); do the same here. A
  // negative crossbridge fraction would give a negative stiffness, which
  // would destabilise the scheme instead of stabilising it.
  const double xs = std::fmax(0.0, XS);
  const double xw = std::fmax(0.0, XW);

  return Lfac * (Tref / dr) * A * (xs + xw);
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
