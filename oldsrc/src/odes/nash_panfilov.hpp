#ifndef NASH_PANFILOV_H_
#define NASH_PANFILOV_H_

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/*!
 *  Nash-Panfilov modified version of the Aliev-Panfilov model
 * 
 *  Original paper: 
 *  "Electromechanical model of excitable tissue to study cardiac arrhythmias"
 *   Progress in Biophysics & Molecular Biology, 85, (2004), 501-522
 *
 *  Parameters for simulation:
 *  - Vm [0,1]
 *  - Stimulus 0.5 should work
 *
 *  Parameters from Physical Review Letters 95, 258104 (2005)
 *   "Self Organized Pacemakers in a Coupled Reaction-Diffusion-Mechanics System"
 *  A. V. Panfilov and R. H. Kelderman
*/
class NashPanfilov : public CellModel
{
 public:
  
  //! Default constructor
  NashPanfilov();

  //! Set initial conditions for NP model
  virtual void init(double * values) const;

  //! Compute RHS equations for the NP model
  virtual void equation(const double t, const double * sv, double * values);

 private:

  //! Rate switch function proposed by Nash & Panfilov (2004)
  double rateswitchNP(double v);

  //! Rate switch function proposed by Goektepe & Kuhl (2010)
  double rateswitchGK(double v);

  //! Rate switch function proposed by Eriksson et. al (2013)
  double rateswitchEr(double v);

};

#endif
