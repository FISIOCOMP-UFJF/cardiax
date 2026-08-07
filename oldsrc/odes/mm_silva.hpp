#ifndef MM_SILVA_H_
#define MM_SILVA_H_

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/*!
 *  Minimal Model + Joao Gabriel TA
*/
class MMSilva : public CellModel
{
 public:
  
  //! Default constructor
  MMSilva();

  //! Set initial conditions for NP model
  virtual void init(double * values) const;

  //! Compute RHS equations for the NP model
  virtual void equation(const double t, const double * sv, double * values);

 private:
 
   double H(double x, double y);

   double vnoinf(double x, double y);

};

#endif
