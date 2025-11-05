#ifndef TENTUSSCHER2004_H
#define TENTUSSCHER2004_H

#include <cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

//#define HF

using namespace std;

class TenTusscher2004 : public CellModel
{
 public:     

  //! Default constructor
  TenTusscher2004();
  
  //! Set initial conditions for TT2 model
  virtual void init(double * values) const;

  //! Compute RHS equations for the TT2 model
  virtual void equation(const double t, const double * sv, double * values);

};

#endif
