#ifndef TENTUSSCHERII_H
#define TENTUSSCHERII_H

#include <cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

//#define HF 

using namespace std;

/*!
 *  ten Tusscher II model
 *  ten Tusscher KH, Panfilov AV. 
 *  Ref.:"Alternans and spiral breakup in a human ventricular tissue model."
 *        Am J Physiol Heart Circ Physiol. 2006; 291: H1088-H1100
*/
class TenTusscher2006 : public CellModel
{
 public:     

  //! Default constructor
  TenTusscher2006();
  
  //! Set initial conditions for TT2 model
  virtual void init(double * values) const;

  //! Compute RHS equations for the TT2 model
  virtual void equation(const double t, const double * sv, double * values);

};

#endif
