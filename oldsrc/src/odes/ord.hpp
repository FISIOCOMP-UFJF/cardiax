#ifndef ORD_HPP
#define	ORD_HPP

#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

class OHaraRudy : public CellModel
{

public:

  //! Default constructor (initializes constants ans some other vars)
  OHaraRudy();

  //! Set initial conditions for ORd model
  virtual void init(double * values) const;

  //! Compute RHS equations for the ORd model
  virtual void equation(const double t, const double * sv, double * values);
 
};

#endif	/* ORD_HPP */

