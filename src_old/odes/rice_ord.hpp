#ifndef RICE_ORD_HPP
#define	RICE_ORD_HPP

#include <cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

class RiceOHaraRudy : public CellModel
{

public:

  //! Default constructor (initializes constants ans some other vars)
  RiceOHaraRudy();

  //! Set initial conditions for ORd model
  virtual void init(double * values) const;

  //! Compute RHS equations for the ORd model
  virtual void equation(const double t, const double * sv, double * values);
  
private:
    
    double active_force;
 
};

#endif

