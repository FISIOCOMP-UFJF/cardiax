#ifndef LUORUDYI_H
#define LUORUDYI_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

using namespace std;

/*! Luo-Rudy Phase I model
 *   Ching-hsing Luo and Yoram Rudy
 *   "A Model of the Ventricular Cardiac Action Potential.
 *   Depolarization, repolarization and their interaction"
 *   Circulation Research  , 68, 1501-1526. PubMed ID: 1709839
*/
class LuoRudy : public CellModel
{
 public:
  
  //! Default constructor
  LuoRudy();

  //! Set initial conditions for LRI model
  virtual void init(double * values) const;
  
  //! Compute RHS equations for the LRI model
  virtual void equation(const double time, const double * sv, double * values);

};

#endif
