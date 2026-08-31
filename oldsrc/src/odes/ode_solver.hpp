#ifndef ODESOLVER_H
#define ODESOLVER_H

#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include "util/util.hpp"

// Forward declaration
class CellModel;

/** 
    Base class for ODE solvers
*/
class ODESolver
{
protected:  
  
  // The ODE model to solve (mainly a CellModel object)
  CellModel * ode;

public:
  
  // Constructor
  ODESolver (CellModel * odep) : ode(odep) {}
  
  // Virtual destructor
  virtual ~ODESolver() {}

  // Interface for ODE solvers
  virtual void advance(double * y, double & t, double & dt) = 0;

  // Factory method to create ODE solvers
  static ODESolver * create(std::string method, CellModel *cell);

};

#endif
