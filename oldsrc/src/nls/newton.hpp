#ifndef NEWTON_H
#define	NEWTON_H

#include "nonlinear_solver.hpp"

class Newton : public NonlinearSolver
{
public:

  Newton(NonlinearProblem * p);
  virtual ~Newton(); 
  
  void init();
  int solve();
  void ViewMatrix();
   
  petsc::Vector & residual() const { return *r; }  
  
protected:
    
  int nits;
  int nmax;
  int pits; 
  petsc::Matrix * K;
  petsc::Vector * u;
  petsc::Vector * r;
};

#endif	/* NEWTON_H */

