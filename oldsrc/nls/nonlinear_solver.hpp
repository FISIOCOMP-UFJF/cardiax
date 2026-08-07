#ifndef NLSOLVER_HPP
#define NLSOLVER_HPP

#include <map>
#include "nonlinear_problem.hpp"
#include "linalg/petsc_linear_solver.hpp"

//#define USE_PARDISO

/*! This class implements the Newton method for solving
 *  system of nonlinear equations. 

    It can use
       -> TODO: PETSC KSP framework 
    for solving the linear systems within the nonlinear problem.
*/

class NonlinearSolver
{
public:
  
  //! Constructor
  NonlinearSolver(){};
    
  //! Constructor with a nonlinear problem
  NonlinearSolver(NonlinearProblem * p);

  //! Destructor 
  ~NonlinearSolver();

  //! Initialization
  virtual void init() = 0;
  
  //! Solve the system
  virtual int solve() = 0;

    virtual void ViewMatrix()=0;
  
  //! Change absolute tolerance
  void set_absolute_tol(double tol) { atol = tol; }
  
  //! Change relative tolerance
  void set_relative_tol(double tol) { rtol = tol; }  

  TimerSection timer;
  
protected:

  //! Relative tolerance
  double rtol;
  
  //! Absolute tolerance
  double atol;
     
  //! The nonlinear problem
  NonlinearProblem * nlp;
  
  //! KSP solver
  petsc::LinearSolver * ls;

};


#endif
