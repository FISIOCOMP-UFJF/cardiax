#ifndef NLPROBLEM_HPP
#define NLPROBLEM_HPP

#include "linalg/petsc_matrix.hpp"
#include "linalg/petsc_vector.hpp"

class NonlinearProblem
{
public:
  
  NonlinearProblem() {}
  virtual ~NonlinearProblem() {};

  //! Compute the residual r (or evaluates F at current point x)
  virtual void evaluate(petsc::Vector & r) = 0;

  //! Compute Jacobian K = F' at current point x
  virtual void jacobian(petsc::Matrix & K) = 0;
  
  //! Update state of system (x_k = x_k-1 + u)
  virtual void update(petsc::Vector & u, double s) = 0;

  //! Decide whether solution has converged or not
  virtual bool converged(petsc::Vector & du) = 0;
  
  //! Return size of the problem
  virtual uint size() = 0;

};

#endif
