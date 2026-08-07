/* 
 * A simple two dimensional example.
 * Taken from 
 *   "Solving nonlinear equations with Newton's method"
 *   C. T. Kelley, chapter 2, page 39
 */

#include "../nonlinear_problem.hpp"

#ifndef SIMPLE_HPP
#define	SIMPLE_HPP

class SimpleNLDemo : public NonlinearProblem
{
public:
 
  SimpleNLDemo();
  ~SimpleNLDemo();
  void evaluate(petsc::Vector & r);
  void jacobian(petsc::Matrix & K);
  void update(petsc::Vector & u, double s);
  bool converged(petsc::Vector & du);
  uint size() { return 2; }
  
private:
  double x[2];
      
};

#endif	/* SIMPLE_HPP */

