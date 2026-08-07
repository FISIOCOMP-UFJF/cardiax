#ifndef BFGS_HPP
#define	BFGS_HPP

#include "nonlinear_solver.hpp"
#include "nonlinear_problem.hpp"
#include "../linalg/pardiso_solver.hpp"

/*!
 * BFGS Solver
 * 
 * Some references and notes about implementation:
 *   1. "The solution of nonlinear finite element equations"
 *       H. Matthies and G. Strang, 1979, IJNME, vol 14.
 */

class BFGSSolver : public NonlinearSolver
{
 public:
    
  BFGSSolver(NonlinearProblem * p);
  ~BFGSSolver();
  
  void init(); // not in use
  void init(int neq, PardisoSolver * ps);
  double line_search(double s);
  petsc::Vector & residual0() { return aux_r0; }
  int solve();
  void solve_eqs(double * x, double * b);
  bool update(double s, arma::vec & UI, arma::vec & R0, arma::vec & R1);

 protected:
     
  int neq;       //!< Number of equations
  int nits;      //!< Number of iterations
  int nmax;      //!< Maximum number of iterations
  
  int maxups;    //!< Max num of QN iterations between updates
  int maxref;    //!< Max num of updates (reformations)
  double cmax;   //!< Maximum value for the condition number
  
  double ls_min; //!< Minimum LS step
  double ls_tol; //!< LS tolerance
  int ls_mxiter; //!< Max number of LS iterations
  
  int nups;      //!< Number of stiffness updates (counter)
  
  arma::vec ui;  //!< Displacement increment 
  arma::vec r0;  //!< Residual at step k-1
  arma::vec r1;  //!< Residual at step k
    
  //! BFGS stuff
  arma::mat V;   //!< BFGS vector
  arma::mat W;   //!< BFGS vector
  arma::vec D;   //!< temp
  arma::vec G;   //!< temp
  arma::vec H;   //!< temp
  
  petsc::Matrix * K0;
  petsc::Vector aux_r0, aux_r1, aux_ui;
  
  PardisoSolver * parsolver;
    
};

#endif	/* BFGS_HPP */

