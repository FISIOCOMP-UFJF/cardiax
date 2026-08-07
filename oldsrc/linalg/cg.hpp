// ----------------------------------------------
//  cg.hpp
// 
//  * provides ConjugateGradient-Solver and Preconditioned CG-Solver 
//                A x = f 
//    with A symmetric positive definite matrix 
//  * provides a diag-preconditioner 
// ----------------------------------------------

#ifndef CG_HPP
#define CG_HPP 
   
#include "vector.hpp" 
#include "spmat.hpp"
#include "array.hpp"

/*   
/// Conjugate Gradient Solver for A x = f
/// return number of iterations    
int cg_solve(const SparseMatrix & A, const Vector & f, Vector & x, 
	     int it, int pr, int maxit=1000, Real tol = 1e-6); 

/// Preconditioned conjugate Gradient Solver for A x = f 
/// with Preconditioner invC 
/// return-value is number of pcg-iterations
int pcg_solve(const SparseMatrix & A, const SparseMatrix & invC, Vector & f, 
	      Vector & x, int it, int pr, int maxit=1000, Real tol = 1e-6); 

/// Calculates Diagonal Preconditioner invC = diag(A)^{-1}
void calc_diag_prec(const SparseMatrix & A, SparseMatrix & invC);
*/

#endif
