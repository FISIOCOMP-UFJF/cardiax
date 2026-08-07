// ----------------------------------------------
//  cg.cpp
// 
//  * Conjugate Gradient-Solver and Preconditioned CG-Solver 
//                         A x = f 
//    with A symmetric positive definite matrix 
//  * provides a diag-preconditioner 
// ----------------------------------------------

#include "cg.hpp"
#include "linalg.h"
#include "spmat.hpp"
#include "matrixops.hpp"

/*

int cg_solve(const SparseMatrix & A, const Vector & f, Vector & x, 
	     int it, int pr, int maxit, Real tol)
{
  int n = A.m();
  Vector r(n), p(n), Ap(n);
  Real alpha, beta;  
  
  // inital value for x
  x = 0;   
  r = f - A * x; 
  p = r;
  
  Real rtr = sprod(r,r), rtr1 = rtr; // dummy init of rtr1
  Real err_0 = sqrt(rtr);
  Real ptAp, err_k;

  int k;

  for(k=0; k<maxit; k++)
  { 
    Ap = A*p;
    ptAp = sprod(p,Ap); 
    alpha = rtr/ptAp; 
    
    x += alpha * p;
    r -= alpha * Ap;
    
    rtr1 =  sprod(r,r);  
    err_k = sqrt(rtr1); 
         
    if(err_k/err_0 < tol) break; 
           
    beta = rtr1/rtr;   
    rtr = rtr1;

    p = beta * p; 
    p += r;
  }

  if(it%pr==0)
    cout << "CG Iterations: " << setw(2) << k+1 
	 << " Residual norm: " <<  scientific << sqrt(rtr1/err_0);

  return (k+1);
}

int pcg_solve(const SparseMatrix & A, const SparseMatrix & invC, 
	      Vector & f, Vector & x, int it, int pr, int maxit,
	      Real tol)
{
  int n = A.n();
  x.resize(n);

  // To use the previous known solution, if transient, comment
  // the line below
  x=0;
   
  Vector z(n), r(n), p(n), Ap(n);
  r = f - A * x;
  z = invC * r; 
  p = z; 
  Real rtz0 = sprod(z,r), rtz = rtz0, rtz1 = rtz0;
  Real alpha, beta;  

  int k;
  for(k=0; k<maxit; k++)
  {    
    Ap = A * p;  
    alpha = rtz/sprod(Ap,p);
    x += alpha * p;
    r -= alpha * Ap; 
    
    z = invC*r;
    rtz1 = sprod(r,z);    
    if(sqrt(rtz1/rtz0) < tol) break;

    beta = rtz1/rtz;

    //p = p * beta + z; ---> This expression template is not working yet.
    
    p *= beta; 
    p += z; 

    rtz=rtz1; 
  }
  
  if(it%pr==0)
    cout << "PCG Iterations: " << setw(2) << k+1 
	 << " Residual norm: " << scientific << sqrt(rtz1/rtz0);

  return (k+1);
}

void calc_diag_prec(const SparseMatrix & A, SparseMatrix & invC)
{
  invC.resize(A.n(),A.m(),1);
  
  for(int i=0; i<invC.n();i++){
    invC(i,i) = 1./A(i,i);
  }
  invC.compress();
}

*/

