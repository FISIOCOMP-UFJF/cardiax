#ifndef PARDISO_SOLVER_HPP
#define	PARDISO_SOLVER_HPP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "petsc_matrix.hpp"
#include "../util/util.hpp"

using namespace std;

//! PARDISO prototypes

extern "C" void pardisoinit(void *, int *, int *, int *, double *, int *);
extern "C" void pardiso(void *, int *, int *, int *, int *, int *, 
                        double *, int *, int *, int *, int *, int *,
                        int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec(int *, int *, double *, int *);

/*!
 * Wrapper class around PARDISO solver
 */
class PardisoSolver
{
public:
    
  PardisoSolver();
  ~PardisoSolver();
      
  typedef petsc::Matrix SparseMatrix;
     
  void init();
  void factor();
  void back_solve(double * x, double * b);
  void destroy();  
  void solve(SparseMatrix & A, double * x, double * b);
  
  SparseMatrix * get_matrix() { return pA; };
  void set_matrix(SparseMatrix * p) { pA = p; }
  void update_matrix(SparseMatrix * p);
  
private:   
  
  //! Pointer to the sparse matrix of the system
  SparseMatrix * pA;
    
  // PARDISO stuff ------------------------------------------------------------
  
  //! Pardiso control parameters
  int iparm[64];  
  int maxfct, mnum, phase, error, msglvl, solver;
  double dparm[64];
  
  //! Matrix data
  int mtype;
  int n, nnz, nrhs;
  int *ia, *ja;
  double *aa;
  
  //! Internal solver memory pointer pt  
  void *pt[64];  

  //! Number of processors
  int num_procs;
  
  // Dummy variables
  int idum;
  double ddum;
  
  //! Check error from PARDISO
  void checkerr(int error)
  {
    if (error != 0)
    {
      switch(-error)
      {
        case 1:   printf("Input inconsistent."); break;
        case 2:   printf("Not enough memory."); break;
        case 3:   printf("Reordering problem."); break;
        case 4:   printf("Zero pivot, numerical fact or IR problem."); break;
        case 5:   printf("Unclassified internal error."); break;
        case 6:   printf("Preordering failed."); break;
        case 7:   printf("Diagonal matrix problem."); break;
        case 8:   printf("32-bit integer overflow problem."); break;
        case 10:  printf("No license file pardiso.lic found."); break;
        case 11:  printf("License is expired."); break;
        case 12:  printf("Wrong username or hostname."); break;
        case 100: printf("Max num Krylov iterations."); break;
        case 101: printf("No sufficient convergence in Krylov."); break;
        case 102: printf("Error in Krylov iteration."); break;
        case 103: printf("Break-down in Krylov iteration."); break;
      }      
      printf("\nERROR in PARDISO Solver: %d", error);
      exit(1);      
    }    
  }
    
};

#endif	/* PARDISO_SOLVER_HPP */

