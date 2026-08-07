#include "pardiso_solver.hpp"

PardisoSolver::PardisoSolver()
{
  // numbers of processors, value of OMP_NUM_THREADS
  char *var;
  var = getenv("OMP_NUM_THREADS");
  if (var != NULL)
    sscanf(var, "%d", &num_procs);
  else
  {
    printf("Set environment OMP_NUM_THREADS to 1");
    exit(1);
  }
}

PardisoSolver::~PardisoSolver()
{
  delete [] ia;
  delete [] ja;
  delete [] aa;
}

void PardisoSolver::back_solve(double * x, double * b)
{ 
   pardiso_chkvec (&n, &nrhs, b, &error);
   checkerr(error); 
 
  // Back substitution and iterative refinement

  phase = 33;
  iparm[7] = 1; // Max numbers of iterative refinement steps

  pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
          &n, aa, ia, ja, &idum, &nrhs,
          iparm, &msglvl, b, x, &error, dparm);
  checkerr(error);

}

void PardisoSolver::destroy()
{
  // Release internal memory
  phase = -1; 

  pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
          &n, &ddum, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error, dparm);
}

void PardisoSolver::factor()
{
  pardiso_chkmatrix  (&mtype, &n, aa, ia, ja, &error);
  checkerr(error); 

  // Reordering and Symbolic Factorization.  This step also allocates
  // all memory that is necessary for the factorization.
  phase = 11;

  pardiso(pt, &maxfct, &mnum, &mtype, &phase,
          &n, aa, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error, dparm);
  checkerr(error);
   
  // Numerical factorization
  phase = 22;

  pardiso(pt, &maxfct, &mnum, &mtype, &phase,
          &n, aa, ia, ja, &idum, &nrhs, 
          iparm, &msglvl, &ddum, &ddum, &error, dparm);
  checkerr(error);   
}

void PardisoSolver::init()
{
  mtype    = -2;        // real symmetric matrix
  iparm[0] = 0;         // default parameters  
  iparm[2] = num_procs; // procs for OpenMP
  solver   = 0;         // use sparse direct solver
  error    = 0;
  
  pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
  checkerr(error); 
  
  maxfct = 1; // Maximum number of numerical factorizations
  mnum   = 1; // Which factorization to use  
  msglvl = 0; // Print statistical information -> 0 No information
  error  = 0; // Initialize error flag
  
  // matrix information
  n    = pA->size();          
	nnz  = pA->get_nnz();
	nrhs = 1;
  
  // define matrix
  ia = new int[n+1];
  ja = new int[nnz];
  aa = new double[nnz];  
  
  pA->get_CSR(&n,ia,ja,aa);

//  printf("n   ---> %d\n", n);
//  printf("nz  ---> %d\n", nnz);
//  printf("ia  ---> ");
//  for(int i=0; i<n+1; i++) printf("%d ", ia[i]);
//  printf("\n");
//  printf("ja  ---> ");
//  for(int i=0; i<nnz; i++) printf("%d ", ja[i]);
//   printf("\n");
//  printf("aa  ---> ");
//  for(int i=0; i<nnz; i++) printf("%lf ", aa[i]);
//  printf("\n");
}

void PardisoSolver::solve(SparseMatrix & A, double *x, double *b)
{
  printf("PARDISO solving ...\n"); 
  set_matrix(&A);
  printf("Creating matrix completed ...\n");
  init(); 
  printf("Initialization completed ...\n");             
  
  pardiso_chkmatrix  (&mtype, &n, aa, ia, ja, &error);
  checkerr(error);
  
  factor();
  printf("Factorization completed ...\n");
  back_solve(x,b);
  printf("Solve completed ... \n");
  
//  printf("The solution of the system is: ");
//  for (int i = 0; i < n; i++) {
//    printf("\n x [%d] = % lf", i, x[i]);
//  }
//  printf("\n");
}

void PardisoSolver::update_matrix(SparseMatrix * p)
{
  int n = pA->size();    
  pA = p;
  pA->get_CSR(&n,ia,ja,aa);
}