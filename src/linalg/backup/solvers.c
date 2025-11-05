#include "solvers.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* CSR solvers */

int
csr_pcg (int n, csr_matrix *A, Real *invM, Real *y, Real *x, 
	  int it, int printRate)
{
  int i;
  Real *r, *p, *Ap, *z;
  Real alfa, beta;
  Real rtz, rtz0, rtz1;
  
  r  = new_array<Real>(n);
  p  = new_array<Real>(n);
  Ap = new_array<Real>(n);
  z  = new_array<Real>(n);
  
  //vectorSet(n, x, 0.0);          // x = 0
  //vectorCpy(n, r, y);            // r = y
  csr_spmv (A, x, r);              // Ap = A * p
  vectorSub(n, y, r, r);
  vectorMul(n, z, invM, r);        // z = invM * r
  vectorCpy(n, p, z);              // p = z

  rtz0 = scalarProd (z,r,n);
  rtz  = rtz0;

  for(i=0; i<n; i++)
  {
    csr_spmv (A, p, Ap);            // Ap = A * p

    alfa = rtz / scalarProd(Ap,p,n);

    vectorSAXPY (n,  alfa, p,  x);  // x = x + alfa * p
    vectorSAXPY (n, -alfa, Ap, r);  // r = r - alfa * Ap

    vectorMul (n, z, invM, r);      // z = invM * r
    
    rtz1 = scalarProd(r,z,n);
    if(sqrt(rtz1/rtz0) < CGTOL) break;

    beta = rtz1 / rtz;

    vectorSXPAY (n, beta, z, p);    // p = z + beta * p
    
    rtz = rtz1;
  }

  if(it%printRate==0)
    fprintf(stdout,"CPU_PCG_CSR: it %d rnorm %6.4e ",i,sqrt(rtz1/rtz0));

  free(z);
  free(r);
  free(p); 
  free(Ap);

  return 0;
}

/* ELLPACK solvers */

int
ell_pcg1d (int n, int max_row, int *Aj, Real *Ax, Real *inv,
	  Real *y, Real *x, int it, int printRate)
{
  int i;
  Real *r, *p, *Ap, *z;
  Real alfa, beta;
  Real rtz, rtz0, rtz1;
  
  r  = new_array<Real>(n);
  p  = new_array<Real>(n);
  Ap = new_array<Real>(n);
  z  = new_array<Real>(n);
  
  //vectorSet(n, x, 0.0);          // x = 0
  //vectorCpy(n, r, y);            // r = y

  ell_spmv1d (n, max_row, Aj, Ax, x, r);    // r = A * x

  vectorSub(n, y, r, r);
  vectorMul(n, z, inv, r);                  // z = invM * r
  vectorCpy(n, p, z);                       // p = z

  rtz0 = scalarProd (z,r,n);
  rtz  = rtz0;

  for(i=0; i<n; i++)
  {
    ell_spmv1d (n, max_row, Aj, Ax, p, Ap); // Ap = A * p

    alfa = rtz / scalarProd(Ap,p,n);

    vectorSAXPY (n,  alfa, p,  x);          // x = x + alfa * p
    vectorSAXPY (n, -alfa, Ap, r);          // r = r - alfa * Ap

    vectorMul (n, z, inv, r);               // z = invM * r
    
    rtz1 = scalarProd(r,z,n);
    if(sqrt(rtz1/rtz0) < CGTOL) break;

    beta = rtz1 / rtz;

    vectorSXPAY (n, beta, z, p);            // p = z + beta * p
    
    rtz = rtz1;
  }

  if(it%printRate==0)
    fprintf(stdout,"CPU_PCG_ELL: it %d rnorm %6.4e ",i,sqrt(rtz1/rtz0));

  free(z); free(r); free(p); free(Ap);

  return 0;
}

int ell_cg (int n, int num_row_per_col, int **Aj, Real **Ax, 
	    Real *y, Real *x, int it, int printRate)
{
  int i,k;
  Real r[n]; // FIX
  Real z[n]; // FIX
  Real w[n]; // FIX
  Real alfa, beta, aux;
  Real deltaNew, deltaOld, deltaO;

  vectorSet(n, x, 0.0);
  
  for(k=0; k<n; k++){
    aux  = y[k];
    r[k] = aux;
    w[k] = aux; 
  }

  deltaNew = scalarProd (r,r,n);
  deltaO   = deltaNew;

  i = 0;

  while( (i<n) && (deltaNew > CGTOL*CGTOL )){
    
    ell_spmv (n, 9, Aj, Ax, w, z); // z = A * w

    alfa = deltaNew / scalarProd(w,z,n);

    for(k=0; k<n; k++){
      x[k] = x[k] + alfa * w[k];
      r[k] = r[k] - alfa * z[k];
    }
   
    deltaOld = deltaNew;
    deltaNew = scalarProd(r,r,n);

    beta = deltaNew / deltaOld;

    for(k=0; k<n; k++){
      w[k] = r[k] + beta * w[k];
    }
    i++;
  }

  if(it%printRate==0) 
    fprintf(stdout,"CPU_CG_ELL: it %d rnorm %6.4e ",i,deltaNew);

  return 0;
}

int ell_cg1d (int n, int num_row_per_col, int *Aj, Real *Ax, 
	      Real *y, Real *x, int it, int printRate)
{
  int i,k;
  Real r[n]; // FIX THIS
  Real z[n]; // FIX
  Real w[n]; // FIX
  Real alfa, beta, aux;
  Real deltaNew, deltaOld, deltaO;

  vectorSet(n, x, 0.0);
  
  for(k=0; k<n; k++){
    aux  = y[k];
    r[k] = aux;
    w[k] = aux; 
  }

  deltaNew = scalarProd (r,r,n);
  deltaO   = deltaNew;

  i = 0;

  while( (i<n) && (deltaNew > CGTOL*CGTOL )){
    
    ell_spmv1d (n, 9, Aj, Ax, w, z); // z = A * w

    alfa = deltaNew / scalarProd(w,z,n);

    for(k=0; k<n; k++){
      x[k] = x[k] + alfa * w[k];
      r[k] = r[k] - alfa * z[k];
    }
   
    deltaOld = deltaNew;
    deltaNew = scalarProd(r,r,n);

    beta = deltaNew / deltaOld;

    for(k=0; k<n; k++){
      w[k] = r[k] + beta * w[k];
    }
    i++;
  }

  if(it%printRate==0) 
    fprintf(stdout,"CPU_CG_ELL_OK: it %d rnorm %6.4e ",i,deltaNew);

  return 0;
}

/* DIAGONAL format solvers */

int dia_cg (int n, int num_row_per_col, int *Aoff, Real **Adiag,
	    Real *y, Real *x, int it, int printRate, 
	    int *num_its, double *error)
{
  int i,k;
  Real r[n];
  Real z[n];
  Real w[n];
  Real alfa, beta, aux;
  Real deltaNew, deltaOld, deltaO;

//  vectorSet(n, x, 0.0);
  
//  for(k=0; k<n; k++){
//    aux  = y[k];
//    r[k] = aux;
//    w[k] = aux;
//  }

  dia_spmv(n, num_row_per_col, Aoff, Adiag, x, w);
  vectorSub(n, y, w, r);

  for(k=0; k<n; k++){
	  w[k] = r[k];
  }


  deltaNew = scalarProd (r,r,n);
  deltaO   = deltaNew;

  i = 0;

//  fprintf(stdout,"Residuo CPU = %e\n", deltaNew);
  while( (i<n) && (deltaNew > CGTOL*CGTOL )){

    dia_spmv (n, 9, Aoff, Adiag, w, z); // z = A * w

    alfa = deltaNew / scalarProd(w,z,n);

    for(k=0; k<n; k++){
      x[k] = x[k] + alfa * w[k];
      r[k] = r[k] - alfa * z[k];
    }
   
    deltaOld = deltaNew;
    deltaNew = scalarProd(r,r,n);

    beta = deltaNew / deltaOld;

    for(k=0; k<n; k++){
      w[k] = r[k] + beta * w[k];
    }
    i++;
  }

  *error = sqrt(deltaNew);
  *num_its = i;

  return 0;
}

int
dia_pcg (int n, int num_row_per_col, int *Aoff, 
	 Real **Adiag, Real *invM, Real *y, Real *x,
	 int it, int printRate, int *num_its, double *error)
{
  int i;
  Real r[n];
  Real z[n];
  Real p[n];
  Real Ap[n];

  Real alfa, beta;
  Real rtz, rtz0, rtz1;

  dia_spmv(n, num_row_per_col, Aoff, Adiag, x, p);

  vectorSub(n, y, p, r);
  vectorMul(n, z, invM, r);      // z = invM * r
  vectorCpy(n, p, z);            // p = z

  rtz0 = scalarProd (z,r,n);
  rtz  = rtz0;

  for(i=0; i<n; i++)
  {
    dia_spmv (n, 9, Aoff, Adiag, p, Ap); // Ap = A * p
    alfa = rtz / scalarProd(Ap,p,n);
    vectorSAXPY (n,  alfa, p,  x);       // x = x + alfa * p
    vectorSAXPY (n, -alfa, Ap, r);       // r = r - alfa * Ap
    vectorMul (n, z, invM, r);           // z = invM * r
    rtz1 = scalarProd(r,z,n);

    if(sqrt(rtz1/rtz0) < CGTOL) break;

    beta = rtz1 / rtz;
    vectorSXPAY (n, beta, z, p);         // p = z + beta * p
    rtz = rtz1;
  }

  *error = sqrt(rtz1/rtz0);
  *num_its = i;

  return 0;
}













