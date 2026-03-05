#ifndef _solvers_h_
#define _solvers_h_

/**
 * @file   solvers.h
 * @author Bernardo M. Rocha <bernardo@lncc.br>
 * @date   Wed Mar  3 13:36:50 2010
 * 
 * @brief  Solvers for the sparse matrices formats
 * 
 * 
 */

#include "formats.h"
#include "../util/util.h"
#include "../util/datatype.h"

#define CGTOL 1.0e-8

void
print_timings_cpu();

void
initialize_solver(int ndof);

void
finalize_solver();


/* CSR matrix */

int
csr_pcg (int n, csr_matrix *A, Real *invM, 
	  Real *y, Real *x, int it, int printRate);

int
csr_cg (int n, int *A_dsp, int *A_col, Real *A_ele, 
	Real *y, Real *x, int printRate, int *num_its, double *error);

/* ELLPACK matrix */


int
ell_cg (int n, int num_row_per_col, int **Aj, Real **Ax, 
	Real *y, Real *x, int it, int printRate);  // 2D matrices

int
ell_cg1d (int n, int num_row_per_col, int  *Aj, Real  *Ax, 
	  Real *y, Real *x, int it, int printRate); // 1D linearized

int
ell_pcg1d (int n, int num_row_per_col, int *Aj, Real *Ax, Real *inv,
	  Real *y, Real *x, int it, int printRate); // 1D linearized

/* DIAgonal matrix */

int
dia_pcg (int n, int num_row_per_col, int *Aoff,	Real **Adiag,
	 Real *invM, Real *y, Real *x, int it, int printRate, 
	 int *num_its, double *error);

int
dia_cg (int n, int num_row_per_col, int *Aoff, Real **Adiag, 
	Real *y, Real *x, int it, int printRate, 
	int *num_its, double *error);

#endif
