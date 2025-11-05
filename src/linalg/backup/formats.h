/**
 * @file   formats.h
 * @author rocha <rocha@rocha-laptop>
 * @date   Tue Feb 23 18:07:03 2010
 * 
 * @brief  Sparse matrix formats
 *           Supported formats: CSR, ELL and DIA
 *
 */

#ifndef FORMATS_H
#define FORMATS_H

#include "../util/util.h"
#include "../util/datatype.h"
#include "../util/malloc.h"

// ===================================================================
// Generic sparse matrix 
// ===================================================================

// supported storage for sparse matrices
enum sparse_format_t { _CSR, _DIA, _ELL };

typedef enum sparse_format_t sparse_format;

struct sparse_matrix_t {
  sparse_format format;
  void *repr;
};

typedef struct sparse_matrix_t sparse_matrix;

sparse_matrix *
create_sparse_matrix (sparse_format format, void* repr);

void
delete_sparse_matrix (sparse_matrix* A);

// ===================================================================
// triplet storage
// ===================================================================

struct triplet_type {
  int size;
  int  *i;
  int  *j;
  Real *x;
};

typedef struct triplet_type triplet;

triplet *
triplet_create (int size);

void
triplet_print  (triplet *T, int nen, int nel);

void
triplet_delete (triplet *T);

// ===================================================================
// CSR sparse matrix
// ===================================================================

struct csr_matrix_type {
  int m;       // number of rows
  int n;       // number of cols
  int nnz;     // number of stored nonzero entries (if compr)
               // otherwise the size of the triplet
  int nzmax;   // maximum of nz allowed
  Real *vals;  // stored values
  int  *cols;  // columns
  int  *ptrs;  // pointers to ax and ja
  int cmpr;    // is compressed ?
  triplet *T;  // triplet if it is not compressed
};

typedef struct csr_matrix_type csr_matrix;

csr_matrix *
csr_alloc (const int m, const int n, const int nnz, int triplet);

csr_matrix *
csr_create (const int m, const int n, const int nnz);

csr_matrix *
csr_pack (const int m, const int n, const int nnz, Real *ax, int *ja, int *ia);

void
csr_delete (csr_matrix *A);

int
csr_count_nnz (int m, int n, Real **K);

void
csr_createmat (int n, Real **A, Real *ele, int *dsp, int *col);

void
csr_inv_diag (int n, csr_matrix *A, Real *invA);

int
csr_entry_ij (csr_matrix *A, int i, int j, Real x);

csr_matrix *
csr_compress (csr_matrix *A);

int
csr_sumdup (csr_matrix *A);

// ===================================================================
// ELLPACK sparse matrix
// ===================================================================

struct ell_matrix_type {
  int m;
  int n;
  int nnz;
  int maxcol;
  int **cols;
  Real **vals;
};

typedef struct ell_matrix_type ell_matrix;

ell_matrix *
create_ell_matrix (const int m, const int n, const int nnz, 
		   const int maxcol);

ell_matrix *
pack_ell_matrix (const int m, const int n, const int nnz, 
		 const int maxcol, int **ja, Real **ax);

void
delete_ell_matrix (ell_matrix *A);

// ===================================================================
// DIA sparse matrix
// ===================================================================

struct dia_matrix_type {
  int m;
  int n;
  int nnz;
  int numdiag;
  int *offs;
  Real **diag;
};

typedef struct dia_matrix_type dia_matrix;

dia_matrix *
create_dia_matrix (const int m, const int n, const int numdiag);

dia_matrix *
pack_dia_matrix (const int m, const int n, const int numdiag, 
		 int *offs, Real **diag);

void
delete_dia_matrix (dia_matrix *A);

// ===================================================================
// Utils - General purporse and conversion routines
// ===================================================================

void
infdia (int n, int *ja, int *ia, int *ind, int *idiag);

int
csrell (int nrow, Real *a, int *ja, int *ia, int maxcol, 
	Real **coef, int **jcoef, int *ndiag);

void
csrcoo (int n, int nz, const int *Ap, const  int *Aj, const Real *Ax, 
	int *row, int *col, Real *data);

csr_matrix *
coocsr (int m, int n, int nz, const int *I, const int *J, const Real *X);

void
csrdia (int n, int idiag, const Real *a, int *ja, int *ia, 
	int ndiag, Real **diag, int *ioff);

// ===================================================================
// SPMV Operations
// ===================================================================

//void
//spmv (sparse_matrix *A, real *x, real *y);

void
csr_spmv (csr_matrix *A, Real *x, Real *y);

void
dia_spmv (int n, int nd, int *ioff, Real **diag, Real *x, Real *y);

void
ell_spmv1d (int num_rows, int num_cols_per_row, int *Aj, Real *Ax, 
	    Real *x, Real *y);

void
ell_spmv (int num_rows, int num_cols_per_row, int **Aj, Real **Ax, 
	  Real *x, Real *y);

// ===================================================================
// Included by Ronan
// ===================================================================

int*
sort_dia (int *dia_off, int numdiags);

#endif
