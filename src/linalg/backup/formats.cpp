// -------------------------------------------------------------------
// SPARSE MATRIX Formats
//
// CSR: compressed sparse row
// ELL: ELLPACK format
// COO: coordinate format
// DIA: diagonal format
// -------------------------------------------------------------------

#include "formats.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int max0 (int a, int b);
int min0 (int a, int b);

// -------------------------------------------------------------------
// SPMV (SParse Matrix Vector) routines
// -------------------------------------------------------------------

/**
 * CSR matrix-vector multiplication (y = Ax)
 *
 */
void csr_spmv (csr_matrix *A, Real *x, Real *y) 
{
  int n = A->n; // m = n 
  int i, j, k1, k2, q;
  int *cols, *ptrs;
  Real *vals, sum;
  cols = A->cols;
  ptrs = A->ptrs;
  vals = A->vals;

#pragma omp parallel for private(i,j,k1,k2,q,sum) \
  shared(n,cols,ptrs,vals,x,y)
  for(i=0; i<n; i++){
    k1 = ptrs[i];
    k2 = ptrs[i+1];
    sum = 0.0;
    for(j=k1; j<k2; j++){
      q = cols[j];
      sum  += vals[j] * x[q];
    }
    y[i] = sum;
  }
  // end of omp parallel for
}

/*!
 DIA format matrix-vector multiplication (y = Ax)
 */
void dia_spmv (int n, int ndiag, int *ioff, Real **diag, Real *x, Real *y)
{
  int j,k,io,i1,i2;

  vectorSet(n,y,0.0);

  for(j=0; j<ndiag; j++){
    io = ioff[j];
    i1 = max0(0,0-io); // 0 for c, 1 fortran
    i2 = min0(n,n-io);
    for(k=i1; k<=i2; k++){
      y[k] = y[k] + diag[j][k] * x[k+io];
    }
  }
}

/*!
 ELLPACK format matrix-vector multiplication (y = Ax)
 */
void ell_spmv (int n, int maxcol, int **Aj, Real **Ax, Real *x, Real *y)
{ 
  int i, j, col, val;
  vectorSet(n,y,0.0);
#pragma omp parallel for private(i,j,col,val)	\
  shared(n,Ax,Aj,x,y) schedule(static)
  for(j=0; j < maxcol; j++){
    for(i=0; i < n; i++){
    
      col = Aj[i][j];
      val = Ax[i][j];

      if (val != 0){
        y[i] += val * x[col];
      }
    }    
  }  
}

/*!
 ELLPACK format matrix-vector multiplication (y = Ax)
 in linearized form using a 1D array (identical to the CUDA Kernel)
 */
void ell_spmv1d (int n, int maxcol, int *Aj, Real *Ax, Real *x, Real *y)
{    
  vectorSet(n,y,0.0);
  for(int j=0; j < maxcol; j++){
    for(int i=0; i < n; i++){
      int   col = Aj[i + n * j];
      Real val = Ax[i + n * j];
      if (val != 0){
        y[i] += val * x[col];
      }
    }   
  }  
}

// -------------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------------

int max0 (int a, int b)
{ 
  if (a > b){ return a;} else { return b; }
}

int min0 (int a, int b)
{
  if (a < b){ return a;} else { return b; }
}

// -------------------------------------------------------------------
// Conversion functions
// -------------------------------------------------------------------

/*!
 Obtains information on the diagonals of A.
 This subroutine finds the lengths of each of the 2*n-1 diagonals of A
 it also outputs the number of nonzero diagonals found. 
 @param[in]  n       Dimension of the matrix A
 @param[in]  a       Matrix stored in CSR format
 @param[in]  ja      Matrix stored in CSR format
 @param[in]  ia      Matrix stored in CSR format
 @param[out] idiag   Integer number of nonzero diagonals found
 @param[out] ind     Integer array of length at least 2*n-1. The k-th 
                     entry in ind contains the number of nonzero 
                     elements in the diagonal number k, the 
                     numbering beeing from the lowermost diagonal
                     (bottom-left). In other words ind(k) = length of diagonal
                     whose offset wrt the main diagonal is = - n + k.
 
 Translated to C from SPARSKIT unary.f 
 SPARSKIT by Y. Saad
 */

void infdia (int n, int *ja, int *ia, int *ind, int *idiag)
{
  int i,j,k;
  int n2 = n+n-1;

  vectorSet(n2, ind, 0);

  for(i=0; i<n; i++){
    for(k=ia[i]; k<=ia[i+1]-1; k++){
      j = ja[k];
      ind[n+j-i] = ind[n+j-i] + 1;
    }
  }

  *idiag = 0;
  for(k=0; k<n2; k++){
    if(ind[k] != 0) *idiag = *idiag + 1;
  }
}

/*!
 Converts compressed sparse row to ellpack format.
 This subroutine converts  matrix stored in the general a, ja, ia 
 format into the coef, jcoef itpack format.
 @param[in]  nrow     Dimension of the matrix A
 @param[in]  a        Matrix stored in CSR format
 @param[in]  ja       Matrix stored in CSR format
 @param[in]  ia       Matrix stored in CSR format
 @param[in]  ncoef    First dimension of arrays coef and jcoef
 @param[in]  maxcol   First dimension of arrays coef and jcoef
 @param[out] coef     real array containing the values of the matrix A in
                      itpack-ellpack format
 @param[out] jcoef    Int array containing the column indices of coef(i,j) in A
 @param[out] ndiag    Number of active 'diagonals' found
 
 Translated to C from SPARSKIT formats.f
 SPARSKIT by Y. Saad
 */

int csrell (int nrow, Real *a, int *ja, int *ia, int maxcol, 
    Real **coef, int **jcoef, int *ndiag)
{
  int i, j, k, k1, k2;
  int ierr;

  // first determine the length of each row of lower-part-of(A)
  ierr = 0;
  *ndiag = 0;
  for(i=0; i<nrow; i++){
    k = ia[i+1]-ia[i];
    *ndiag = max0(*ndiag,k);
  }
  
  // check whether sufficient columns are available
  if(*ndiag > maxcol){
    ierr = 1; 
    return ierr;    
  }
  
  // fill coef with zero elements and jcoef with row numbers
  for(i=0; i<nrow; i++){
    for(j=0; j<*ndiag; j++){
      coef[i][j] = 0.0;
    }
    jcoef[i][j] = i;
  }
  
  // copy elements row by row
  for(i=0; i<nrow; i++){
    k1 = ia[i];
    k2 = ia[i+1];
    for(k=k1; k<k2; k++){ 
      coef[i][k-k1] = a[k];
      jcoef[i][k-k1] = ja[k];
    }
  }
  return 0;
}

/*!
 Converts CSR to COO format.
 @param[in]  n      Dimension of the matrix A
 @param[in]  nz     Number of non-zero entries
 @param[in]  Ap     Matrix stored in CSR format
 @param[in]  Aj     Matrix stored in CSR format
 @param[in]  Ax     Matrix stored in CSR format
 @param[out] row    Row array of COO format
 @param[out] col    Column array of COO format
 @param[out] data   Data array with non-zero entries of COO format
 */
void csrcoo (int n, int nz, const int *Ap, const  int *Aj, 
	     const Real *Ax, int *row, int *col, Real *data)
{
  int i, j, rowStart, rowStop, rowCount;
  // copy Aj and Ax to col and data
  for(i=0; i<nz; i++){
    col[i]  = Aj[i];
    data[i] = Ax[i];
  }

  // expand Ap into row
  rowCount = 0;
  for(i=0; i<n; i++){
    rowStart = Ap[i];
    rowStop  = Ap[i+1];
    for(j=rowStart; j<rowStop; j++){ 
      row[j] = rowCount;
    }
    rowCount++;
  }
}

/*!
 Converts COO to CSR format.
 */
csr_matrix * coocsr(int m, int n, int nz, const int *I, const int *J,
		    const Real *X)
{
  int i,k,*w;
  csr_matrix *A;
  
  // prepare workspace
  w = new_array<int> (m);

  // allocate C
  A = csr_alloc(m,n,nz,0);

  // copy values and columns to CSR
  for(i=0; i<nz; i++)
    {
      A->cols[i] = J[i];
      A->vals[i] = X[i];
      A->nnz++;
    }

  // row counts
  for(k=0; k<nz; k++)
  {
    i = I[k];
    w[i]++;    
  }

  // compute row pointers
  A->ptrs[0] = 0;
  for(i=0; i<m; i++)
  {
    A->ptrs[i+1] = A->ptrs[i] + w[i];
    w[i] = A->ptrs[i];
  }

  // free workspace memory
  delete_array<int> (w);

  return A;
}

/*!
 Converts CSR to DIA format
 @param[in]  n       dimension of A
 @param[in]  idiag   number of diagonals to be extracted
 @param[in]  a       CSR matrix array data
 @param[in]  ja      CSR matrix array column
 @param[in]  ia      CSR matrix array displacement
 @param[out] diag    real array of size (ndiag x idiag)
 */
void csrdia (int n, int idiag, const Real *a, int *ja, int *ia,
	     int ndiag, Real **diag, int *ioff)
{
  int i,j,k,l;
  int numdiag;
  int *ind;
  int n2 = n+n-1;

  ind = new_array<int> (n2+1);

  // extracts information about diagonals
  infdia (n, ja, ia, ind, &numdiag);

  // determine diagonals to extract
  int ii = 0;
  int jmax;
  
  while(ii < idiag){

    jmax = 0;

    for(k=0; k<=n2; k++){ // find the max of entries in the diagonals
      j = ind[k];
      if (j > jmax){
	i = k;
	jmax = j;
      }
    }

    if (jmax <= 0){
      ii = ii - 1;
      idiag = ii;
      break; // exit while loop
    }

    ioff[ii] = i - n;
    ind[i] = - jmax;
    ii++;
  }

  // initialize diag to zero
  for(i=0; i<n; i++){
    for(j=0; j<idiag; j++){
      diag[j][i] = 0.0;
    }
  }

  // extracts diagonals
  for(i=0; i<n; i++){
    for(k=ia[i]; k<ia[i+1]; k++){
      j = ja[k];
      for(l=0; l<idiag; l++){
	if(j-i == ioff[l]){
	  diag[l][i] = a[k]; // assumes that diag is (ndiag,n)
	}
      }
    }
  }
}

// -------------------------------------------------------------------
// CSR sparse format
// -------------------------------------------------------------------

/*!
 Count NNZ elements in a dense matrix K
 */
int csr_count_nnz (int m, int n, Real **K)
{
  int i,j,numNnz, total;
  total = 0;
  for(i=0; i<m; i++){
    numNnz = 0;
    for(j=0; j<n; j++){
      if(K[i][j] != 0.0) numNnz++;
    }
    total += numNnz;
  }
  return total;
}

/*!
 Create matrix in CSR format (cnt, col, ele) given a full matrix
 */
void csr_create_mat (int n, Real **A, Real *ele, int *dsp, int *col)
{
  int i, j, k, disp;
  k = 0;
  for (i=0; i<n; i++){
    disp = -1;
    for(j=0; j<n; j++){

      if(A[i][j] != 0.0) {
        ele[k] = A[i][j];
        col[k] = j;
        if(disp==-1) { dsp[i] = k; disp = 0; }
        k++;
      }
    }
  }
  dsp[i] = k;
}

/*!
  Extracts the inverse of the diagonal of the a CSR matrix
 */
void csr_inv_diag (int n, csr_matrix *A, Real *invA)
{
  int i, j;
  int row_start, row_stop;
  int *Ap, *Aj; 
  Real *Ax;
  Ap = A->ptrs;
  Aj = A->cols;
  Ax = A->vals;

#pragma omp parallel for private(i,j,row_start,row_stop)\
  shared(n,Ap,Aj,Ax,invA)
  for(i=0; i<n; i++){
    row_start = Ap[i];
    row_stop  = Ap[i+1]-1;
    for(j=row_start; j<=row_stop; j++){
      if(i == Aj[j]){
	invA[i] = 1./Ax[j];
	break;
      }
    }    
  }
}

/*!
 * 1 - csr
 * 0 - triplet
 */
csr_matrix * csr_alloc (const int m, const int n, const int nzmax, int triplet)
{
  // alloc matrix
  csr_matrix *A;
  A = (csr_matrix *) malloc(sizeof(csr_matrix));
  if (A == NULL)
  {
    fprintf(stderr,"csr_alloc: error allocating A\n");
    exit(1);
  }

  // init matrix
  A->m = m;
  A->n = n;
  A->nnz = 0;   // no nonzero entries so far
  A->nzmax = nzmax;

  if(!triplet){
    A->cmpr = 0; // compressed
    A->ptrs = new_array<int> (m+1);
    A->cols = new_array<int> (nzmax);
    A->vals = new_array<Real> (nzmax);
  } else {
    A->cmpr = -1; 
    A->T = triplet_create(nzmax);
  }
 
  return A;  
}

/**
 * Free allocated csr_matrix structure
 */
void csr_delete (csr_matrix *A)
{
  if(A->cmpr == -1){
    triplet_delete(A->T);
  } else { 
    delete_array<int> (A->ptrs);
    delete_array<int> (A->cols);
    delete_array<Real> (A->vals);
  }

  free(A);
}

/**
 *  Set the entry (i,j) of the sparse matrix (not compressed yet)
 */
int csr_entry_ij (csr_matrix *A, int i, int j, Real x)
{  
  if (A->cmpr == 0){
    fprintf(stdout,"csr_entry_ij: A is not a triplet.\n");
    exit(0);
  }

  int p, *Ti, *Tj;
  Real *Tx;

  Ti = A->T->i;
  Tj = A->T->j;
  Tx = A->T->x; 
  
  p = A->nnz;
  Ti[p] = i;
  Tj[p] = j;
  Tx[p] = x;
  A->nnz = A->nnz + 1;

  return 0; // no error
}

/**
 * A is in triplet form
 */
csr_matrix * csr_compress (csr_matrix *A)
{
  if (A->cmpr != -1) 
    error("csr_compress: the matrix is already compressed.\n");

  // local variables
  int i, k, p, m, n, nzmax, *Ti, *Tj, *Cj, *Cp, *w, serr;
  Real *Tx, *Cx;
  csr_matrix *C;
  
  m = A->m; n = A->n; nzmax = A->nzmax;
  Ti = A->T->i;
  Tj = A->T->j;
  Tx = A->T->x;

  // prepare workspace
  w = new_array<int> (m);

  // allocate C
  C = csr_alloc(m,n,nzmax,0);
  Cp = C->ptrs; 
  Cj = C->cols;
  Cx = C->vals;

  // row counts
  for(k=0; k<nzmax; k++)
  {
    i = Ti[k];
    w[i]++;    
  }

  // compute row pointers
  Cp[0] = 0;
  for(i=0; i<m; i++)
  {
    Cp[i+1] = Cp[i] + w[i];
    w[i] = Cp[i];
  }
 
  // construct the row form
  for(k=0; k<nzmax; k++)
  {
    p = w[Ti[k]]++;
    Cj[p] = Tj[k];
    Cx[p] = Tx[k];
    C->nnz = C->nnz + 1;
  }

  // sum up duplicates
  serr = csr_sumdup(C);
  if(serr) error("csr_compress: error after csr_cumdup.\n");
 
  // free workspace and temporary csr_matrix
  delete_array<int> (w);
  csr_delete(A);

  return C;
}

/**
 * Sum up duplicated entries
 */
int csr_sumdup (csr_matrix *A)
{
  int i, j, p, q, nz = 0 , n, m, *Ap, *Aj, *w ;
  Real *Ax ;
 
  m = A->m; n = A->n;
  Ap = A->ptrs;
  Aj = A->cols;
  Ax = A->vals;

  // prepare workspace
  w = new_array<int> (n);
  for(j=0; j<n; j++) w[j] = -1;

  for(i=0; i<m; i++)
  {
    q = nz;
    for(p=Ap[i]; p<Ap[i+1]; p++)
    {
      j = Aj[p];
      if(w[j] >= q){
        // A(i,j) is duplicate
        Ax[w[j]] += Ax[p];
      } else {
        // keep A(i,j)
        w[j] = nz; 
        Aj[nz] = j;
        Ax[nz++] = Ax[p];
      }
    }
    Ap[i] = q;    
  }
  Ap[n] = nz;

  A->nnz = nz;

  delete_array<int> (w);
 
  // remove extra space from A
  // need to change Aj and Ax only, Ap remains the same
  A->cols = (int *) realloc(A->cols, nz*sizeof(int));
  if(A->cols == NULL){
    fprintf(stdout,"csr_sumdup: error reallocating memory!\n");
    exit(1);
  }

  A->vals = (Real *) realloc(A->vals, nz*sizeof(Real));
  if(A->vals == NULL){
    fprintf(stdout,"csr_sumdup: error reallocating memory!\n");
    exit(1);
  }
  
  // no error
  return 0;
}

/**
 * Triplet format (aka COOrdinate format)
 */
triplet * triplet_create (int size)
{
  triplet *T;
  T = (triplet *) malloc(sizeof(triplet));

  T->size = size;
  T->i = new_array<int> (size);
  T->j = new_array<int> (size);
  T->x = new_array<Real> (size);

  return T;
}

void print_triplet (triplet *T, int nen, int nel)
{
  int a = 0;
  for(int k=0; k<nel; k++){
    for(int i=0; i<nen; i++){
      for(int j=0; j<nen; j++){
        Real aux = T->x[a];
        fprintf(stdout,"%lf ",aux);
        a++;
      }
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n\n");
  }
}

void triplet_delete (triplet *T)
{
  T-> size = 0;
  delete_array<int> (T->i);
  delete_array<int> (T->j);
  delete_array<Real> (T->x);
  free(T);
}


// -------------------------------------------------------------------
// sort_dia - Included by Ronan Amorim
// -------------------------------------------------------------------

/**
 * returns the offsets for the DIAgonal matrix in the right sequence
 * starting from the lowest diagonal
 */

int *
sort_dia(int *dia_off, int numdiags)
{
  int i,j;

  int *sort_idx = (int*) malloc(sizeof(int)*numdiags);
  int sorted_off[9];
  
  for(i=0; i<numdiags; i++){
    sorted_off[i] = dia_off[i];
    sort_idx[i] = i;
  }
  
  for(i=0; i<numdiags; i++){
    
    int min = sorted_off[i];
    int idx = i;
    
    for(j=i+1; j<numdiags; j++){
      
      if(sorted_off[j] < min){
	idx = j;
	min = sorted_off[j];
      }
    }
    
    sorted_off[idx] = sorted_off[i];
    sorted_off[i] = min;
    
    int swap_idx = sort_idx[i];
    sort_idx[i] = sort_idx[idx];
    sort_idx[idx] = swap_idx;
  }
  return sort_idx;
}

