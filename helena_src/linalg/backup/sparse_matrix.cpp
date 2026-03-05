#include <iostream>
#include <string>

#include "sparse_matrix.h"

extern "C"{
#include "../util/mmio.h"
}

using namespace std;

//
// Triplet format
//
Triplet::Triplet(int m, int n, int s) : SparseMatrix(), cont_nz(0)
{
  nrows = m;
  ncols = n;
  nnz = s;
  Ai = new_array<int>(s);
  Aj = new_array<int>(s);
  Ax = new_array<Real>(s);
}

Triplet::Triplet(int m, int n, int s, int *ptr_i, int *ptr_j, 
		 Real *ptr_val)
{
  nrows = m;
  ncols = n;
  nnz = s;
  Ai = ptr_i;
  Aj = ptr_j;
  Ax = ptr_val;
}

Triplet::~Triplet()
{
  delete_array<Real>(Ax);
  delete_array<int>(Ai);
  delete_array<int>(Aj);
}

void Triplet::info()
{
  cout << "Triplet matrix information" << endl;
  cout << " Number of rows:     " << nrows << endl;
  cout << " Number of cols:     " << ncols << endl;
  cout << " Number of non-zero: " << nnz << endl;
}

void Triplet::set_entry_ij(int i, int j, Real x)
{
  int p = cont_nz;
  Ai[p] = i;
  Aj[p] = j;
  Ax[p] = x;

  cont_nz++;  
}

//
// Compressed sparse row format member functions
//
CompressedMatrix::CompressedMatrix() 
  : SparseMatrix(), vals(0), ptrs(0), cols(0)
{
}

CompressedMatrix::CompressedMatrix(int m, int n) :  
  SparseMatrix()
{
  nrows = m;
  ncols = n;
  nnz = m*n;
  vals = NULL;
  ptrs = NULL;
  cols = NULL;
}

CompressedMatrix::CompressedMatrix(int m, int n, int nz, Real *Ax, 
				   int *Aj, int *Ap) :
  SparseMatrix(m,n,nz)
{
  vals = Ax;
  cols = Aj;
  ptrs = Ap;
}

void CompressedMatrix::info()
{
  cout << "Compressed sparse row matrix information" << endl;
  cout << " Number of rows:     " << nrows << endl;
  cout << " Number of cols:     " << ncols << endl;
  cout << " Number of non-zero: " << nnz << endl;
}

void CompressedMatrix::convert(const Triplet & tri)
{
  // local variables
  int nzmax;
  int *Ti, *Tj, *Cj, *Cp, *w;
  Real *Tx, *Cx;

  nrows = tri.rows();
  ncols = tri.cols();
  nzmax = tri.nz();

  Ti = tri.Ai; 
  Tj = tri.Aj; 
  Tx = tri.Ax;

  // prepare workspace
  w = new_array<int>(nrows);

  // allocate CSR arrays
  Cp = new_array<int>(nrows+1);
  Cj = new_array<int>(nzmax);
  Cx = new_array<Real>(nzmax);

  // rows counts
  for(int k=0; k<nzmax; k++)
  {
    int i = Ti[k];
    w[i]++;
  }

  // compute row pointers
  Cp[0] = 0;
  for(int i=0; i<nrows; i++)
  {
    Cp[i+1] = Cp[i] + w[i];
    w[i] = Cp[i];
  }

  // construct the row form
  for(int k=0; k<nzmax; k++)
  {
    int p = w[Ti[k]]++;
    Cj[p] = Tj[k];
    Cx[p] = Tx[k];
    nnz = nnz + 1;
  }

  ptrs = Cp;
  cols = Cj;
  vals = Cx;

  sum_duplicates();  
}

void CompressedMatrix::sum_duplicates()
{
  int i, j, p, q, nz = 0, *Ap, *Aj, *w;
  Real *Ax;
 
  Ap = ptrs;
  Aj = cols;
  Ax = vals;

  // prepare workspace
  w = new_array<int>(ncols);
  for(j=0; j<ncols; j++) w[j] = -1;
  
  for(i=0; i<nrows; i++)
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
  Ap[ncols] = nz;

  nnz = nz;
  
  delete_array(w);
  
  // remove extra space from A
  // need to change Aj and Ax only, Ap remains the same
  Aj = resize_array<int>(Aj,nz);
  Ax = resize_array<Real>(Ax,nz);
}

const Real & CompressedMatrix::operator()(int i, int j) const
{
  for(int k=ptrs[i]; k<ptrs[i+1]; k++)
  {
    if(cols[k] == j)
      return vals[k];
  }
  
  // return the last non-zero entry
  return vals[nzmax];
}

Real & CompressedMatrix::operator()(int i, int j)
{
  for(int k=ptrs[i]; k<ptrs[i+1]; k++)
  {
    if(cols[k] == j)
      return vals[k];
  }
  
  // return the last non-zero entry
  return vals[nzmax];
}

Vector<Real> CompressedMatrix::operator*(const Vector<Real> & x) const
{
  int n = nrows;
  if(x.size() != n)
    error(" matrix vector sizes do not match");

  Vector<Real> y(n);
  int i, j, k1, k2, q;
  Real sum;

#pragma omp parallel for private(i,j,k1,k2,q,sum)	\
                         shared(n,cols,ptrs,vals,x,y)
  for(i=0; i<n; i++)
  {
    k1 = ptrs[i];
    k2 = ptrs[i+1];
    sum = 0.0;
    for(j=k1; j<k2; j++)
    {
      q = cols[j];
      sum += vals[j] * x[q];
    }
    y[i] = sum;
  }
  return y;
}

void CompressedMatrix::spmv(const Real *x, Real *y)
{
  int n = nrows;
  int i, j, k1, k2, q;
  Real sum;

#pragma omp parallel for private(i,j,k1,k2,q,sum)	\
                         shared(n,cols,ptrs,vals,x,y)
  for(i=0; i<n; i++)
  {
    k1 = ptrs[i];
    k2 = ptrs[i+1];
    sum = 0.0;
    for(j=k1; j<k2; j++)
    {
      q = cols[j];
      sum += vals[j] * x[q];
    }
    y[i] = sum;
  }
}

CompressedMatrix::CompressedMatrix(char const * mm_filename) : 
  SparseMatrix()
{
  /*
  int M, N, nz, ret_code;
  int i, *I, *J;
  real *val;
  FILE * fid;
  MM_typecode matcode;
  
  // read from matrix market format
  fid = fopen(mm_filename, "r");
  
  if (fid == NULL){
    printf("Unable to open file %s\n", mm_filename);
    exit(1);
  }
  
  if (mm_read_banner(fid, &matcode) != 0){
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }
  
  if (!mm_is_valid(matcode)){
    printf("Invalid Matrix Market file.\n");
    exit(1);
  }
  
  if (!((mm_is_real(matcode) || mm_is_integer(matcode) || 
	 mm_is_pattern(matcode)) && mm_is_coordinate(matcode) && 
	mm_is_sparse(matcode) ) ){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    printf("Only sparse real-valued or pattern coordinate matrices are supported\n");
    exit(1);
  }

  // find out sizes
  if ((ret_code = mm_read_mtx_crd_size(fid, &M, &N, &nz)) !=0)
    exit(1);

  nrows = M;
  ncols = N;
  nnz = nz;
  nzmax = nz;

  // allocate arrays and setup
  isCompressed = false;

  // allocate triplet
  I = new_array<int>(nz);
  J = new_array<int>(nz);
  val = new_array<real>(nz);

  // start reading
  for (i=0; i<nz; i++)
  {
    fscanf(fid, "%d %d %f\n", &I[i], &J[i], &val[i]);
    I[i]--; 
    J[i]--; 
  } 

  // build temporary triplet structure
  t = new Triplet(M,N,nz,I,J,val);

  // now compress to CSR format
  compress();
  */
}


