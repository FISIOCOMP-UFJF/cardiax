#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include <iostream>
#include <fstream>
#include "vector.h"
#include "../util/datatype.h"
#include "../util/malloc.h"

//
// Sparse matrix abstract base class
//
class SparseMatrix
{
 protected:
  int nrows;
  int ncols;
  int nnz;
  
 public:
  SparseMatrix(){};
  SparseMatrix(int m, int n, int nz=0) : nrows(m), ncols(n), nnz(nz) {}
  ~SparseMatrix(){};

  int rows() const {return nrows;};
  int cols() const {return ncols;};
  int nz() const {return nnz;};

  void set_nz(int nz) {nnz = nz;}

  // Sparse matrix vector multiplitcation interface
  virtual Vector<Real> operator*(const Vector<Real> &) const = 0;
  virtual void spmv(const Real *x, Real *y) = 0;
  virtual void info() = 0;

};

//
// COO/Triplet matrix
//
class Triplet : public SparseMatrix
{
 public:
  int * Ai;
  int * Aj;
  int cont_nz;
  Real * Ax;

 public:
  Triplet(int m, int n, int s);
  Triplet(int m, int n, int s, int *ptr_i, int *ptr_j, Real *ptr_val);
  ~Triplet();

  Vector<Real> operator*(const Vector<Real> &) const {};
  void spmv(const Real *x, Real *y){};
  void info();
  void set_entry_ij(int i, int j, Real x);
};


//
// CSR matrix (compressed sparse row matrix)
//
class CompressedMatrix : public SparseMatrix
{
 private:
  int nzmax;
  int *cols;
  int *ptrs;
  Real *vals;

  bool is_compressed;
  void sum_duplicates();

 public:
  CompressedMatrix();
  CompressedMatrix(int m, int n);
  CompressedMatrix(int m, int n, int nz, Real *Ax, int *Aj, int *Ap);
  CompressedMatrix(const Triplet & tri);
  CompressedMatrix(char const * mm_filename);
  ~CompressedMatrix(){};

  Vector<Real> operator*(const Vector<Real> &) const;
  void spmv(const Real *x, Real *y);
  void info();
  void compress();
  void convert(const Triplet & tri);

  // Operators
  const Real & operator()(int, int) const;
  Real & operator()(int, int); 

  int * get_cols() const { return cols; }
  int * get_ptrs() const { return ptrs; }
  Real * get_vals() const { return vals; }

  // Output
  std::ostream& print (std::ostream & o) const;
};

#endif
