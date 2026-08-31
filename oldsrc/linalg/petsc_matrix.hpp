#ifndef PETSC_MATRIX_HPP_
#define PETSC_MATRIX_HPP_

#include <cassert>
#include "petscksp.h"
#include "petsc_vector.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace petsc
{

class Matrix
{
public:

  Matrix();
  Matrix(const Matrix & m);
  ~Matrix();

  void assemble ();
  void add (const int i, const int j, double s);
  void add (const int m, const int n, const int * idxm, const int * idxn,
	    const double * v);
  void create (int r, int c, int nz);
  void destroy ();
  uint get_nnz();
  void get_CSR(int *n, int *idxrow, int *idxcol, double *v);
  void store_CSR();
  bool is_null() const;
  bool is_sym() const;
  Mat mat() const;
  void mult (const Vector & x, Vector & y);
  void set (const int i, const int j, double s);
  void set (const int m, const int n, const int * idxm, const int * idxn,
	    const double * v);
  uint size();
  void set_symmetric();  
  void view();
  void setNullSpace(petsc::Vector *coord);
  

  /*!zero_rows_cols(int num_rows, int * rows, double diag)
   * Wrapper for PETSc MatZeroRowsColumns
   * Requires PETSc Version >= 3.2
   */
  void zero_rows_cols(int num_rows, int * rows, double diag);
  void zero_rows_cols(int num_rows, int * rows, double diag, 
		      petsc::Vector & x, petsc::Vector & b);
  Matrix & operator = (const double s);

private:

  Mat _mat;
  PetscErrorCode ierr;

};

// ------------------- inline and template functions --------------

inline Matrix::Matrix() : _mat(0)
{
  // do nothing
}

inline Matrix::~Matrix()
{
  destroy();
}

inline void Matrix::destroy ()
{
  if (_mat != NULL)
  {
    ierr = MatDestroy(&_mat);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
}

inline Mat Matrix::mat() const
{
  return _mat;
}

inline bool Matrix::is_null() const
{
  if (_mat == NULL)
    return true;
  else
    return false;
}

inline bool Matrix::is_sym() const
{
  PetscBool sym;
  MatIsSymmetric(_mat, 1.0e-8, &sym);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  return sym;
}

}

#endif /* PETSC_MATRIX_HPP_ */
