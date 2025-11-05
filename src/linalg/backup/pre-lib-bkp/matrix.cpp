#include "matrix.h"

inline void error(string v)
{
  cout << v << ". program exited\n";
  exit(1);
}

Matrix::Matrix(uint n, uint m, Real ** dbp)
{
  nrows = n;
  ncols = m;
  ets = new Real*[nrows];
  for(uint i=0; i<nrows; i++)
  {
    ets[i] = new Real[ncols];
    for(uint j=0; j<ncols; j++) ets[i][j] = dbp[i][j];
  }    
}

Matrix::Matrix(uint n, uint m, Real a)
{
  nrows = n;
  ncols = m;
  ets = new Real*[nrows];
  for(uint i=0; i<nrows; i++)
  {
    ets[i] = new Real[ncols];
    for(uint j=0; j<ncols; j++) ets[i][j] = a;
  }    
}

Matrix::Matrix(const Matrix& mat)
{
  nrows = mat.nrows;
  ncols = mat.ncols;
  ets = new double* [nrows];
  for(uint i=0; i<nrows; i++)
  {
    ets[i] = new Real[ncols];
    for(uint j=0; j<ncols; j++) ets[i][j] = mat[i][j];
  }    
}

inline Matrix::~Matrix()
{
  for(uint i=0; i<nrows; i++) delete [] ets[i];
  delete [] ets;
}

// Operators

Matrix& Matrix::operator=(const Matrix& mat)
{
  if(this != &mat)
  {
    if(nrows != mat.nrows || ncols != mat.ncols)
      error(" bad matrix sizes");
    for(uint i=0; i<nrows; i++)
      for(uint j=0; j<ncols; j++)
	ets[i][j] = mat[i][j];
  }
  return *this;
}

Matrix& Matrix::operator+=(const Matrix& mat)
{
  if(nrows != mat.nrows || ncols != mat.ncols)
    error(" bad matrix sizes");
  for(uint i=0; i<nrows; i++)
    for(uint j=0; j<ncols; j++)
      ets[i][j] += mat[i][j];
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& mat)
{
  if(nrows != mat.nrows || ncols != mat.ncols)
    error(" bad matrix sizes");
  for(uint i=0; i<nrows; i++)
    for(uint j=0; j<ncols; j++)
      ets[i][j] -= mat[i][j];
  return *this;
}

inline Matrix operator+(const Matrix& mat) 
{
  return mat;
}

inline Matrix operator-(const Matrix& mat)
{
  return Matrix(mat.nrows, mat.ncols) - mat;
}

// m = m1 + m2
Matrix operator+(const Matrix& m1, const Matrix& m2)
{
  if(m1.nrows != m2.nrows || m1.ncols != m2.ncols)
    error(" bad matrix sizes");
  Matrix sum = m1;
  sum += m2;
  return sum;
}

// m = m1 - m2
Matrix operator-(const Matrix& m1, const Matrix& m2)
{
  if(m1.nrows != m2.nrows || m1.ncols != m2.ncols)
    error(" bad matrix sizes");
  Matrix sum = m1;
  sum -= m2;
  return sum;
}

// u = M*v
Vector Matrix::operator*(const Vector& v) const
{
  if(ncols != v.size())
    error(" matrix and vector sizes do not match");
  Vector tmp(nrows);
  for(uint i=0; i<nrows; i++)
    for(uint j=0; j<ncols; j++)
      tmp[i] += ets[i][j]*v[j];
  return tmp;
}


