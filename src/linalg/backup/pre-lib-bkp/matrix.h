#ifndef MATRIX_H
#define MATRIX_H

#include "array.hpp"
#include "vector.h"
#include "datatype.h"

using namespace std;

class Matrix
{
 private:
  uint nrows, ncols;
  Real **ets;
  
 public:
  Matrix(uint n, uint m, Real**);
  Matrix(uint n, uint m, Real d=0);
  Matrix(const Matrix &);
  ~Matrix();

  uint numrows() const {return nrows;}
  uint numcols() const {return ncols;}

  // Operators
  Matrix& operator=(const Matrix&);
  Matrix& operator+=(const Matrix&);
  Matrix& operator-=(const Matrix&);
  Vector<Real> operator*(const Vector<Real>&) const;
  Real* operator[](uint i) const { return ets[i]; }
  Real& operator()(uint i, uint j) { return ets[i][j]; }

  // Friends
  friend Matrix operator+(const Matrix&);
  friend Matrix operator-(const Matrix&);
  friend Matrix operator+(const Matrix&, const Matrix&);
  friend Matrix operator-(const Matrix&, const Matrix&);

};

#endif
