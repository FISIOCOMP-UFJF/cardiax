#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iomanip>
#include <iostream>
#include "array.hpp"
#include "util/util.hpp"

/*
class Matrix
{
protected:
  int m_, n_;
  Real * data;

  int index(int,int) const;
public:
  Matrix();
  Matrix(int m, int n);
  Matrix(const Matrix&);
  ~Matrix();

  int height() const { return(m_); }
  int width() const { return(n_); } 

  Real max();
  int size(int d) const;
  Matrix& resize(int,int);
  const Real & operator()(int, int) const;
  Real& operator()(int, int);
  Matrix& operator=(Real d);   
  Matrix& operator=(const Matrix& M2);
  Real* ref() const;
};

std::ostream& operator<<(std::ostream& o, const Matrix& M);
std::istream& operator>>(std::istream& i, Matrix& M);

// --------------------------------------------------------------
// IMPLEMENTATION
// --------------------------------------------------------------

inline int Matrix::index(int i,int j) const
{
#ifdef DEBUG
  if (i<0 || i>=size(0) || j<0 || j>= size(1))
    {
      std::cerr << "error in int Matrix::index(int i, int j) const: index out of bounds\n";
      exit(1);
    }
#endif
  return i*n_+j;
}

inline Matrix::Matrix() : m_(0), n_(0), data(NULL)
{}

inline Matrix::Matrix(int m, int n) : m_(m), n_(n)
{
  data = new Real[m*n];
}

inline Matrix::Matrix(const Matrix& M) : m_(M.m_), n_(M.n_)
{
  data = new Real[M.m_*M.n_];
  for(int i=0; i<m_*n_; i++)
    data[i] = M.data[i];
}

inline Matrix::~Matrix()
{
  //cout << "DESTRUCTORRRR "  << m_ << " " << n_ << " -> " << data <<  endl;
  delete [] data;
}

inline Real Matrix::max()
{
  Real m = data[0];
  for(int i=1; i<m_*n_; i++)
    if (data[i] > m) m = data[i];
  return m;
  //return data.max();
}

inline int Matrix::size(int d) const
{
  if (d==0)
    return m_;
  else if (d==1)
    return n_;
  else
    {
      std::cerr << " Error in int Matrix::size(int d) const:"
		<< " Matrix has only two dimension, (0,1)\n";
      exit(1);
    }
}

inline const Real & Matrix::operator()(int i,int j) const
{
  // fast acces ok, since checked inside matrix
  //return data[index(i,j)];
  return data[i*n_+j];
}

inline Real& Matrix::operator()(int i, int j)
{
  //return data[index(i,j)];
  return data[i*n_+j];
}

inline Matrix& Matrix::operator=(const Matrix& M2)
{
  resize(M2.m_, M2.n_);
  //data = M2.data; ////////////// AQUI
  m_ = M2.m_;
  n_ = M2.n_;
  data = new Real[m_*n_];
  for(int i=0; i<m_*n_; i++)
    data[i] = M2.data[i];
  return *this;
}

inline Matrix& Matrix::operator=(Real d)
{
  for(int i=0; i<m_*n_; i++) data[i] = d;
  //data = d;
  return *this;
}

inline Matrix& Matrix::resize(int m, int n)
{
  if (m < 0 || n < 0)
  {
    std::cerr << " Error in Matrix::resize(int m, int n):"
	      << " size negative!\n";
    exit(1);
  }

  if (m*n != m_*n_)
  {
    delete [] data;
    data = new Real[m*n];
    //data.resize(m*n);
  }

  m_ = m;
  n_ = n;
  return *this;
}

inline Real* Matrix::ref() const
{
  //return data.ref();
  return data;
}

// ------------------- input and output functions --------------

inline std::ostream& operator<<(std::ostream& o, const Matrix& M)
{
  for (int i=0; i<M.size(0); i++)
    {
      for (int j=0; j<M.size(1); j++)
	{
	  o << std::setprecision(8);
	  o << std::fixed;
	  o << (M(i,j)) << "\t";
	}
      o << "\n";
    }
  o << std::flush;
  return o;
}

inline std::istream& operator>>(std::istream& is, Matrix& M)
{
  // todo: add checks!!
  int m, n;
  is >> m >> n;
  M.resize(m,n);
  for (int i=0; i<m; i++)
    {
      for (int j=0; j<n; j++)
	is >> M(i,j);
    }
  return is;
}

//
// Some useful routines for elasticity problems
//

inline Real det_matrix_3x3(const Matrix & mat)
{
  if (mat.height() != 3 || mat.width() != 3)
    error("matrix.hpp", "calc det matrix 3x3", "not a 3x3 matrix");

  Real det = 0.0;
  det = mat(0,0)*mat(1,1)*mat(2,2)
      + mat(0,1)*mat(1,2)*mat(2,0)
      + mat(0,2)*mat(1,0)*mat(2,1)
      - mat(0,2)*mat(1,1)*mat(2,0)
      - mat(0,1)*mat(1,0)*mat(2,2)
      - mat(0,0)*mat(1,2)*mat(2,1);
  return det;
}

inline void inv_matrix_3x3(const Matrix & mat, Matrix & inv)
{
  if (mat.height() != 3 || mat.width() != 3)
    error("matrix.hpp", "calc det matrix 3x3", "not a 3x3 matrix");

  Real det = det_matrix_3x3(mat);

  if (det == 0.0) det = 1.0;
  det = 1.0/(det);

  inv(0,0) = (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) * (det);
  inv(0,1) = (mat(0,2)*mat(2,1) - mat(0,1)*mat(2,2)) * (det);
  inv(0,2) = (mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1)) * (det);

  inv(1,0) = (mat(1,2)*mat(2,0) - mat(1,0)*mat(2,2)) * (det);
  inv(1,1) = (mat(0,0)*mat(2,2) - mat(0,2)*mat(2,0)) * (det);
  inv(1,2) = (mat(0,2)*mat(1,0) - mat(0,0)*mat(1,2)) * (det);

  inv(2,0) = (mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0)) * (det);
  inv(2,1) = (mat(0,1)*mat(2,0) - mat(0,0)*mat(2,1)) * (det);
  inv(2,2) = (mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0)) * (det);
}

*/

#endif

