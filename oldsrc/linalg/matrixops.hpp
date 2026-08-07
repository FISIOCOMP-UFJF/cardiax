#ifndef MATOPS_HPP
#define MATOPS_HPP

#include <cmath>
#include <assert.h>

#include "matrix.hpp"
#include "vector.hpp" 
#include "util/datatype.hpp"

/*
// Matrix operations
Matrix operator-(const Matrix& A);

// Matrix real operations --------------------------
Matrix& operator+=(Matrix& A, Real a);
Matrix& operator-=(Matrix& A, Real a);
Matrix& operator*=(Matrix& A, Real a); 
Matrix& operator/=(Matrix& A, Real a);

Matrix operator+(Real a, const Matrix& A); 
Matrix operator+(const Matrix& A, Real a); 
Matrix operator-(Real a, const Matrix& A); 
Matrix operator-(const Matrix& A, Real a); 
Matrix operator*(const Real a, const Matrix& A); 
Matrix operator*(const Matrix& A, const Real a); 
Matrix operator/(const Matrix& A, Real a);


// Matrix Matrix operations ---------------------------
Matrix& operator+=(Matrix& A, const Matrix& B);
Matrix& operator-=(Matrix& A, const Matrix& B);

// Matrix operator/(const Matrix& A, const Matrix& B); // inv(A)*B

// Matrix functions ----------------------------------
Matrix transpose(const Matrix& M); 

// Matrix inverse(const Matrix& M); 
Matrix eye(int n);

// inverse of a Matrix
// Matrix inv(const Matrix& M);

// pointwise functions, sin, cos, exp, .. comming soon
Matrix abs(const Matrix& M);
Real max(const Matrix& M);

Real sprod(const Vector& v, const Vector & w);
Vector cross(const Vector & a, const Vector & b);
Real l2Norm(const Vector & v);   

//inline Matrix mult(const Matrix& A, const Matrix& B) { return A*B; }

// -------------------------------------------------------------------
// IMPLEMENTATION
// -------------------------------------------------------------------

// Matrix functions
inline Matrix operator-(const Matrix& A)
{
  Matrix R(A);
  Real *r=R.ref();
  for (int i=0; i<R.size(0)*R.size(1); i++)
    r[i] = -r[i];
  return R;
}  

// Matrix real ;operations --------------------------
inline Matrix& operator+=(Matrix& A, Real a)
{
  Real *r=A.ref();
  for (int i=0; i<A.size(0)*A.size(1); i++)
    r[i] += a;
  return A;
}

inline Matrix& operator-=(Matrix& A, Real a)
{
  return A+=(-a);
}

inline Matrix& operator*=(Matrix& A, Real a)
{
  Real *r=A.ref();
  for (int i=0; i<A.size(0)*A.size(1); i++){
    r[i] *= a;
    //if (r[i] < 1.0e-10) {
    //  cout << "OPPAAAAAAAAAAAAAA" << endl;
    //  cout << r[i] << endl;
    //}
  }
  return A;
}

inline Matrix& operator/=(Matrix& A, Real a)
{
  return A*=(1.0/a);
}



inline Matrix operator+(Real a, const Matrix& A)
{
  Matrix R(A);
  R += a;
  return R;
}

inline Matrix operator+(const Matrix& A, Real a)
{
  Matrix R(A);
  R += a;
  return R;
}

inline Matrix operator-(Real a, const Matrix& A)
{
  Matrix R(A);
  R = a;
  R -= A;
  return R;
}

inline Matrix operator-(const Matrix& A, Real a)
{
  Matrix R(A);
  R -= a;
  return R;
}

inline Matrix operator*(const Real a, const Matrix& A)
{
  Matrix R(A);
  R *= a;
  return R;
}

inline Matrix operator*(const Matrix& A, const Real a)
{
  Matrix R(A);
  R *= a;
  return R;

}

inline Matrix operator/(const Matrix& A, Real a)
{
  Matrix R(A);
  R *= (1.0/a);
  return R;
}

// Matrix Matrix operations ---------------------------

inline Matrix& operator+=(Matrix& A, const Matrix& B)
{
  if (A.size(0) != B.size(0) || A.size(1) != B.size(1))
  {
    std::cerr << " Error in Matrix& opertor+=(Matrix& A, const Matrix& B):"
	      << " dimensions do not match" << std::endl;
    exit(0);
  }
  Real *a=A.ref();
  Real *b=B.ref();
  for (int i=0; i<A.size(0)*A.size(1); i++)
    a[i] += b[i];
  return A;
}

inline Matrix& operator-=(Matrix& A, const Matrix& B)
{
  if (A.size(0) != B.size(0) || A.size(1) != B.size(1))
  {
    std::cerr << " Error in Matrix& opertor+=(Matrix& A, const Matrix& B):"
	      << " dimensions do not match" << std::endl;
    exit(0);
  }
  Real *a=A.ref();
  Real *b=B.ref();
  for (int i=0; i<A.size(0)*A.size(1); i++)
    a[i] -= b[i];
  return A;
}


inline Matrix operator+(const Matrix& A, const Matrix& B)
{
  Matrix R(A);
  R += B;
  return R;
}

inline Matrix operator-(const Matrix& A, const Matrix& B)
{
  Matrix R(A);
  R -= B;
  return R;
}

inline Matrix operator*(const Matrix& A, const Matrix& B)
{
  if (A.size(1) != B.size(0))
  {
    std::cerr << " Error in Matrix operator*(const Matrix& A, const Matrix& B):"
	      << " dimensions do not match" << std::endl;
    exit(0);
  }

  Matrix R(A.size(0),B.size(1));
  for (int i=0; i<A.size(0); i++)
  {
    for (int j=0; j<B.size(1); j++)
    {
      Real r=0.0;
      for (int k=0; k<A.size(1); k++)
        r += A(i,k)*B(k,j);
      R(i,j) = r;
      //cout << " val=" << r << endl;
    }
  }
  return R;
}

// Matrix functions ----------------------------------

inline Matrix transpose(const Matrix& M)
{
  Matrix R(M.size(1),M.size(0));
  
  for (int i=0; i<M.size(0); i++)
    for (int j=0; j<M.size(1); j++)
      R(j,i) = M(i,j);
  return R;
}

inline Matrix eye(int n)
{
  Matrix R(n,n);
  for (int i=0; i<n; i++)
    R(i,i) = 1;
  return R;
}
  
// pointwise functions, sin, cos, exp, .. comming soon
inline Matrix abs(const Matrix& M)
{
  Matrix R(M);
  for (int i=0; i<M.size(0); i++)
  {
    for (int j=0; j<M.size(1); j++)
    {
      R(i,j) = abs( M(i,j));
    }
  }
  return R;
}

inline Real max(const Matrix& M)
{
  Real r=M(0,0);
  for (int i=0; i<M.size(0); i++)
  {
    for (int j=0; j<M.size(1); j++)
    {
      r = std::max(r, M(i,j));
    }
  }
  return r;
}

inline Real sprod(const Vector& v, const Vector& w)
{
  assert(v.size()==w.size());
  Real sp = 0; 
  for (int k=0;k<v.size();k++)
      sp += v(k)*w(k);
       
  return(sp);
}   

inline Real l2Norm(const Vector & v)
{
  Real sp = sprod(v,v);
  return(sqrt(sp)); 
}   


inline Vector cross(const Vector & a, const Vector & b)
{
	Vector c(3);
	c(0) = a(1)*b(2) - a(2)*b(1);
	c(1) = a(2)*b(0) - a(0)*b(2);
	c(2) = a(0)*b(1) - a(1)*b(0);
	return c;
}
*/


/*inline Matrix inverse(const Matrix& M)
{
  return gauss(M,eye(M.size(1));
}*/

// inline Matrix operator/(const Matrix& A, const Matrix& B) // inv(A)*B
// {
//   return gauss(A,B);
// }      

/*
#include "gauss.hpp"
             
Matrix inv(const Matrix& M)
{
  Matrix I=eye(M.size(0));
  cout << "I=" << I;
  return gauss(M,I);
}             
*/

#endif


