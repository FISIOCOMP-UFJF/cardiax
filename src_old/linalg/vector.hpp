#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <iostream>
#include "matrix.hpp"

/*
template<class Left, class Op, class Right>
struct VectorExpression {
  const Left & left;
  const Right & right;
  VectorExpression(const Left & l, const Right & r) : left(l), right(r) { }

  inline Real operator[](int i) const {
    return Op::apply(left[i], right[i]);
  }
};

template<class Vexp>
struct VectorScalarExpression {
  const Vexp & a;
  Real s;
  VectorScalarExpression(const Vexp & v, Real scalar) : a(v), s(scalar) {}
  
  inline Real operator[](int i) const {
    return s * a[i];
  }
};

class Vector : public Matrix
{
public:

  Vector();
  Vector(int n);
  Vector(const Vector&);
  Vector(const Matrix&); // with size check!!
  ~Vector();
		
  int size() const;   

  Real operator()(int) const;        // element access:   d=v(n)
  Real& operator()(int);             // element access:   v(n)=d
  Real operator[](int) const;        // element access:   d=v(n)
  Real& operator[](int);             // element access:   v(n)=d
  Vector& operator=(const Vector& w);  // assignment copy 
  Vector& operator=(const Matrix& M);  // include size check
  Vector& operator=(Real d);

  //Vector& operator+=(const Vector &v);
  // friend Vector operator+(const Vector &);
  //friend Vector operator*(const real scalar, const Vector &v);

  // Resize internal array for data
  Vector& resize(int);
		
  // Return reference to storage    
  Real* ref() const; 

  // Expression template for vectorial expressions
  template<class Left, class Op, class Right>
  void operator=(VectorExpression<Left,Op,Right> expr) {
    for(int i=0; i<size(); i++) data[i] = expr[i];
  }

  template<class Left, class Op, class Right>
  void operator+=(VectorExpression<Left,Op,Right> expr) {
    for(int i=0; i<size(); i++) data[i] += expr[i];
  }

  template<class Left, class Op, class Right>
  void operator-=(VectorExpression<Left,Op,Right> expr) {
    for(int i=0; i<size(); i++) data[i] -= expr[i];
  }

  // Expression template for vector-scalar operations
  template<class Vexp>
  void operator=(VectorScalarExpression<Vexp> expr) {
    for(int i=0; i<size(); i++) data[i] = expr[i];
  }

  template<class Vexp>
  void operator+=(VectorScalarExpression<Vexp> expr) {
    for(int i=0; i<size(); i++) data[i] += expr[i];
  }

  template<class Vexp>
  void operator-=(VectorScalarExpression<Vexp> expr) {
    for(int i=0; i<size(); i++) data[i] -= expr[i];
  }


};

// Vector expressions

struct Add {
  static Real apply(Real a, Real b){
    return a+b;
  }
};

template<class Left> 
VectorExpression<Left, Add, Vector> operator+(Left a, Vector b) {
  return VectorExpression<Left, Add, Vector>(a,b);
}

struct Minus {
  static Real apply(Real a, Real b){
    return a-b;
  }
};

template<class Left> 
VectorExpression<Left, Minus, Vector> operator-(Left a, Vector b) {
  return VectorExpression<Left, Minus, Vector>(a,b);
}

// Vector-scalar expressions
template<class L> 
VectorScalarExpression<L> operator*(const L & v, Real a) {
  return VectorScalarExpression<L>(v,a);
}

template<class R> 
VectorScalarExpression<R> operator*(Real a, const R & v) {
  return VectorScalarExpression<R>(v,a);
} 

// output funktion: inherited from matrix
//std::ostream& operator<<(std::ostream& o, const Vector& v);

// input oprator
std::istream& operator>>(std::istream& is, Vector& v);


// --------------------------------------------------
// IMPLEMENTIERUNG 
// --------------------------------------------------
// Bem: momentan werden alle funktionen inline definiert!!

// Vector v(5);
inline Vector::Vector() : Matrix()
{}

// Vector v(5);
inline Vector::Vector(int n) : Matrix(n,1)
{}

// copy-Konstruktor: wichtig!
inline Vector::Vector(const Vector& v) : Matrix(v)
{}

// copy-Konstruktor: wichtig!
inline Vector::Vector(const Matrix& M) : Matrix(M)
{
  if (M.size(1) != 1)
  {
    std::cerr << "*** error in Vector::Vector(const Matrix& M) : Matrix(M): M must have exactly 1 column" << std::endl;
    exit(0);
  }
}

inline Vector::~Vector()
{}


inline int Vector::size() const
{
  return m_;
}

// r=v(i)
inline Real Vector::operator()(int i) const
{
#ifdef DEBUG  
  if (i<0 || i>=m_)
  {
    std::cerr << "*** error in real Vector::operator()(int i) const: index " << i << " out of range" << std::endl;
    exit(0);
  }
#endif
  return data[i];
}

// v(i)=r 
inline Real& Vector::operator()(int i)
{
#ifdef DEBUG  
  if (i<0 || i>=m_)
{
  std::cerr << "*** error in real& Vector::operator()(int i): index out of range" << std::endl;
  exit(0);
}
#endif
  return data[i];
}

// r=v(i)
inline Real Vector::operator[](int i) const
{
#ifdef DEBUG  
  if (i<0 || i>=m_)
  {
    std::cerr << "*** error in real Vector::operator()(int i) const: index " << i << " out of range" << std::endl;
    exit(0);
  }
#endif
  return data[i];
}

// v(i)=r 
inline Real& Vector::operator[](int i)
{
#ifdef DEBUG  
  if (i<0 || i>=m_)
{
  std::cerr << "*** error in real& Vector::operator()(int i): index out of range" << std::endl;
  exit(0);
}
#endif
  return data[i];
}
		

// Operators
		
inline Vector& Vector::operator=(const Vector& w)
{
  Matrix::operator=(w);
  return *this;
}

inline Vector& Vector::operator=(const Matrix& M)
{
  if (M.size(1) != 1)
  {
    std::cerr << "*** error in Vector& Vector::operator=(const Matrix& M): M must have exactly one column" << std::endl;
    exit(0);
  }
  Matrix::operator=(M);
  return *this;
}

// v=d; all elements to d
inline Vector& Vector::operator=(Real d)
{
  Matrix::operator=(d);
  return *this;
}


// inline Vector& Vector::operator+=(const Vector & v)
// {
//   if (size() != v.size()) std::cerr << "bad vec sizes" << endl;
//   for (int i=0; i<size(); i++) data(i) += v(i);
//   return *this;
// }


// Vector operator*(const real scalar, const Vector & v)
// {
//   Vector tmp(v.size());
//   for (int i=0; i<v.size(); i++) tmp(i) = scalar * v(i);
//   return tmp;
// }

// Vector operator+(const Vector & v1, const Vector & v2)
// {
//   if(v1.size() != v2.size()) std::cerr << "erro in op+" << endl;
//   Vector sum = v1;
//   for (int i=0; i<v1.size(); i++) sum(i) += v2(i);
//   return sum;
// }


inline Vector& Vector::resize(int n)
{
  Matrix::resize(n,1);
  return *this;
}

inline Real* Vector::ref() const
{
  return Matrix::ref();
}

// output: cout << v << endl;
inline std::ostream& operator<<(std::ostream& o, const Vector& v)
{
  o << std::setprecision(8);
  o << std::fixed;
  for (int i=0; i<v.size(); i++)
    o << v(i) << endl;
  return o;
}

// input
inline std::istream& operator>>(std::istream& is, Vector& v)
{
  int n;
  is >> n;
  if (n<0)
  {
    std::cerr << "Error in std::istream& operator>>(std::istream& is, Vector& v):"
	      << "size must be positive" << std::endl;
    exit(0);
  }
    
  v.resize(n);
  for (int i=0; i<n; i++)
    is >> v(i);
	
  return is;
}
*/
	
#endif
