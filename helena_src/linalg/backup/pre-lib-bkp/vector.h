#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <string>
#include <algorithm>
#include <math.h>
#include "array.hpp"
#include "saxpy.h"
#include "exptempl.h"
#include "datatype.h"

using namespace std;

// Forward declarations
template<class T> struct Sax;
template<class T> struct Saxpy;

/// Vector class definition using Expression templates
template<class T> class Vector
{
 private:
  uint _n;
  T * _vec;

 public:
  Vector(uint, const T *);
  Vector(uint n=0, T d=0);
  Vector(const Vector<T>&);
  Vector(const Saxpy<T>&); // constructor from Saxpy
  ~Vector(){ delete[] _vec; }

  uint size() const;

  T operator()(uint) const; // d = v(n)  
  T& operator()(uint);      // v(n) = d
  T operator[](uint) const; // d = v(n)  
  T& operator[](uint);      // v(n) = d

  Vector& operator=(const T);
  Vector& operator=(const Vector &);
  Vector& operator=(const Saxpy<T> &);
  Vector& operator+=(const Vector &);
  Vector& operator-=(const Vector &);

  Vector& operator+=(const Sax<T> &);
  Vector& operator-=(const Sax<T> &);
  Vector& operator*=(const T);

  T maxnorm() const;
  T twonorm() const;

  T *ref() const; // return reference to internal storage array       

  // Friends
  template<class S>
  friend Vector operator+(const Vector&); // v = +v2

  template<class S>
  friend Vector operator-(const Vector&); // v = -v2

  template<class S>
  friend Vector operator+(const Vector&, const Vector&); // v = v1+v2

  template<class S>
  friend Vector operator-(const Vector&, const Vector&); // v = v1-v2

  //friend Vector operator*(T, const Vector&); // vec-scalar multiply

  template<class S>
  friend Vector operator/(const Vector&, T); // vec-scalar divide

  template<class S>
  friend Vector operator*(const Vector&, const Vector&); // vec multiply

  template<class S>
  friend T dot(const Vector&, const Vector&); // dot (inner) product

  template<class S>
  friend ostream& operator<<(ostream&, const Vector&); // output operator

};

// Definitions

inline void error(string v)
{
  cout << v << ". program exited\n";
  exit(1);
}

template<class T>
Vector<T>::Vector(uint n, const T * w)
{
  _n = n;
  _vec = new T[n];
  for(uint i=0; i<n; i++) 
    _vec[i] = *(w+i);
}

template<class T>
Vector<T>::Vector(uint n, T d)
{
  _n = n;
  _vec = new T[n];
  for(uint i=0; i<n; i++)
    _vec[i] = d;
}

template<class T>
Vector<T>::Vector(const Vector & v)
{
  _n = v.size();
  _vec = new T[_n];
  for(uint i=0; i<v.size(); i++)
    _vec[i] = v[i];
}

template<class T>
Vector<T>::Vector(const Saxpy<T>& sp)
{
  _n = sp.x.size();
  _vec = new T[_n];
  for(uint i=0; i<_n; i++)
    _vec[i] = sp.a*sp.x[i] + sp.y[i];
}

template<class T>
inline uint Vector<T>::size() const
{
  return _n;
}


template<class T>
inline T Vector<T>::operator()(uint i) const
{
#ifdef DEBUG  
  if (i<0 || i>=_n)
  {
    cerr << "*** error in Vector<T>::operator()(uint i) const: index " 
	 << i << " out of range" << endl;
    exit(0);
  }
#endif
  return _vec[i];
}

template<class T>
inline T& Vector<T>::operator()(uint i)
{
#ifdef DEBUG  
  if (i<0 || i>=_n)
{
  cerr << "*** error in T& Vector<T>::operator()(uint i): index out of range" << endl;
  exit(0);
}
#endif
  return _vec[i];
}

template<class T>
inline T Vector<T>::operator[](uint i) const
{
  return _vec[i];
}

template<class T>
inline T& Vector<T>::operator[](uint i)
{
  return _vec[i];
}

template<class T>
Vector<T>& Vector<T>::operator=(const T val)
{
  for(uint i=0; i<_n; i++) 
    _vec[i] = val;
}

template<class T>
Vector<T>& Vector<T>::operator=(const Vector& w)
{
  if (this != &w)
  {
    if (_n != w.size()) error(" bad vector sizes 1");
    for(uint i=0; i<w.size(); i++) _vec[i] = w[i];
  }  
}


template<class T>
Vector<T>& Vector<T>::operator=(const Saxpy<T>& sp)
{
  for(uint i=0; i<_n; i++)
    _vec[i] = sp.a*sp.x[i] + sp.y[i];
  return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector& v)
{
  if (size() != v.size()) error(" bad vector sizes 2");
  for(uint i=0; i<size(); i++) _vec[i] += v[i];
  return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Sax<T>& s)
{
  for(uint i=0; i<size(); i++) 
    _vec[i] += s.a*s.x[i];
  return *this;
}

template<class T>
Vector<T>& Vector<T>::operator-=(const Sax<T>& s)
{
  for(uint i=0; i<size(); i++) 
    _vec[i] -= s.a*s.x[i];
  return *this;
}

template<class T>
Vector<T>& Vector<T>::operator*=(const T val)
{
  for(uint i=0; i<size(); i++) 
    _vec[i] *= val;
  return *this;
}


template<class T>
Vector<T>& Vector<T>::operator-=(const Vector& v)
{
  if (size() != v.size()) error(" bad vector sizes 3");
  for(uint i=0; i<size(); i++) _vec[i] -= v[i];
  return *this;
}

template<class T>
inline T* Vector<T>::ref() const
{
  return _vec;
}

//
// Declaration of the friends
//

template<class T>
inline Vector<T> operator+(const Vector<T> &v)
{
  return v;
}

template<class T>
inline Vector<T> operator-(const Vector<T> &v)
{
  return Vector<T>(v.size()) - v;
}

template<class T>
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2)
{
  if(v1.size() != v2.size()) error(" bad vector sizes 4");
  Vector<T> sum = v1;
  sum += v2;
  return sum;
}

template<class T>
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2)
{
  if(v1.size() != v2.size()) error(" bad vector sizes 5");
  Vector<T> sub = v1;
  sub -= v2;
  return sub;
}

/*
template<class T>
Vector<T> operator*(T scalar, const Vector<T>& v)
{
  Vector<T> tmp(v.size());
  for(uint i=0; i<v.size(); i++)
    tmp(i) = scalar*v[i];
  return tmp;
}
*/

template<class T>
Vector<T> operator*(const Vector<T>& v1, const Vector<T>& v2)
{
  if(v1.size() != v2.size()) error(" bad vector sizes 6");
  uint n = v1.size();
  Vector<T> tmp(n);
  for(uint i=0; i<n; i++) 
    tmp[i]= v1[i]*v2[i];
  return tmp;
}

template<class T>
Vector<T> operator/(const Vector<T>& v, T scalar)
{
  if(!scalar)
    error(" division by zero in vector-scalar division");
  return(1.0/scalar)*v;
}

template<class T>
T Vector<T>::twonorm() const
{
  T norm = _vec[0]*_vec[0];
  for(uint i=1; i<_n; i++)
    norm += _vec[i]*_vec[i];
  return sqrt(norm);
}

template<class T>
T Vector<T>::maxnorm() const
{
  T norm = fabs(_vec[0]);
  for(uint i=0; i<_n; i++)
    norm = max(norm, fabs(_vec[i]));
  return norm;
}

template<class T>
T dot(const Vector<T>& v1, const Vector<T>& v2)
{
  if(v1.size() != v2.size()) error(" bad vector sizes 7");
  T tmp = v1[0] * v2[0];
  for(uint i=1; i<v1.size(); i++)
    tmp += v1[i]*v2[i];
  return sqrt(tmp);
}

template<class T>
ostream& operator<<(ostream& s, const Vector<T>& v)
{
  for(uint i=0; i<15; i++)
    //  for(uint i=0; i<v.size(); i++)
  {
    s << v[i] << "  ";
    if(i%10==9) s << "\n"; // print ten entrie on a line
  }
  return s;    
}

#endif
