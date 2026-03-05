// Auxiliary classes for efficient implementation of the SAXPY operation
// using an expression similar to the mathematical one.
// This technique is called deferred-evaluation.

#ifndef SAXPY_H
#define SAXPY_H

#include "vector.h"
#include "datatype.h"

// Forward declaration
template<class T> class Vector;

template<class T>
struct Sax
{
  const Real& a;
  const Vector<T>& x;
  Sax(const Real& d, const Vector<T>& v) : a(d), x(v) {}
};

template<class T>
inline Sax<T> operator*(const Real& d, const Vector<T>& v)
{
  return Sax<T>(d,v);
}

template<class T>
struct Saxpy
{
  const Real& a;
  const Vector<T>& x;
  const Vector<T>& y;
  Saxpy(const Sax<T>& s, const Vector<T>& u) : a(s.a), x(s.x), y(u) {}
};

template<class T>
inline Saxpy<T> operator+(const Sax<T>& s, const Vector<T>& u)
{
  return Saxpy<T>(s,u);
}

template<class T>
inline Saxpy<T> operator+(const Vector<T>& u, const Sax<T>& s)
{
  return Saxpy<T>(s,u);
}

#endif
