#ifndef EXPTEMPL_H
#define EXPTEMPL_H

#include "vector.h"

// Forward declarations
template<class T> class Vector;

// Expression templates for effiecient vector operations
template<class LeftOpd, class Op, class RightOpd>
struct LOR
{
  LeftOpd lod;
  RightOpd rod;

  LOR(LeftOpd l, RightOpd r) : lod(l), rod(r) {}

  double operator[](int i){
    return Op::apply(lod[i],rod[i]);
  }
};

// Multiplication
struct Multiply
{
  static double apply(double a, double b)
  {
    return a*b;
  }
};

#endif
