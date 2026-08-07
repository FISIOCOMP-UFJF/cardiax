#ifndef VEC_HPP
#define VEC_HPP

#include "mat.hpp"

template<class T, int N>
class TinyVec
{
private:
  
  T data[N];

public:
  
  TinyVec() { for(int i=0; i<N; i++) data[i] = 0.0;}
  TinyVec(T v) { for(int i=0; i<N; i++) data[i] = v; }

  // For Point3d
  TinyVec(T a, T b, T c) { data[0] = a; data[1] = b; data[2] = c; }

  TinyVec(const TinyVec<T,N> & vec) { *this = vec; }

  T & operator[](int i) { return data[i]; }
  const T operator[](int i) const { return data[i]; }

  T & operator()(int i) { return data[i]; }
  const T & operator()(int i) const { return data[i]; }

  //  TinyVec<T,N> & operator=(const T & val)
  //{ for (int i=0; i<N; i++) data[i] = val; return *this; }

  TinyVec<T,N> & operator=(const TinyVec<T,N> & vec)
  { for (int i=0; i<N; i++) data[i] = vec(i); return *this; }
    
  TinyVec<T,N> & operator=(const TinyMat<T,N,1> & mat) 
  { for (int i=0; i<N; i++) data[i] = mat(i,0); return *this; }
  
  TinyVec<T,N> & operator=(const T & d)
  { for (int i=0; i<N; i++) data[i] = d; return *this; }

};

// Useful type definition
typedef TinyVec<double,3> Point3d;
typedef TinyVec<double,2> Point2d;
typedef double Point1d;

// Overload << operator for output
template<class T, int N>
std::ostream & operator<<(std::ostream & out, const TinyVec<T,N> & vec)
{
  for (int i=0; i<N; i++)
    out << vec(i) << " ";
  out << std::endl;
  return out;
}

#endif
