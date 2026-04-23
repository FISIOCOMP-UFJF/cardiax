#ifndef VECOPS_HPP
#define VECOPS_HPP

// -------------------------------------------------------------------
// Vec operations
// -------------------------------------------------------------------

// assignment operations  
// *=
template<class T, int N>
TinyVec<T,N> & operator*=(TinyVec<T,N>& A, T d)
{
  for (int i=0; i<N; i++) A(i) *= d;
  return A;
}    

// /=
template<class T, int N>
TinyVec<T,N> & operator/=(TinyVec<T,N>& A, T d)
{
  for (int i=0; i<N; i++) A(i) /= d;
  return A;
}    

// +=
template<class T, int N>
TinyVec<T,N> & operator+=(TinyVec<T,N>& A, T d)
{
  for (int i=0; i<N; i++) A(i) += d;
  return A;
}    

// +=
template<class T, int N>
TinyVec<T,N> & operator+=(TinyVec<T,N>& A, const TinyVec<T,N>& B)
{
  for (int i=0; i<N; i++) A(i) += B(i);
  return A;
}    

// -=
template<class T, int N>
TinyVec<T,N> & operator-=(TinyVec<T,N>& A, const TinyVec<T,N>& B)
{
  for (int i=0; i<N; i++) A(i) -= B(i);
  return A;
}    


// minus
template<class T, int N>
TinyVec<T,N> operator-(const TinyVec<T,N>& A)
{
  TinyVec<T,N> vec(A);
  for (int i=0; i<N; i++) vec(i)*=-1;
  return vec;
}    

/*
// trans
template<class T, int N>
HTinyVec<T,N> trans(const TinyVec<T,N>& A)
{
  HTinyVec<T,N> hvec;
  for (int i=0; i<N; i++)
      hvec(i)=A(i);
  return hvec;
} 
*/

// matrix T operations
// +
template<class T, int N>
TinyVec<T,N> operator+(const TinyVec<T,N>& A, T d)
{
  TinyVec<T,N> vec(A);
  return vec+= d;
}    

// +
template<class T, int N>
TinyVec<T,N> operator+(T d, const TinyVec<T,N>& A)
{
  TinyVec<T,N> vec(A);
  return vec+= d;
}    

// -
template<class T, int N>
TinyVec<T,N> operator-(const TinyVec<T,N>& A, T d)
{
  TinyVec<T,N> vec(A);
  return vec-= d;
}    

// -
template<class T, int N>
TinyVec<T,N> operator-(T d, const TinyVec<T,N>& A)
{
  TinyVec<T,N> vec;
  for (int i=0; i<N; i++) vec(i) = d-A(i);
  return vec;
}    

// *
template<class T, int N>
TinyVec<T,N> operator*(const TinyVec<T,N>& A, T d)
{
  TinyVec<T,N> vec(A);
  vec *=d;
  return vec;
}    

// *
template<class T, int N>
TinyVec<T,N> operator*(T d, const TinyVec<T,N>& A)
{
  TinyVec<T,N> vec(A);
  vec *=d; 
  return vec;
}    

// /
template<class T, int N>
TinyVec<T,N> operator/(const TinyVec<T,N>& A, T d)
{
  TinyVec<T,N> vec(A);
  return A/= d;
}    

// Matrix - Matrix operations
// +
template<class T, int N>
TinyVec<T,N> operator+(const TinyVec<T,N>& A, const TinyVec<T,N>& B)
{
  TinyVec<T,N> vec(A);
  vec += B;
  return vec;
}    

// -
template<class T, int N>
TinyVec<T,N> operator-(const TinyVec<T,N>& A, const TinyVec<T,N>& B)
{
  TinyVec<T,N> vec(A);
  return vec -= B;
}    

/*
// A*B
template<int M, int N> 
TinyVec<T,N> operator*(const TinyMat<M,N>& A, const TinyVec<N>& v)
{
  TinyVec<T,N> r;
  for (int i=0; i<M; i++)
  {
    T sum=0.0; 
    for (int j=0; j<N; j++)
        sum += A(i,j)*v(j);
      r(i)=sum;
  }
  return r;
}

// A*B
template<int N> 
T operator*(const HTinyVec<N>& v, const TinyVec<N>& w)
{
  T sum=0.0; 
  for (int j=0; j<N; j++)
    sum += v(j)*w(j);
  return sum; 
}
*/

// primary template
template <int I>
struct meta_dot {
  template<class T, int N>
  static T result (TinyVec<T,N> & a, TinyVec<T,N> & b) {
    return a[I]*b[I]  +  meta_dot<I-1>::result(a,b);
  }
};

// partial specialization as end criteria
template <>
struct meta_dot<0> {
  template<class T, int N>
  static T result (TinyVec<T,N> & a, TinyVec<T,N> & b) {
    return a[0]*b[0];
  }
};

// convenience function
template <int N, class T>
inline T dot (TinyVec<T,N> & a, TinyVec<T,N> & b)
{
  return meta_dot<N-1>::result(a,b);
}

// A*v - TinyMat times TinyVec
template<class T, int M, int N> 
TinyVec<T,M> operator*(const TinyMat<T,M,N>& A, const TinyVec<T,N>& v)
{
  TinyVec<T,M> r;
  for (int i=0; i<M; i++)
  {
    T sum=0.0; 
    for (int j=0; j<N; j++)
        sum += A(i,j)*v(j);
      r(i)=sum;
  }
  return r;
}


#endif
