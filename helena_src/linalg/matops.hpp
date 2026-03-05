#ifndef MATOPS_HPP
#define MATOPS_HPP

// -------------------------------------------------------------------
// Mat operations
// -------------------------------------------------------------------
  
// assignment operations  
// *=
template<class T, int M, int N>
TinyMat<T,M,N> & operator*=(TinyMat<T,M,N>& A, T d)
{
  for (int i=0; i<M*N; i++) A(i) *= d;
  return A;
}    

// /=
template<class T, int M, int N>
TinyMat<T,M,N> & operator/=(TinyMat<T,M,N>& A, T d)
{
  for (int i=0; i<M*N; i++) A(i) /= d;
  return A;
}    

// +=
template<class T, int M, int N>
TinyMat<T,M,N> & operator+=(TinyMat<T,M,N>& A, T d)
{
  for (int i=0; i<M*N; i++) A(i) += d;
  return A;
}    

// +=
template<class T, int M, int N>
TinyMat<T,M,N> & operator+=(TinyMat<T,M,N>& A, const TinyMat<T,M,N>& B)
{
  for (int i=0; i<M*N; i++) A(i) += B(i);
  return A;
}    

// -=
template<class T, int M, int N>
TinyMat<T,M,N> & operator-=(TinyMat<T,M,N>& A, const TinyMat<T,M,N>& B)
{
  for (int i=0; i<M*N; i++) A(i) -= B(i);
  return A;
}    


// minus
template<class T, int M, int N>
TinyMat<T,M,N> operator-(const TinyMat<T,M,N>& A)
{
  TinyMat<T,M,N> mat(A);
  for (int i=0; i<M*N; i++) mat(i)*=-1;
  return mat;
}    

// trans
template<class T, int M, int N>
TinyMat<T,N,M> trans(const TinyMat<T,M,N>& A)
{
  TinyMat<T,N,M> mat;
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)
      mat(j,i)=A(i,j);
  return mat;
}    

// matrix T operations
// +
template<class T, int M, int N>
TinyMat<T,M,N> operator+(const TinyMat<T,M,N>& A, T d)
{
  TinyMat<T,M,N> mat(A);
  return mat+= d;
}    

// +
template<class T, int M, int N>
TinyMat<T,M,N> operator+(T d, const TinyMat<T,M,N>& A)
{
  TinyMat<T,M,N> mat(A);
  return mat+= d;
}    

// -
template<class T, int M, int N>
TinyMat<T,M,N> operator-(const TinyMat<T,M,N>& A, T d)
{
  TinyMat<T,M,N> mat(A);
  return mat-= d;
}    

// -
template<class T, int M, int N>
TinyMat<T,M,N> operator-(T d, const TinyMat<T,M,N>& A)
{
  TinyMat<T,M,N> mat;
  for (int i=0; i<M*N; i++) mat(i) = d-A(i);
  return mat;
}    

template<class T, int M, int N>
TinyMat<T,M,N> operator*(const TinyMat<T,M,N>& A, T d)
{
  TinyMat<T,M,N> mat(A);
  return A*= d;
}    

// *
template<class T, int M, int N>
TinyMat<T,M,N> operator*(T d, const TinyMat<T,M,N>& A)
{
  TinyMat<T,M,N> mat(A);
  return A*= d;
}    

// /
template<class T, int M, int N>
TinyMat<T,M,N> operator/(const TinyMat<T,M,N>& A, T d)
{
  TinyMat<T,M,N> mat(A);
  return A/= d;
}    

// TinyMatrix - TinyMatrix operations
// +
template<class T, int M, int N>
TinyMat<T,M,N> operator+(const TinyMat<T,M,N>& A, const TinyMat<T,M,N>& B)
{
  TinyMat<T,M,N> mat(A);
  mat += B;
  return mat;
}    

// -
template<class T, int M, int N>
TinyMat<T,M,N> operator-(const TinyMat<T,M,N>& A, const TinyMat<T,M,N>& B)
{
  TinyMat<T,M,N> mat(A);
  return mat -= B;
}    

// A*B
template<class T, int M, int K, int N> 
TinyMat<T,M,N> operator*(const TinyMat<T,M,K>& A, const TinyMat<T,K,N>& B)
{
  TinyMat<T,M,N> mat;
  for (int i=0; i<M; i++)
  {
    for (int j=0; j<N; j++)
    {
      T sum=0.0; 
      for (int k=0; k<K; k++)
        sum += A(i,k)*B(k,j);
      mat(i,j)=sum;
    }
  }
  return mat;
}

#endif
