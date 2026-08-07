#ifndef MAT_HPP
#define MAT_HPP

template<class T, int NumRows, int NumCols>
class TinyMat
{
private:
  T data[NumRows*NumCols];
    
public:
  TinyMat() {};
  TinyMat(const TinyMat<T,NumRows,NumCols>& mat) { *this=mat; }
    
  T * operator[](int i) { return data+i*NumCols; }
  const T * operator[](int i) const { return data+i*NumCols; }
    
  T & operator()(int i, int j) { return data[i*NumCols+j]; }
  const T & operator()(int i, int j) const { return data[i*NumCols+j];}
    
  T & operator()(int i) { return data[i]; }
  const T & operator()(int i) const { return data[i]; }
    
  TinyMat<T,NumRows,NumCols> & operator=(const TinyMat<T,NumRows,NumCols>& mat) 
  { for (int i=0; i<NumRows*NumCols; i++) data[i] = mat(i); return *this; }

  TinyMat<T,NumRows,NumCols> & operator=(const T & d)
  { for (int i=0; i<NumRows*NumCols; i++) data[i] = d; return *this; }
};


// Compute determinant

template<class T>
inline T det(const TinyMat<T,1,1>& A)
{ return A(0,0); }

template<class T>
inline T det(const TinyMat<T,2,2>& A)
{ return A(0,0)*A(1,1)-A(1,0)*A(0,1); }

template<class T>
inline T det(const TinyMat<T,3,3>& A)
{
  return A(0,0)*A(1,1)*A(2,2) + A(0,1)*A(1,2)*A(2,0) 
       + A(0,2)*A(1,0)*A(2,1) - A(0,2)*A(1,1)*A(2,0)
       - A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,2)*A(2,1); 
}

// Compute inverse

template<class T>
inline T inv(const TinyMat<T,1,1>& A)
{ return 1.0/A(0,0); }

template<class T>
inline TinyMat<T,2,2> inv(const TinyMat<T,2,2>& A)
{ 
  TinyMat<T,2,2> mat;
  mat(0,0) =  A(1,1); mat(1,1) =  A(0,0);
  mat(0,1) = -A(0,1); mat(1,0) = -A(1,0);
      
  return mat/=det(A);
}

template<class T>
inline TinyMat<T,3,3> inv(const TinyMat<T,3,3>& A)
{ 
  T detA = det(A);
  TinyMat<T,3,3> mat;
  
  mat[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1])/detA;
  mat[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2])/detA;
  mat[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1])/detA;

  mat[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2])/detA;
  mat[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0])/detA;
  mat[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2])/detA;

  mat[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0])/detA;
  mat[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1])/detA;
  mat[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0])/detA;

  return mat;
}

// *********************************************
//     Matrix A * Mat<2,2> Operation
// *********************************************
/*
// A*B<2,2> 
template<class T>
inline Matrix operator*(const Matrix & A, const TinyMat<T,2,2>& B)
{
  if(A.width() != 2){ 
    cout << "vecmatops.hpp: operator* for Matrix * Mat<2,2> error in dimension\n";
    exit(0);
  } 
  
  Matrix mat(A.height(),A.width()); 

  for (int i=0; i<A.height(); i++){
    for (int j=0; j<2; j++){
      double sum=0.0; 
      for (int k=0; k<2; k++)
        sum += A(i,k)*B(k,j);
      mat(i,j)=sum;
    }
  }
  return mat;
}

// A*B<3,3>
inline Matrix operator*(const Matrix & A, const TinyMat<double,3,3>& B)
{
  if(A.width() != 3){ 
    cout << "vecmatops.hpp: operator* for Matrix * Mat<3,3> error in dimension\n";
    exit(0);
  } 
  
  Matrix mat(A.height(),A.width()); 

  for (int i=0; i<A.height(); i++){
    for (int j=0; j<3; j++){
      double sum=0.0; 
      for (int k=0; k<3; k++)
        sum += A(i,k)*B(k,j);
      mat(i,j)=sum;
    }
  }
  return mat;
}

// output: cout << M << endl;
template<class T, int NRow, int NCol>
inline std::ostream& operator<<(std::ostream& o, const TinyMat<T,NRow,NCol> & M)
{
  for (int i=0; i<NRow; i++)
  {
    for (int j=0; j<NCol; j++)
    {
      o << setprecision(8) << (M(i,j)) << "   ";
    }
    o << "\n";
  }
  o << std::flush;
  return o;
}
*/

#endif


