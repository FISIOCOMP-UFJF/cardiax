#ifndef SPMAT_HPP
#define SPMAT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

//#include "array.hpp"
//#include "vector.hpp"

/*
class SparseMatrix
{
public:
  SparseMatrix();
  SparseMatrix(int m, int n, int elprow=10);
  SparseMatrix(const SparseMatrix& B);
  ~SparseMatrix(){};

  // attention:
  //   if access (i,j) with no entry returns 0, but entry (i,j) is generated
  const double & operator()(int, int) const;
  double& operator()(int, int);

  // Access-Routine
  // where if access to (i,j) with nor entry returns 0, but no entry (i,j) is generated
  double get(int,int) const;

  int m() const {return m_;};
  // Columns
  int n() const {return n_;};
  // Number of Zero Elements
  int nze() const;

  // Output
  std::ostream& Print (std::ostream & o) const;

  // Save arrays in binary format
  void saveBinaryArrays();
 
  // Matlab-Print
  std::ostream& MatlabPrint (std::ostream & o) const;   

  // ----------------------
  // Matrix-Operations
  // ----------------------
  // Resize SparseMatrix
  SparseMatrix& resize(int m, int n, int elpr=10);

  //void clear(){resize(0,0,0)};

  // Sparse Matrix assignment
  SparseMatrix& operator=(const SparseMatrix&);

  // SparseMatrix Multiplication
  Vector operator*(const Vector &x) const;

  // SparseMatrix Transpose Multiplication   
  Vector trans_mult(const Vector &x) const;   

  // compress SparseMatrix graph
  void compress();

  // Is Matrix compressed?
  bool isCompressed() const { return (compressed_);};

  // Add Small Matrix with indices dnums
  void Add(const Matrix & elmat, const vector<int> & dnums);

  // Set the i position of vals to v
  void set_value(int i, double v) {vals[i]=v;}

  // Return the value of vals in position i
  double get_value(int i) const {return vals[i];}

  // Set the i position of vals to v
  void set_col(int i, int c) {cols[i]=c;}

  // Return the value of vals in position i
  int get_col(int i) const {return cols[i];}

  // Set the i position of vals to v
  void set_ptr(int i, int p) {ptrs[i]=p;}

  // Return the value of vals in position i
  int get_ptr(int i) const {return ptrs[i];}    

  void Info() const
  {
    cout << "m_ : " << m_ << endl; 
    cout << "n_ : " << n_ << endl; 
    cout << "ptrs " << endl << ptrs << endl; 
    cout << "cols " << endl << cols << endl;
    cout << " vals " << endl << vals << endl;                           
  }     

protected:
  Array<double> vals;
  Array<int> ptrs;
  Array<int> cols;
  int m_,n_;
  bool compressed_;

};

// output function
std::ostream& operator<<(std::ostream & o, const SparseMatrix & M);

// --------------------------------------------------------------
// IMPLEMENTATION  (-> better to *.cpp)
// --------------------------------------------------------------

// Constructors
inline SparseMatrix::SparseMatrix()
  : vals(), ptrs(), cols(), m_(0), n_(0), compressed_(0)
{}

inline SparseMatrix::SparseMatrix(int m, int n, int elprow)
  : vals(m*elprow+1), ptrs(m+1), cols(m*elprow+1), m_(m), n_(n), compressed_(0)
{
  for(int i=0;i<m+1;i++)
    ptrs[i] = i*elprow;

  cols=-1;

  // vals=0.; // this initialization is not necessary    
  vals[ptrs[m_]] = 0; // only this is used as dummy value    
}

// copy-Constructor
inline SparseMatrix::SparseMatrix(const SparseMatrix& M)
  : vals(M.vals),  ptrs(M.ptrs), cols(M.cols), m_(M.m_), n_(M.n_),
    compressed_(M.compressed_)
{}


inline double SparseMatrix::get(int i, int j) const
{
  if( i<0 || j <0 || i >= m_ || j >= n_)
    {
      std::cerr << " Index ( " << i << "," << j << ") in Element access is out of range --> access ignored" << endl;

      return(0);
    }

  for(int k=ptrs[i]; k<ptrs[i+1]; k++)
    { 
      if(cols[k]==j)
        return(vals[k]);
      else if(cols[k]==-1) break;   
    }
    
  return(0);
}

inline const double & SparseMatrix::operator()(int i, int j) const
{
  if( i<0 || j <0 || i >= m_ || j >= n_)
  {
    std::cerr << " Index ( " << i << "," 
	      << j << ") in Element access is out of range"
	      << " --> access ignored -- AQUI" << endl;
    return(vals[ptrs[m_]]);
  }

  for(int k=ptrs[i]; k<ptrs[i+1]; k++)   
  {
    if(cols[k]==j) 
      return(vals[k]);
    else if (cols[k]==-1) break;   
  }
    
  return(vals[ptrs[m_]]);
}


// Element Access
inline double& SparseMatrix::operator()(int i, int j)
{
  if( i<0 || j <0 || i >= m_ || j >= n_)
    {
      cerr << " Index ( " << i << "," << j << ") in Element assignement is out of range --> assignement ignored -- OPERATOR()" << endl;
      cerr << "m_=" << m_ << "\tn_=" << n_ << endl;

      return(vals[ptrs[m_]]); // assigned to last dummy value
    }

  for (int k=ptrs[i]; k<ptrs[i+1]; k++)
    {
      if (cols[k] == j)
	return (vals[k]);
      else
	if (cols[k] == -1)
	  {
	    cols[k] = j;
	    vals[k]=0.0;      
	    return (vals[k]);
	  }
    }
  this->Print(cout);     
  std::cerr << "Array element (" << i << "," << j << ") cannot be directly accessed/assigned to sparse structure -- row full."
	    << endl;
 
  exit(1);
  return vals[ptrs[m_]];   // return to suppress compiler warning message and exit
}


// ReSize
inline SparseMatrix& SparseMatrix::resize(int m, int n, int elpr)
{
  if (m < 0 || n < 0 || elpr <0)
    {
      std::cerr << "error in SparseMatrix::resize(int m, int n, int elpr): size negative!\n";
      exit(1);
    }
  if (vals.size() != m*elpr+1)
    {
      vals.resize(m*elpr+1);
      cols.resize(m*elpr+1);
      ptrs.resize(m+1);
    }

  for(int i=0;i<m+1;i++)
    ptrs[i] = i*elpr;

  cols=-1;

  //vals=0.;         
  vals[ptrs[m]]=0.0; // this value is used if an entry is not yet set      

  m_ = m;
  n_ = n;
  return *this;
}


// ----------------------------------
// in and output functions
// ----------------------------------

inline std::ostream & SparseMatrix::Print (std::ostream & o) const
{
  o << "SparseMatrix: " << endl;
  for(int i = 0; i<m_; i++)
    {
      o << "Row " << i << "::" ;
      for(int k=ptrs[i]; k<ptrs[i+1]; k++)
	{
	  if(cols[k] != -1)
	    o << "(" << i << "," << cols[k] << "):"  << vals[k] << "\t" ;
	}
      o << endl;
    }
  return(o);
}

  
// Printout for matlab (indices 1-based, A  = sparse(i,j,vals) )
inline std::ostream & SparseMatrix::MatlabPrint (std::ostream & o) const
{
  for(int i = 0; i<m_; i++)
    for(int k=ptrs[i]; k<ptrs[i+1]; k++)
      if(cols[k] != -1)
	o  << i+1 << " \t " << cols[k]+1 << " \t "  << vals[k] << endl ;
  return(o);
}

  
inline void SparseMatrix::saveBinaryArrays ()
{
  std::ofstream array_vals;
  std::ofstream array_cols;
  std::ofstream array_ptrs;

  array_vals.open("sparse_vals.bin",ios::out|ios::binary);
  array_cols.open("sparse_cols.bin",ios::out|ios::binary);
  array_ptrs.open("sparse_ptrs.bin",ios::out|ios::binary);

  int size_nzd = sizeof(double) * nze();
  int size_nzi = sizeof(int)    * nze();
  int size_ptr = sizeof(int)    * n()+1;

  int nz = nze();
  int nn = n()+1;

  if (array_vals.is_open()){
    array_vals.write ((char*)&nz, sizeof(int));
    array_vals.write ((char*)vals.ref(), size_nzd);
    array_vals.close ();
  } else {
    cout << "Error saving binary of SparseMatrix." << endl; exit(1);
  }

  if (array_cols.is_open()){
    array_cols.write ((char*)&nz, sizeof(int));
    array_cols.write ((char*)cols.ref(), size_nzi);
    array_cols.close ();
  } else {
    cout << "Error saving binary of SparseMatrix." << endl; exit(1);
  }

  if (array_ptrs.is_open()){
    array_ptrs.write ((char*)&nn, sizeof(int));
    array_ptrs.write ((char*)ptrs.ref(), size_ptr);
    array_ptrs.close ();
  } else {
    cout << "Error saving binary of SparseMatrix." << endl; exit(1);
  }
}


// output: cout << v << endl;
inline std::ostream& operator<<(std::ostream& o, const SparseMatrix& M)
{
  return(M.Print(o));
}

// ----------------------------
// MathOps
// ----------------------------

inline int SparseMatrix::nze() const
{
  if(compressed_) return(vals.size()-1);

  int nz = 0;
  for(int k=0;k<vals.size()-1;k++)
    if(cols[k] != -1) nz++;

  return(nz);
}


inline void SparseMatrix::compress()
{
  if(compressed_) return;

  int nz = nze();

  Array<double> valsc(nz+1);
  Array<int> colsc(nz+1);
  int knew = 0;
  int ptr = ptrs[0];
  for(int i=0;i<m_;i++)
    {
      for(int k=ptr;k<ptrs[i+1];k++)
	{
	  if(cols[k] == -1) break;
	  colsc[knew] = cols[k];
	  valsc[knew] = vals[k];
	  knew++;
	}
      ptr = ptrs[i+1];
      ptrs[i+1] = knew;
    }
  // automatically resized by copy of array-class
  vals = valsc;
  cols = colsc;

  compressed_ = 1;
  return;
}


// y=M*x
inline Vector SparseMatrix::operator*(const Vector &x) const
{
  if(x.size() != n_)
    {
      std::cerr << "ERROR: SparseMatrix::operator* dimensions not matching : x.size() = " << x.size() << " while A.n() = " << n_ << "\n ";
      exit(1);
    }

  Vector y(m_);
  y=0.;

  if(!compressed_)
    {
#ifdef DEBUG
      std::cerr <<"WARNING: It should be cheaper to compress SparseMatrix before Multiplication \n";
#endif
      for(int i=0;i<m_;i++)
	for(int k=ptrs[i]; k<ptrs[i+1]; k++)
	  {
	    if(cols[k]==-1) break;
	    y(i)+= vals[k] * x(cols[k]);
	  }
      return (y);
    }

  for(int i=0;i<m_;i++)
    for(int k=ptrs[i]; k<ptrs[i+1]; k++)
      y(i)+= vals[k] * x(cols[k]);


  return (y);
}

// y=M^T*x
inline Vector SparseMatrix::trans_mult(const Vector &x) const
{
  if(x.size() != m_)
    {
      std::cerr << "ERROR: SparseMatrix::trans_mult dimensions not matching : x.size() = " << x.size() << " while A.m() = " << n_ << "\n ";
      exit(1);
    }

  Vector y(n_);
  y=0.;

  if(!compressed_)
    {
#ifdef DEBUG
      std::cerr <<"WARNING: It should be cheaper to compress SparseMatrix before Multiplication \n";
#endif
      for(int i=0;i<m_;i++)
	for(int k=ptrs[i]; k<ptrs[i+1]; k++)
	  {
	    if(cols[k]==-1) break;
	    y(cols[k])+= vals[k] * x(i);
	  }
      return (y);
    }

  for(int i=0;i<m_;i++)
    for(int k=ptrs[i]; k<ptrs[i+1]; k++)
      y(cols[k])+= vals[k] * x(i);


  return (y);
}   


inline void SparseMatrix::Add(const Matrix & elmat, const vector<int> & dnums)
{
  for(int k=0;k< (int)dnums.size(); k++)
    for(int l=0; l< (int)dnums.size(); l++)
      (*this)(dnums[k],dnums[l])+=elmat(k,l);
}
*/

#endif
