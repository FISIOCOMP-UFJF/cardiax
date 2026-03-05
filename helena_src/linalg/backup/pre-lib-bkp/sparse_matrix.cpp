#include "sparse_matrix.hpp"

// Constructors
SparseMatrix::SparseMatrix() 
  : vals(), ptrs(), cols(), m_(0), n_(0), compressed_(0)
{}

SparseMatrix::SparseMatrix(int m, int n, int elprow) 
  : vals(m*elprow+1), ptrs(m+1), cols(m*elprow+1), 
    m_(m), n_(n), compressed_(0)
{
  for(int i=0;i<m+1;i++)
    ptrs[i] = i*elprow;

  cols=-1;

  // vals=0.; // this initialization is not necessary    
  vals[ptrs[m_]] = 0; // only this is used as dummy value    
}

// copy-Constructor
SparseMatrix::SparseMatrix(const SparseMatrix& M) 
  : vals(M.vals),  ptrs(M.ptrs), cols(M.cols), 
    m_(M.m_), n_(M.n_), compressed_(M.compressed_)
{}


double SparseMatrix::get(int i, int j) const
{
  if( i<0 || j <0 || i >= m_ || j >= n_)
  {
    std::cerr << " Index ( " << i << "," << j << 
      ") in Element access is out of range --> access ignored" << endl;
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

const double & SparseMatrix::operator()(int i, int j) const
{
  if( i<0 || j <0 || i >= m_ || j >= n_)
  {
    std::cerr << " Index ( " << i << "," << j 
	      << ") in Element access is out of range --> access ignored\n";
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
double& SparseMatrix::operator()(int i, int j)
{
  if( i<0 || j <0 || i >= m_ || j >= n_)
  {
    std::cerr << " Index ( " << i << "," << j << 
      ") in Element assignement is out of range --> assignement ignored\n";
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

  this->print(cout);     
  std::cerr << "Array element (" << i << "," << j 
	    << ") cannot be directly accessed/assigned to sparse structure"
	    << "-- row full."
	    << endl;
  exit(1);

  // return to suppress compiler warning message and exit
  return vals[ptrs[m_]]; 
}

SparseMatrix& SparseMatrix::resize(int m, int n, int elpr)
{
  if (m < 0 || n < 0 || elpr <0)
  {
    std::cerr << "error in SparseMatrix::resize(int m, int n, int elpr):"
	      << " size negative!\n";
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

// -------------------------------------------------------------------
// I/O functions
// -------------------------------------------------------------------

std::ostream & SparseMatrix::print (std::ostream & o) const
{
  /*
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
  }*/
  for (int i=0; i<30; i++)
    cout << vals[i] << " ";
  cout << endl;
  return(o);
}

  
// Printout for matlab (indices 1-based, A  = sparse(i,j,vals) )
std::ostream & SparseMatrix::matlab_print (std::ostream & o) const
{
  for(int i = 0; i<m_; i++)
    for(int k=ptrs[i]; k<ptrs[i+1]; k++)
      if(cols[k] != -1)
	o  << i+1 << " \t " << cols[k]+1 << " \t "  << vals[k] << endl ;
  return(o);
}

  
void SparseMatrix::save_bin_arrays ()
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
std::ostream& operator<<(std::ostream& o, const SparseMatrix& M)
{
  return(M.print(o));
}

// -------------------------------------------------------------------
// MathOps
// -------------------------------------------------------------------

int SparseMatrix::nze() const
{
  if(compressed_) return(vals.size()-1);

  int nz = 0;
  for(int k=0;k<vals.size()-1;k++)
    if(cols[k] != -1) nz++;

  return(nz);
}


void SparseMatrix::compress()
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

// A = M
SparseMatrix & SparseMatrix::operator=(const SparseMatrix &M)
{
  m_ = M.m();
  n_ = M.n();

  for(int i=0; i<M.nze(); i++)
  {
    cols[i] = M.get_col(i);
    vals[i] = M.get_value(i);
  }

  for(int i=0; i<M.m(); i++)
    ptrs[i] = M.get_ptr(i);
  
  return *this;
}


// y=scalar*M
/*
SparseMatrix& SparseMatrix::operator*(const real val)
{ 
  for(int i=0; i<nze(); i++)
      vals[i] = val * vals[i];

  return *this;
}*/

// y=M*x
Vector<Real> SparseMatrix::operator*(const Vector<Real> &x) const
{
  if(x.size() != n_)
  {
    std::cerr << "ERROR: SparseMatrix::operator* dimensions not matching :"
	      << "x.size() = " << x.size() << " while A.n() = " << n_ 
	      << "\n ";
    exit(1);
  }

  Vector<Real> y(m_, 0.0);

  if(!compressed_)
  {
#ifdef DEBUG
    std::cerr <<"WARNING: It should be cheaper to compress SparseMatrix"
	      << "before Multiplication \n";
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

  return y;
}

// y=M^T*x
Vector<Real> SparseMatrix::trans_mult(const Vector<Real> &x) const
{
  if(x.size() != m_)
  {
    std::cerr << "ERROR: SparseMatrix::trans_mult dimensions not matching"
	      << " : x.size() = " << x.size() << " while A.m() = " 
	      << n_ << "\n ";
    exit(1);
  }

  Vector<Real> y(n_, 0.0);

  if(!compressed_)
  {
#ifdef DEBUG
    std::cerr <<"WARNING: It should be cheaper to compress SparseMatrix"
	      << "before Multiplication \n";
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

  return y;
}   


void SparseMatrix::add(const Matrix & elmat, const vector<int> & dnums)
{
  /*
    for(int k=0;k< (int)dnums.size(); k++)
    for(int l=0; l< (int)dnums.size(); l++)
    (*this)(dnums[k],dnums[l])+=elmat(k,l);
  */
}
