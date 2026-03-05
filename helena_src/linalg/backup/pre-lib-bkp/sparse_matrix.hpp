// A simple sparse matrix class with some matops

#ifndef SPMAT_HPP
#define SPMAT_HPP

#include "array.hpp"
#include "vector.h"
#include "matrix.h"
#include "datatype.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

class SparseMatrix
{
  protected:
    Array<double> vals;
    Array<int> ptrs;
    Array<int> cols;
    int m_,n_;
    bool compressed_;
  public:
    // Constructors
    SparseMatrix();
    // Constructor (Rows, Columns, maximal elements_per_row)
    SparseMatrix(int m, int n, int elprow=10);
    SparseMatrix(const SparseMatrix& B);
    ~SparseMatrix(){};

    // --------------------------------
    // Assignment,Access-Operatoren
    // --------------------------------
    // attention:
    //   if access (i,j) with no entry returns 0, but entry (i,j) is generated
    const double & operator()(int, int) const;
    double& operator()(int, int);
    // Access-Routine
    // where if access to (i,j) with nor entry returns 0, but no entry (i,j) is generated
    double get(int,int) const;


    // -----------------------
    // Info
    // -----------------------
    // Rows
    int m() const {return m_;};
    // Columns
    int n() const {return n_;};
    // Number of Zero Elements
    int nze() const;

    // Output
    std::ostream& print (std::ostream & o) const;

    // Save arrays in binary format
    void save_bin_arrays();
 
    // Matlab-Print
    std::ostream& matlab_print (std::ostream & o) const;

    // Resize SparseMatrix
    SparseMatrix& resize(int m, int n, int elpr=10);

    // Sparse Matrix assignment
    SparseMatrix& operator=(const SparseMatrix&);

    // Multiplication by a scalar
    SparseMatrix& operator*(const Real val);

    // SparseMatrix Multiplication
    Vector<Real> operator*(const Vector<Real> &x) const;

    // SparseMatrix Transpose Multiplication   
    Vector<Real> trans_mult(const Vector<Real> &x) const;   

    // compress SparseMatrix graph
    void compress();

    // Is Matrix compressed?
    bool is_compressed() const { return (compressed_);};

    // Add Small Matrix with indices dnums
    void add(const Matrix & elmat, const vector<int> & dnums);

    void info() const
    {
      cout << "m_ : " << m_ << endl; 
      cout << "n_ : " << n_ << endl; 
      cout << "ptrs " << endl << ptrs << endl; 
      cout << "cols " << endl << cols << endl;
      cout << "vals " << endl << vals << endl;                           
    }     

    // Set the i position of vals to v
    void set_value(int i, Real v) {vals[i]=v;}

    // Return the value of vals in position i
    Real get_value(int i) const {return vals[i];}

    // Set the i position of vals to v
    void set_col(int i, int c) {cols[i]=c;}

    // Return the value of vals in position i
    int get_col(int i) const {return cols[i];}

    // Set the i position of vals to v
    void set_ptr(int i, int p) {ptrs[i]=p;}

    // Return the value of vals in position i
    int get_ptr(int i) const {return ptrs[i];}


};

// output function
std::ostream& operator<<(std::ostream & o, const SparseMatrix & M);

#endif
