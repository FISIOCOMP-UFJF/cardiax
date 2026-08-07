#ifndef TENSOR4_H
#define TENSOR4_H

#include <cassert>
#include <iostream>
#include <armadillo>

/** A class for fourth order tensors of arbitrary dimensions.
    Important for Elasticity and NonlinearElasticity problems
    where we need to compute the C_ijkl tensor.
*/
class Tensor4
{
public:

  /// Constructor.
  Tensor4(){}

  /// Destructor
  ~Tensor4(){}

    // : data(new double[81])
    //delete [] data;

  /// Get the corresponding component of the tensor.
  double& operator() (int i, int j, int k, int l)
  {
    int index = get_array_index(i,j,k,l);
    return data[index];
  }

  const double& operator() (int i, int j, int k, int l) const
  {
    int index = get_array_index(i,j,k,l);
    return data[index];
  }

  /// Return the matrix (Voigt) representation of the tensor
  void get_matrix(arma::mat & D);

  /// Print all tensor information
  void view();

  /// Set all components of the tensor to zero.
  void zero();

private:

  // Array to hold the entries of the tensor.
  double data[81];

  // Compute the index of the (i,j,k,l) entry in the array.
  inline int get_array_index(int i, int j, int k, int l) const
  {
    assert(i<3 && j<3 && k<3 && l<3);
    // Row-major order
    return i*27 + j*9 + 3*k + l;
    // Column-major order
    //return i + 3*j + 9*k + 27*l;
  }

};

// Some operations
Tensor4 tensor_product(const arma::mat33 & a, const arma::mat33 & b);
Tensor4 unit_tensor();

/** A class for Symmetric 4th order tensors.

    Due to the major symmetry we can store this tensor as a 6x6 matrix.
    The tensor is stored in column major order:
    
        / 0   1   3   6   10   15  \
        |     2   4   7   11   16  |
        |         5   8   12   17  |
    A = |             9   13   18  |
        |                 14   19  |
        \                      20  /

    Due to the minor symmetry we only store the upper matrix of this tensor.
    Hand notes about this in my notepad.
*/

class Tensor4s 
{
public:

  /// Constructor.
  Tensor4s();

  /// Destructor.
  ~Tensor4s();

  /// Acessor
  double& operator () (int i, int j, int k, int l)
  {    
    const int mm[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
    Tensor4s& T = (*this);
    return T(mm[i][j], mm[k][l]);
  }

  /// Print all tensor information
  void view();
 
private:

  /// Array to hold the entries of the tensor.
  double data[21];

  /// Access to internal data (i,j,k,l) -> (I,J)
  double& operator () (int i, int j)
  {
    const int m[6] = {0, 1, 3, 6, 10, 15};
    return (i<=j) ? data[m[j]+i] : data[m[i]+j];
  }
  
  /// Access to internal data (I,J) -> (k)
  double operator () (int i, int j) const
  {
    const int m[6] = {0, 1, 3, 6, 10, 15};
    return (i<=j) ? data[m[j]+i] : data[m[i]+j];
  }

};

#endif
