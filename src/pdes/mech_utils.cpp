#include "mech_utils.hpp"

// Auxiliary functions for computing finite difference approximations of 
// second order and fourth order tensors
void delta_tensor(int k, int l, double eps, arma::mat33 & d)
{
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      d(i,j) = 0.5 * (delta_kronecker(i,k) * delta_kronecker(j,l) + 
		                  delta_kronecker(i,l) * delta_kronecker(j,k) ) * eps;
}

arma::mat33 contract(const arma::mat33 & A, const arma::mat33 & B)
{
  arma::mat33 R;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      R(i,j) = A(i,j) * B(i,j);
    }
  return R;
}

arma::mat33 mt_dev_operator(const arma::mat33 & T, const arma::mat33 & C)
{
  return T - (1./3.) * (contract(T,C) * C.i());
}

arma::mat33 sp_dev_operator(const arma::mat33 & T)
{
  arma::mat33 I = arma::eye(3,3);
  return T - (1./3.) * (contract(T,I) * I);
}

void polar_decomposition(const arma::mat33 & F,
			 arma::mat33 & R, 
			 arma::mat33 & U)
{
  arma::mat33 Ui, U2, eigvec;
  arma::vec3 eigval;

  // Right Cauchy-Green deformation C (or U^2)
  U2 = F.t() * F;
  
  // Spectral decomposition
  eig_sym(eigval, eigvec, U2);
  eigval = arma::sqrt(eigval);

  // Compute U = sum_i lambda_i Ni Nj
  U.zeros();
  for(int i=0; i<3; i++)
    U += eigval(i) * ( eigvec.col(i) * eigvec.col(i).t() );

  // F=RU ==> R = F U^{-1}
  Ui = arma::inv(U);
  R  = F * Ui;
}
