#ifndef MECH_UTILS_HPP
#define MECH_UTILS_HPP

#include <iostream>
#include <armadillo>
#include "linalg/tensor.hpp"

//! Computes the Kronecker delta
inline int delta_kronecker(int i, int j)
{
	return (i==j) ? 1 : 0;
}

//! Levi-Cevita third order alternating (permutation) tensor
//!   1 if (0,1,2), (1,2,0) or (2,1,0)
//!  -1 if (2,1,0), (1,0,2) or (0,2,1)
//!   0 otherwise
inline int levi(int i, int j, int k)
{
	return (j-i)*(k-i)*(k-j)/2;
}

//! Symmetric second order tensor pertubation
void delta_tensor(int k, int l, double eps, arma::mat33 & d);

//! Contraction of two tensors
arma::mat33 contract(const arma::mat33 & A, const arma::mat33 & B);

//! Calculates the deviatoric operator of a second order tensor
arma::mat33 mt_dev_operator(const arma::mat33 & T, const arma::mat33 & C);
arma::mat33 sp_dev_operator(const arma::mat33 & T);

//! Calculates the polar decomposition of the deformation gradient F
void polar_decomposition(const arma::mat33 & F,
						 arma::mat33 & R,
						 arma::mat33 & U);

#endif
