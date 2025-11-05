#ifndef FETOOLS_HPP
#define FETOOLS_HPP

#include <iostream>
#include <map>
#include <vector>
#include <armadillo>
#include "linalg/linalg.h"
#include "fespace.hpp"
#include "mesh/mesh.hpp"

using namespace std;

/*!
 *  Class to hold utilitary static functions for FEM computations
 *  like: computing L2 norm, apply boundary conditions and so on..
*/

class FETools
{
public:

  //! Compute L2 norm of u using the mass matrix M
  static double calc_L2_norm(const arma::sp_mat & M, const arma::vec & u);

  //! Apply boundary conditions
  static void apply_boundary_values(std::map<int,double> nodes,
																		arma::sp_mat & A, arma::vec & b);

  //! Apply boundary conditions
  static void apply_boundary_values(std::map<int,double> nodes,
																		petsc::Matrix & A, petsc::Vector & b);

  //! Apply boundary conditions with given x and update RHS
  static void apply_boundary_values(std::map<int,double> nodes,
																		petsc::Matrix & A, petsc::Vector & b,
																		petsc::Vector & xbar, bool update_rhs);

  //! Apply boundary conditions to matrix A
  static void apply_boundary_values(std::map<int,double> nodes,
																		petsc::Matrix & A);

  //! Post-processing step in order to compute physicald derivative
  //! For more details check the reference:
  //!   "Vigmond and Clements-Construction of a computer model to investigate
  //!    the sawtooth effects in the Purkinje system"
  //!   IEEE Transaction on Biomedical Engineering, Vol.54, No.3, 2007.
	static void cubic_hermite_deriv(const FESpace & fespace, const Mesh & msh,
																	arma::vec & deriv);
};

#endif
