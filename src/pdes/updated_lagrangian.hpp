#ifndef UPDATED_LAGRANGIAN_HPP
#define UPDATED_LAGRANGIAN_HPP

#include "nonlinear_elasticity.hpp"
// #include "materials/neo_hookean.hpp"
// #include "materials/holzapfel_ogden.hpp"

//#define USE_AL

/** 
 * This class implements the Updated Lagrangian formulation
 * for the solution of nonlinear elasticity problems.
 */
class UpdatedLagrangian : public NonlinearElasticity
{
public:

  //! Constructor
  UpdatedLagrangian() {};

  //! Destructor
  ~UpdatedLagrangian() {};

  //! Assemble stiffness matrix
  void assemble_stiff();

  //! Assemble stiffness matrix for the first step (constitutive part only)
  void assemble_const();

  //! Assemble load vector associated with an active (initial) stress state
  //! The assembly of the active stress is performed on the reference
  //! configuration since it is easier....therefore we should use the
  //! reference fiber,sheet and normal directions (f0,s0,n0) and
  //! use the second Piola Kirchhoff stress tensor
  void assemble_active(const arma::vec & is,
                       std::vector<arma::mat33*> & vstrs,
                       std::vector<arma::mat33*> & vfib);

  //! Solve the nonlinear problem using Newton's method
  void solve();

  void pre_solve();

  void solve_old();

protected:


  //! Define if we use Augmented Lagrangian
  bool use_alg;

  //! Store Lambda and (optionally update penalty)
  void al_update(int increment);

  //! Do the augmentation step (update Lagrange multipliers Lambda)
  double al_augment(double tol);

  //! Compute element (initial) stiffness matrix
  void calc_elmat_const (const int e, const MxFE * fe, const Quadrature * qd,
                         arma::mat & elmat);

  //! Calculates the internal nodal forces, geometrical and constitutive
  //! parts of the stiffness matrix (internal, kconst, ksigma)
  void calc_elmatvec (const int e, const MxFE * fe, const Quadrature * qd,
                      arma::mat & Ke, arma::vec & Re);

  //! Computes the B operator
  void calc_B_matrix (const arma::mat & gradn, arma::mat & B) const;

  //! Computes H operator
  void calc_H_matrix (const arma::mat & gradn, arma::mat & H) const;  
  
  //! Compute residual at the element
  void elem_resid (const int iel, const MxFE * fe, const Quadrature * qd,
                   arma::vec & Re);

  //! Compute element stiffness matrix
  void elem_stiff (const int iel, const MxFE * fe, const Quadrature * qd,
                   arma::mat & Ke);

  //! Implementation of the mean dilatation method
  void mean_dilatation(const int iel, const MxFE * fe, const Quadrature * qd,
		       double & p, arma::mat & Ke);

  //! Computes the left Cauchy-Green deformation tensor b = F F^T
  void lcg_tensor (const arma::mat & gradn, const std::vector<arma::vec3> & x0,
                   double * J, arma::mat33 * F, arma::mat & btens);

};

#endif
