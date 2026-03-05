#ifndef TOTAL_LAGRANGIAN_HPP
#define TOTAL_LAGRANGIAN_HPP

#include "nonlinear_elasticity.hpp"

/** This class implements the Total Lagrangian formulation
    for the solution of nonlinear elasticity problems
 */
class TotalLagrangian : public NonlinearElasticity
{
public:

  //! Default constructor
  TotalLagrangian();

  //! Default destructor
  ~TotalLagrangian();

  //! Assemble stiffness matrix
  void assemble_stiff();

  //! Assemble stiffness matrix for the first step (constitutive part only)
  void assemble_const();

  //! Assemble load vector associated with an active (initial) stress state
  void assemble_active(const arma::vec & is,
                       std::vector<arma::mat33*> & vstrs,
                       std::vector<arma::mat33*> & vfib);

  //! Solve the nonlinear problem using Newton's method
  void solve();

  void pre_solve();

  
protected:

  //! Vector of Mixed FE Data handlers
  std::vector<MxFEData*> fedata;  

  //! Compute FE nonlinear matrix
  /*! This is a block matrix with dimension (ndim^2, nnode)
      and it has the following form
             [dN 0   0]
       Bnl = [0  dN  0]
             [0  0  dN]
      where dN is the gradient (gradn) of the shape functions.
   */
  void Bnl_matrix (const arma::mat & gradn, arma::mat & H) const;

  //! Compute FE matrix
  void Bl_matrix (const arma::mat & gradn, const arma::mat * F,
                  arma::mat & Bl) const;

  //! Compute element (initial) stiffness matrix
  void calc_elmat_const (const int e, const MxFE * fe, const Quadrature * qd,
                         arma::mat & elmat);

  //! Compute element stiffness matrix and load vector
  void calc_elmatvec (const int e, const MxFE * fe, const Quadrature * qd,
                      arma::mat & Ke, arma::vec & Re);
  
  //! Compute element residual/reaction vectors
  void elem_resid (const int e, const MxFE * fe, const Quadrature * qd, 
                   arma::vec & Re);
  
  //! Compute element stiffness matrix
  void elem_stiff (const int e, const MxFE * fe, const Quadrature * qd, 
                   arma::mat & Ke);   
  
  //! Pre-compute all the mapping to reference element (one time only)
  void pre_computation();  

  //! Compute the right Cauchy-Green deformation tensor C = F^T F
  void rcg_tensor(const arma::mat & gradn, const vector<arma::vec3> & el_x,
                  double & detf, arma::mat33 * F, arma::mat & C);

};

#endif
