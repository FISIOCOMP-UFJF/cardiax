#ifndef MATERIAL_DATA_HPP
#define MATERIAL_DATA_HPP

#include <armadillo>
#include "mesh/mesh.hpp"

/** This class holds information for stress and tangent stiffness 
    computations, such as local deformation, strain, fiber directions
    etc.
*/

class MaterialData
{
public: 

  MaterialData(const Element & el, const arma::mat33 & ftens);

  arma::mat33 get_F() const { return F; }  
  const double get_J()  const { return J; }
  const int get_marker()  const { return marker; }
  const arma::vec3 & fiber()  const { return f; }
  const arma::vec3 & sheet()  const { return s; }
  const arma::vec3 & normal() const { return n; }
  void set_active_stress(double act_stress) {active_stress = act_stress; };
  double get_active_stress() { return active_stress; };

  //! Stabilisation of the velocity-dependent active tension, at this
  //! integration point: the active stiffness Ka and the fibre stretch at the
  //! previous converged mechanics solve. Ka = 0 disables it.
  void set_active_stabilization(double Ka, double lam_prev)
  { active_stiffness = Ka; lambda_prev = lam_prev; };
  double get_active_stiffness() const { return active_stiffness; };
  double get_lambda_prev() const { return lambda_prev; };

  // Kinematics ----------------------------------------------------------------

  //! Calculates Left Cauchy-Green tensor b = F F^T
  arma::mat33 left_cauchy_green(); 

  //! Calculates Right Cauchy-Green tensor C = F^T F
  arma::mat33 right_cauchy_green();

  //! Calculates the Isochoric deformation gradient
  arma::mat33 isochoric_def_grad();

  //! Calculates the Isochoric Left Cauchy-Green tensor
  arma::mat33 isochoric_lcg();

  //! Calculates the Isochoric Right Cauchy-Green tensor
  arma::mat33 isochoric_rcg();
  
  //! Calculates the Lagrangian strain tensor E
  arma::mat33 lagrangian_strain(); 

  
private:
  
  //! Number of dimensions
  int ndim, marker;

  //! Jacobian (det F)
  double J;

  //! Deformation gradient
  arma::mat33 F;

  //! Local material orientation (orthotropic materials)
  arma::vec3 f, s, n;
  double active_stress;

  //! Stabilisation fields; zero by default so that every material that never
  //! calls set_active_stabilization behaves exactly as before.
  double active_stiffness = 0.0;
  double lambda_prev = 0.0;

};

#endif
