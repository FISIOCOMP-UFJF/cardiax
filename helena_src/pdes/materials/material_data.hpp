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

};

#endif
