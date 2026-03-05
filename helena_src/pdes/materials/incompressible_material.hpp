#ifndef INCOMPRESSIBLE_MATERIAL_HPP
#define INCOMPRESSIBLE_MATERIAL_HPP

#include <string>
#include <stdexcept>
#include "util/util.hpp"
#include "linalg/tensor.hpp"
#include "pdes/elasticity.hpp"
#include "pdes/mech_utils.hpp"
#include "material_data.hpp"
#include "hyperelastic_material.hpp"

class IncompressibleMaterial : public HyperelasticMaterial
{
public:

  //! Constructor
  IncompressibleMaterial(double K);

  //! Destructor
  virtual ~IncompressibleMaterial(){};

  //! Return the bulk modulus
  double get_kappa() const { return kappa; }

  //! Strain energy function
  double strain_energy(MaterialData * md, const arma::mat &) const;

  double active_strain_energy(MaterialData * md, const arma::mat & E) const;

  //! Add pressure to the stress tensor 
  void add_pressure(double press, arma::mat & sigma);
  
  //! Cauchy stress tensor
  void cauchy_stress(MaterialData * md, arma::mat & sigma) const;

  //! Second Piola-Kirchhoff stress tensor
  void piola2_stress(MaterialData * md, arma::mat & S) const;

  //! Spatial elasticity (4th order) tensor
  void sp_elastensor(MaterialData * md, arma::mat & D) const;

  //! Material elasticity (4th order) tensor
  void mt_elastensor(MaterialData * md, arma::mat & D) const;

  //! Volumetric component of the elasticity tensor matrix
  void sp_volumetric_elastensor(const double pressure, Tensor4 & A) const;

  void copy_elastensor(Tensor4 B, Tensor4 & A) const;

  void copy_stress(arma::mat A, arma::mat & S) const;


  // Interface -----------------------------------------------------------------

  virtual void deviatoric_stress(MaterialData * md, arma::mat & stress) const=0;

  virtual void deviatoric_elastensor(MaterialData * md, Tensor4 & A) const=0;

  // Some helpful functions ----------------------------------------------------

  //double dUdJ   (double J) { return log(J)/J; }
  //double d2UdJJ (double J) { return (1.0-log(J))/(J*J); }
  //double h      (double J) { return log(J); }
  //double dhdp   (double J) { return 1.0 / J; }
  //double d2hdpp (double J) { return -1.0 / (J*J); }

  double dUdJ   (double J) { return (J-1.0); }
  double d2UdJJ (double J) { return 1; }
  double h      (double J) { return (J-1.0); }
  double dhdp   (double J) { return 1.0; }
  double d2hdpp (double J) { return 0.0; }

protected:

  //! Bulk modulus
  double kappa;

};

#endif
