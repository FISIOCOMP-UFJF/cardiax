#ifndef NEO_HOOKEAN_HPP
#define NEO_HOOKEAN_HPP

#include "hyperelastic_material.hpp"
#include "incompressible_material.hpp"

/** 
    Nearly incompressible Neo-Hookean hyperelastic material.
 */

class NeoHookean : public IncompressibleMaterial
{
protected:

  //! Lame parameter mu (shear modulus)
  double mu;

  //! Lambda (bulk modulus)
  double lambda;
  
public:

  NeoHookean(const std::vector<double> & prm) :
    IncompressibleMaterial(prm[2]), // kappa
    mu(prm[0]), lambda(prm[1])
  {
    name = "NeoHookean";
    parameters = prm;
  }
  
  double strain_energy(MaterialData * md, const arma::mat &) const;
  double active_strain_energy(MaterialData * md, const arma::mat & E) const;
  void deviatoric_stress(MaterialData * md, arma::mat & stress) const;
  void deviatoric_elastensor(MaterialData * md, Tensor4 & A) const;

};

#endif
