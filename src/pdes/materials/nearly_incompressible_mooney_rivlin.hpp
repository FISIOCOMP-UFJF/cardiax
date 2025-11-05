#ifndef NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_HPP
#define NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_HPP

#include "hyperelastic_material.hpp"

/** Nearly incompressible Mooney-Rivlin hyperelastic material.  */

class NearlyIncompressibleMooneyRivlin : public HyperelasticMaterial
{
protected:

  Real c1;
  Real c2;
  Real K;

public:

  NearlyIncompressibleMooneyRivlin(const std::vector<Real> & prm) :
    HyperelasticMaterial(prm), c1(prm[0]), c2(prm[1]), K(prm[2])
  {
    assert(c1 >= 0);
    assert(c2 >= 0);
    assert(c1 > 5*c2);
    name = "NearlyIncompressibleMooneyRivlin";
  }

  Real strain_energy(MaterialData * md, const arma::mat &) const;

  void cauchy_stress (int ndim, Real detf, const arma::mat & B,
		      arma::mat & sigma) const;

  void pk2_stress(Real detf, const arma::mat & C,  arma::mat & S) const;

  void elas_tensor (ElasticityType elastype, Real detf, arma::mat & D) const;
  
  void mat_elas_tensor (ElasticityType elastype, Real detf,
			const arma::mat & C, arma::mat & D) const;
  
};

#endif
