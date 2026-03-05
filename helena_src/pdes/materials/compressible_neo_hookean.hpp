#ifndef COMPRESSIBLE_NEO_HOOKEAN_HPP
#define COMPRESSIBLE_NEO_HOOKEAN_HPP

#include "hyperelastic_material.hpp"

/** 3D or plane strain compressible Neo-Hookean hyperelastic material.
    Reference: J. Bonet, R. Wood, "Nonlinear continuum mechanics for
    finite element analysis", Chapter 8, Page 201
 */

class CompressibleNeoHookean : public HyperelasticMaterial
{
protected:

  // Shear modulus
  double mu;

  // Bulk modulus
  double lambda;

public:

  CompressibleNeoHookean(const std::vector<double> & prm) :
    HyperelasticMaterial(prm), mu(prm[0]), lambda(prm[1])
  {
    assert(mu >= 0);
    assert(lambda >= 0);
    name = "CompressibleNeoHookean";
  }

  double strain_energy(MaterialData * md, const arma::mat &) const;
  double active_strain_energy(MaterialData * md, const arma::mat & E) const;
  void cauchy_stress(MaterialData * md, arma::mat & sigma) const;
  void piola2_stress(MaterialData * md, arma::mat & S) const;
  void sp_elastensor(MaterialData * md, arma::mat & D) const;
  void mt_elastensor(MaterialData * md, arma::mat & D) const;

};

#endif /* COMPRESSIBLE_NEO_HOOKEAN_HPP */
