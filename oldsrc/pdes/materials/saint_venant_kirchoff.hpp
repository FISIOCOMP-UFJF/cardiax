#ifndef SAINT_VENANT_KIRCHOFF_HPP
#define SAINT_VENANT_KIRCHOFF_HPP

#include "hyperelastic_material.hpp"

/** Saint Venant Kirchoff material */

class SaintVenantKirchoff : public HyperelasticMaterial
{
protected:

  // Shear modulus
  double mu;

  // Bulk modulus
  double lambda;

public:

  SaintVenantKirchoff(const std::vector<double> & prm) :
    HyperelasticMaterial(prm), mu(prm[0]), lambda(prm[1])
  {
    name = "SaintVenantKirchhoff";
  }

  double strain_energy(MaterialData * md, const arma::mat &) const;
  double active_strain_energy(MaterialData * md, const arma::mat & E) const;
  void cauchy_stress(MaterialData * md, arma::mat & sigma) const;
  void piola2_stress(MaterialData * md, arma::mat & S) const;
  void sp_elastensor(MaterialData * md, arma::mat & D) const;
  void mt_elastensor(MaterialData * md, arma::mat & D) const;

};


#endif /* SAINT_VENANT_KIRCHOFF_HPP_ */
