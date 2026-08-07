#ifndef EXPONENTIAL_HPP
#define EXPONENTIAL_HPP

#include "hyperelastic_material.hpp"

/** 
    Exponential material.

    1. Simplification of Holzapfel & Ogden model for passive myocardium
    Ref: "Constitutive modelling of passive myocardium:
          a structurally based framework for 
	  material characterization"
	 Philosophical Transactions of the Royal Society A, 367, 2009.    

    2. Exponential law, eq (3.5), S. Rossi et al
    Ref: "Orthotropic active strain models for the numerical simulation
          of cardiac biomechanics"
	 Int. J. for Numerical Methods in Biomedical Engineering, 28, 2012.
 */

class Exponential : public HyperelasticMaterial
{
protected:

  // Material constants (dimension of stress)
  double a, af, as, afs;

  // Material constants (dimensionless)
  double b, bf, bs, bfs;     

public:

  Exponential (const std::vector<double> & prm) :
    HyperelasticMaterial(prm), 
    a(prm[0]), af(prm[1]),
    as(prm[2]), afs(prm[3]),
    b(prm[4]), bf(prm[5]), 
    bs(prm[6]), bfs(prm[7])
  {
    name = "Exponential";
  }

  double strain_energy(MaterialData * md, const arma::mat &) const;
  double active_strain_energy(MaterialData * md, const arma::mat & E) const;
  void cauchy_stress(MaterialData * md, arma::mat & sigma) const;
  void piola2_stress(MaterialData * md, arma::mat & S) const;
  void sp_elastensor(MaterialData * md, arma::mat & D) const;
  void mt_elastensor(MaterialData * md, arma::mat & D) const;

};


#endif /* EXPONENTIAL_HPP */
