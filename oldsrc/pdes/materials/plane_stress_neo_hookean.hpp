#ifndef PLANE_STRESS_NEO_HOOKEAN_HPP
#define PLANE_STRESS_NEO_HOOKEAN_HPP

#include "hyperelastic_material.hpp"

/** 2D plane stress Neo-Hookean material.
    Reference: J. Bonet, R. Wood, "Nonlinear continuum mechanics for
    finite element analysis", Chapter 8, Page 203
 */
class PlaneStressNeoHookean : public HyperelasticMaterial
{
protected:

  double rho;
  double mu;
  double mu2;
  double thick;

public:

  PlaneStressNeoHookean(const std::vector<double> & prm)
    : HyperelasticMaterial(prm), rho(prm[0]), mu(prm[1]), mu2(0), thick(prm[2])
  {
    name = "PlaneStressNeoHookean";
  }

  double strain_energy(MaterialData * md, const arma::mat &) const { return -1; };

  double get_thick() const { return thick; }

  void cauchy_stress (int ndim, double detf, const arma::mat & B,
		      arma::mat & sigma) const;
  
  void pk2_stress(double detf, const arma::mat & C, arma::mat & S) const;
  
  void elas_tensor (ElasticityType elastype, double detf,
                    arma::mat & D) const;
  
  void mat_elas_tensor (ElasticityType elastype, double detf,
                        const arma::mat & C, arma::mat & D) const;
  
};

#endif
