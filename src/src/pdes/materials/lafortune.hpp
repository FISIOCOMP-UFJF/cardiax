#ifndef LAFORTUNE_HPP
#define LAFORTUNE_HPP

#include "hyperelastic_material.hpp"

/** LaFortune et al model
 */

class LaFortune : public HyperelasticMaterial
{
protected:

  /// Material constants (dimension of stress)
  double a, af, as, afs;

  /// Material constants (dimensionless)
  double b, bf, bs, bfs;     

  /// Auxiliary
  const arma::mat33 I;

public:

  LaFortune (const std::vector<double> & prm) :
    HyperelasticMaterial(prm), 
    a(prm[0]), af(prm[1]),
    as(prm[2]), afs(prm[3]),
    b(prm[4]), bf(prm[5]), 
    bs(prm[6]), bfs(prm[7]),
    I(arma::eye<arma::mat>(3,3))
  {
    assert(a >= 0);
    assert(b >= 0);
    assert(af >= 0);
    assert(bf >= 0);
    assert(as >= 0);
    assert(bs >= 0);
    assert(afs >= 0);
    assert(bfs >= 0);
    name = "LaFortune";
  }
    
  double strain_energy(MaterialData * md, const arma::mat &) const;
  void cauchy_stress(MaterialData * md, arma::mat & sigma) const;
  void piola2_stress(MaterialData * md, arma::mat & S) const;
  void sp_elastensor(MaterialData * md, arma::mat & D) const;
  void mt_elastensor(MaterialData * md, arma::mat & D) const;

};


#endif /* LAFORTUNE_HPP */
