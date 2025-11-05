#ifndef HOLZAPFEL_OGDEN_HPP
#define HOLZAPFEL_OGDEN_HPP

#include "hyperelastic_material.hpp"
#include "incompressible_material.hpp"

/*!
 * Holzapfel & Ogden model for passive myocardium
 * 
 * Ref: "Constitutive modelling of passive myocardium:
 *       a structurally based framework for material characterization"
 *       Philosophical Transactions of the Royal Society A, 367, 2009.
 *   
*/

class HolzapfelOgden : public IncompressibleMaterial
{
protected:

  //! Material constants (dimension of stress)
  std::vector<double> v_a, v_af, v_as, v_afs;

  //! Material constants (dimensionless)
  std::vector<double> v_b, v_bf, v_bs, v_bfs;     

  //! Auxiliary
  const arma::mat33 I;

  std::vector<int> map_mat;

public:
  HolzapfelOgden (const std::vector<std::vector<double>> & prm, int num_materials, std::vector<int> & map_AHAmat) : IncompressibleMaterial(prm[0][9]), I(arma::eye<arma::mat>(3,3))
  {
//    std::cout << " Passou \n";
    for (int ii=0; ii<num_materials; ii++)
    {
      v_a.push_back(prm[ii][0]);
      v_af.push_back(prm[ii][1]);
      v_as.push_back(prm[ii][2]);
      v_afs.push_back(prm[ii][3]);
      v_b.push_back(prm[ii][4]);
      v_bf.push_back(prm[ii][5]);
      v_bs.push_back(prm[ii][6]);
      v_bfs.push_back(prm[ii][7]);
    }

    name = "HolzapfelOgden";
    active_stress = prm[0][8];

    //parameters = prm[0];
    parameters.push_back(prm[0][0]);
    parameters.push_back(prm[0][1]);
    parameters.push_back(prm[0][2]);
    parameters.push_back(prm[0][3]);
    parameters.push_back(prm[0][4]);
    parameters.push_back(prm[0][5]);
    parameters.push_back(prm[0][6]);
    parameters.push_back(prm[0][7]);
    parameters.push_back(prm[0][8]);
    parameters.push_back(prm[0][9]);

    map_mat = map_AHAmat;

  }

public:

  HolzapfelOgden (const std::vector<double> & prm) :    
    IncompressibleMaterial(prm[9]), // kappa
    I(arma::eye<arma::mat>(3,3))
  {

    v_a.push_back(prm[0]);
    v_af.push_back(prm[1]);
    v_as.push_back(prm[2]);
    v_afs.push_back(prm[3]);
    v_b.push_back(prm[4]);
    v_bf.push_back(prm[5]);
    v_bs.push_back(prm[6]);
    v_bfs.push_back(prm[7]); 
    name = "HolzapfelOgden";
    active_stress = prm[8];
    parameters = prm;
    map_mat.push_back(0);
  }

  double strain_energy(MaterialData * md, const arma::mat &) const;
  double active_strain_energy(int iel, MaterialData * md, const arma::mat & E) const;
  void deviatoric_stress(MaterialData * md, arma::mat & stress) const;
  void deviatoric_elastensor(MaterialData * md, Tensor4 & A) const;
  void piola2_stress(MaterialData * md, arma::mat & pk2) const;
  
};

#endif
