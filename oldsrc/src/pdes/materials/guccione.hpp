#ifndef GUCCIONE_HPP
#define GUCCIONE_HPP

#include "hyperelastic_material.hpp"
#include "incompressible_material.hpp"

/**
 * Guccione et al model
 */

class Guccione : public IncompressibleMaterial
{
protected:

    //! Material constants (dimension of stress)
    std::vector<double> v_Cg;

    //! Material constants (dimensionless)
    std::vector<double> v_bf, v_bt, v_bfs;

    std::vector<int> map_AHA;

    //! Auxiliary
    const arma::mat33 I;

  const int getMaterial(int marker) {return map_AHA[marker];}

public:
  Guccione (const std::vector<std::vector<double>> & prm, int num_materials, std::vector<int> & map_AHAmat) : IncompressibleMaterial(prm[0][5]), I(arma::eye<arma::mat>(3,3))
  {
    for (int ii=0; ii<num_materials; ii++)
    {
      v_Cg.push_back(prm[ii][0]);
      v_bf.push_back(prm[ii][1]);
      v_bt.push_back(prm[ii][2]);
      v_bfs.push_back(prm[ii][3]);
    }

    name = "Guccione";
    active_stress = prm[0][4];

    //parameters = prm[0];
    parameters.push_back(prm[0][0]);
    parameters.push_back(prm[0][1]);
    parameters.push_back(prm[0][2]);
    parameters.push_back(prm[0][3]);

    map_AHA = map_AHAmat;

  }

  Guccione (const std::vector<double> & prm) :
          IncompressibleMaterial(prm[5]), I(arma::eye<arma::mat>(3,3))
  {
    v_Cg.push_back(prm[0]);
    v_bf.push_back(prm[1]);
    v_bt.push_back(prm[2]);
    v_bfs.push_back(prm[3]);
    name = "Guccione";
    active_stress = prm[4];
    parameters = prm;
    for(int ii=0; ii<17; ii++)
      map_AHA.push_back(0);
  }



  double strain_energy(MaterialData * md, const arma::mat &) const;

    double active_strain_energy(MaterialData * md, const arma::mat & E) const;

    void piola2_stress(MaterialData * md, arma::mat & S) const;

    void cauchy_stress(MaterialData * md, arma::mat & sigma) const;

    void sp_elastensor(MaterialData * md, arma::mat & D) const;

    void mt_elastensor(MaterialData * md, arma::mat & D) const;

    void deviatoric_stress(MaterialData * md, arma::mat & stress) const;

    void deviatoric_elastensor(MaterialData * md, Tensor4 & A) const;

};


#endif /* GUCCIONE_HPP */
