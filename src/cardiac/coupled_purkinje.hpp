#ifndef CARDIAX_COUPLED_PURKINJE_HPP
#define CARDIAX_COUPLED_PURKINJE_HPP

#include "monodomain.hpp"
#include "monodomain_purkinje.hpp"

class CoupledPurkinje
{

public:

  //! Constructor
  CoupledPurkinje();

  //! Destructor
  ~CoupledPurkinje();

  //! Configure PMJ map data
  void create_pmj_map();

  //! Inject current from the Purkinje System to Tissue at PMJs
  void inject_pmj_current(double t, double dt);

  // //! Setup the coupled problem
  // void setup(std::string & mshname, std::string & pkmshname,
  //            std::string & cellmodel, std::string & odesolver,
  //            double dt, double T, double pr, double pa);

  //! Solve the coupled problem
  void solve();

  //! Setup the coupled problem (no mesh reading / assembly yet)
  void setup(std::string & mshname, std::string & pkmshname,
             std::string & cellmodel, std::string & odesolver,
             double dt, double T, double pr, double pa);

  //! Read meshes, assemble matrices, build PMJ map
  void init();

  //! Access to the tissue solver (e.g. to call set_parameters)
  Monodomain * get_monodomain() { return tissue; }

  MonodomainPurkinje * get_purkinje() { return pktree; }

private:

  //! Monodomain model for cardiac tissue
  Monodomain * tissue;

  //! Monodomain model for Purkinje system
  MonodomainPurkinje * pktree;

  //! PMJ map links pktree nodes to tissue nodes
  std::multimap<int,int> pmj_map;

  //! PMJ currents
  //std::map<int,double> pmj_myo_currents;

  //! PMJ volumetric region of influence
  double pmj_region_spacing;

  //! Coupling PMJ resistance
  double pmj_coupling_resistance;

};


#endif //CARDIAX_COUPLED_PURKINJE_HPP
