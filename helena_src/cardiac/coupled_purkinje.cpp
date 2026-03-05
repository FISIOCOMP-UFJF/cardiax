#include "coupled_purkinje.hpp"

CoupledPurkinje::CoupledPurkinje() :
    pmj_region_spacing(200.0), pmj_coupling_resistance(1.0)
{
  tissue = new Monodomain();
  pktree = new MonodomainPurkinje();
}

CoupledPurkinje::~CoupledPurkinje()
{
  delete tissue;
  delete pktree;
}

void CoupledPurkinje::create_pmj_map()
{
  const double tol = 1.0e-8;
  const Mesh * tmesh = tissue->get_mesh_ptr();
  const PurkinjeMesh * pkmesh = dynamic_cast<const PurkinjeMesh*>( pktree->get_mesh_ptr() );

  std::vector<int>::const_iterator it;
  std::multimap<int,int>::const_iterator mit;
  const std::vector<int> pmj_pk = pkmesh->get_pmj_nodes();

  for(it = pmj_pk.begin(); it != pmj_pk.end(); ++it)
    cout << "PMJ Node = " << *it << endl;

  for(it = pmj_pk.begin(); it != pmj_pk.end(); ++it)
  {
    cout << "PMJ Node = " << *it << endl;
    const arma::vec3 x = pkmesh->get_point(*it);

    // now search the tissue mesh for the same point
    for(size_t i=0; i < tmesh->get_n_points(); i++)
    {
      const arma::vec3 mx = tmesh->get_point(i);
      const arma::vec3 dx = mx-x;
      if(arma::norm(dx,2) < tol) // found
        pmj_map.insert(std::pair<int,int>(*it,i));
    }
  }

  cout << "Initial pmj_map size: " << pmj_map.size() << endl;
  cout << "Adding neighboring nodes of the PMJ's region of influence" << endl;

  // copy all the tissue nodes to a vector
  std::vector<int> tmp_vec1; // purkinje
  std::vector<int> tmp_vec2; // tissue
  for(mit = pmj_map.begin(); mit != pmj_map.end(); ++mit)
  {
    tmp_vec1.push_back(mit->first);
    tmp_vec2.push_back(mit->second);
  }

  // now add the neighbors based on distance
  for(size_t j=0; j<tmp_vec1.size(); j++)
  {
    int p_node = tmp_vec1[j];
    int t_node = tmp_vec2[j];
    const arma::vec3 x = tmesh->get_point(t_node);
    for(size_t i=0; i < tmesh->get_n_points(); i++)
    {
      const arma::vec3 neighbour_x = tmesh->get_point(i);
      const arma::vec3 dx = neighbour_x - x;
      if(arma::norm(dx,2) < pmj_region_spacing)
        pmj_map.insert(std::pair<int,int>(p_node,i));
    }
  }

  // add the current values
  //for(mit = pmj_map.begin(); mit != pmj_map.end(); ++mit)
  //{
  //  pmj_myo_currents.insert(std::pair<int,double>(mit->second, stim_value));
  //}

  cout << "Final pmj_map size: " << pmj_map.size() << endl;


  // *********************************************
  // HARDCODED FOR TESTING!!!!
  // *********************************************
  /*
  //3 - 30401
  //4 - 10201
  pmj_map.insert(std::pair<int,int>(3,30401));

  pmj_map.insert(std::pair<int,int>(4,10201));
  pmj_map.insert(std::pair<int,int>(4,9793));
  pmj_map.insert(std::pair<int,int>(4,9797));
  pmj_map.insert(std::pair<int,int>(4,9801));
  pmj_map.insert(std::pair<int,int>(4,10197));
  pmj_map.insert(std::pair<int,int>(4,10205));
  pmj_map.insert(std::pair<int,int>(4,10601));
  pmj_map.insert(std::pair<int,int>(4,10605));
  pmj_map.insert(std::pair<int,int>(4,10609));

  pmj_myo_currents.insert(std::pair<int,double>(30401,1.0));

  pmj_myo_currents.insert(std::pair<int,double>(10201,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(9793,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(9797,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(9801,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(10197,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(10205,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(10601,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(10605,1.0));
  pmj_myo_currents.insert(std::pair<int,double>(10609,1.0));
  */


}

void CoupledPurkinje::inject_pmj_current(double time, double step)
{
  std::map<int,int>::const_iterator it;
  for(it = pmj_map.begin(); it != pmj_map.end(); ++it)
  {
    // compute v_pk - v_cm to impose I_pmj = 1/R_pmj (v_pk - v_cm)
    double v_pk = pktree->get_cells().get_state(it->first, 0);
    double v_cm = tissue->get_cells().get_state(it->second, 0);
    double i_pmj = (1.0/pmj_coupling_resistance) * (v_pk - v_cm);

    // impose current to myocardial tissue
    tissue->set_stimulus_value(it->second, i_pmj);
  }
}

void CoupledPurkinje::setup(std::string & mshname,
                            std::string & pkmshname,
                            std::string & cellmodel,
                            std::string & odesolver,
                            double dt, double T,
                            double pr, double pa)
{
  cout << "Setup of coupled (Purkinje/Tissue) cardiac problem" << endl;

  cout << "\nCardiac tissue" << endl;
  tissue->setup(mshname, cellmodel, odesolver, dt, T, pr, pa);
  tissue->init();

  cout << "\nPurkinje system" << endl;
  pktree->setup(pkmshname, cellmodel, odesolver, dt, T, pr, pa);
  pktree->init();

  cout << "\nConnecting PMJ nodes" << endl;
  create_pmj_map();

}

void CoupledPurkinje::solve()
{
  cout << "Simulating..." << endl;
  TimeParameters tip(tissue->get_time_parameters());

  tissue->initial_conditions();
  pktree->initial_conditions();

  // loop in time
  while (!tip.finished())
  {
    tip.increase_time();
    tip.show_time();

    inject_pmj_current(tip.time(), tip.get_dt());

    if(tip.time2print()) cout << "\n Tissue";
    tissue->advance();

    if(tip.time2print()) cout << " PKTree";
    pktree->advance();

  }

}