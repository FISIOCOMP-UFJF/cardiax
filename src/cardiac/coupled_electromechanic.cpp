#include <armadillo>
#include "util/command_line_args.h"
#include "coupled_electromechanic.hpp"
#include "util/pugixml.hpp"


Electromechanic::Electromechanic(const std::string &epmodel) : dt_mech(1.0), timer(), ep_type(epmodel)
{
  // ~~
}

Electromechanic::~Electromechanic()
{
  // ~~
}

void Electromechanic::config(const string &basename)
{
  double T, dt;
  string parfile = basename;
  string mshfile = basename;
  string monofile = basename;
  string odesolver;
  filename = basename;

  cout << endl
       << "Starting coupled electromechanics problem" << endl;

  T = CommandLineArgs::read("-t", 1.0);
  dt = CommandLineArgs::read("-dt", 0.001);
  dt_mech = CommandLineArgs::read("-dt_mech", 0.1);
  cell_model = CommandLineArgs::read("-c", "ToRORdLand"); // using ToRORdLand as a standard value
  odesolver = CommandLineArgs::read("-m", "ExplicitEuler");

  double pr = 1.0;
  double pa = 1.0;

  ifstream inp;

  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());

  if (!result)
  {
    std::cout << "XML parsed with errors\n";
    std::cout << "Error description: " << result.description() << "\n";
    std::cout << "Error offset: " << result.offset;
    std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
  }

  cout << "total time: " << T << endl;
  cout << "Mechanical time step: " << dt_mech << endl;

  if (ep_type == "mono") {
    ephy = new Monodomain();
    cout << "Using the Monodomain electrophysiology model." << endl;
  } else if (ep_type == "monodef") {
    ephy = new MonodomainDeformation();
    cout << "Using the MonodomainDeformation electrophysiology model." << endl;

  } else {
    cout << "Error: Unknown electrophysiology model '" << ep_type << "'!" << endl;
    cout << "Available options: mono (Monodomain) or monodef (MonodomainDeformation)"<<endl; 
    exit(1);
  }

  ephy->setup(monofile, cell_model, odesolver, dt, T, pr, pa);
  ephy->update_matrix(false);
  ephy->init();

  elas.config(mshfile, parfile);
  elas.set_output_step(false);
  elas.init();
  elas.setup_data_writer( static_cast<int> (T / dt_mech) + 1 );

  // Configure fibers for mechanical problem
  for (int i = 0; i < elas.get_mesh().get_n_elements(); i++)
  {
    arma::vec3 f0 = ephy->get_fiber(i);
    arma::vec3 s0 = ephy->get_trans(i);
    arma::vec3 n0 = ephy->get_normal(i);
    elas.get_mesh().set_element_fiber(i, f0);
    elas.get_mesh().set_element_trans(i, s0);
    elas.get_mesh().set_element_normal(i, n0);
  }

  // Configure fiber-sheet-normal directions
  for (int i = 0; i < ephy->get_mesh().get_n_elements(); i++)
  {
    arma::mat33 *R = new arma::mat33();
    arma::mat33 *R0 = new arma::mat33();

    R->col(0) = ephy->get_fiber(i);
    R->col(1) = ephy->get_trans(i);
    R->col(2) = ephy->get_normal(i);

    R0->col(0) = ephy->get_fiber(i);
    R0->col(1) = ephy->get_trans(i);
    R0->col(2) = ephy->get_normal(i);

    vec_fib.push_back(R);
    vec_fib0.push_back(R0);
  }
  
  int neln = ephy->get_mesh().get_nen();
  int nelem = elas.get_mesh().get_n_elements();
  int nnodes  = ephy->get_mesh().get_n_points(); 
  for (int i = 0; i < nelem; i++) 
    for (int j = 0; j < neln; j++)
      vec_stress.push_back(new arma::mat33());

  // Check if mechanical and electrical meshes are the same
  assert(elas.get_mesh().get_n_points() == ephy->get_mesh().get_n_points());
  assert(elas.get_mesh().get_n_elements() == ephy->get_mesh().get_n_elements());

  ta.zeros(nnodes);

}


void Electromechanic::solve()
{
  cout << "Solving coupled electromechanical problem" << endl;

  int nstep;
  int i = 0, ii = 0;
  int size = elas.get_mesh().get_n_points();
  double dt = CommandLineArgs::read("-dt", 0.1);
 
  const std::vector<arma::mat33*> & vec_ftens = elas.get_vec_F();

  arma::vec vm(size), u_field(3 * size);
  vm.zeros();

  cout << "\nTime step of mechanics: " << dt_mech << endl;
  TimeParameters tip(ephy->get_time_parameters());
  elas.pre_solve();
  ephy->initial_conditions(); 
  
  tip.reset();

  while (!tip.finished())
  {
    tip.increase_time();
    cout << "Time: " << tip.time() << endl;

    i += 1;
    
    timer.enter("Monodomain");
    //Update the active stress value
    if(ep_type == "monodef")
    {
      MonodomainDeformation* monodef = static_cast<MonodomainDeformation*>(ephy);
      monodef->advance(vec_ftens);
    }
    else
    {
      ephy->advance();
    }

    ephy->get_cells().get_var(0, vm);
    timer.leave(); 

    // TODO: colocar um wrapper aqui
    if (cell_model == "NP")         ephy->get_cells().get_var(2, ta);
    else if (cell_model == "MINI")   ephy->get_cells().get_var(5, ta);
    else if (cell_model == "MMSilva")   ephy->get_cells().get_var(5, ta);
    else if (cell_model == "MNP")   ephy->get_cells().get_var(2, ta);
    else if (cell_model == "TT2Ta") ephy->get_cells().get_var(19, ta);
    else if (cell_model == "MS")
    {
      ephy->get_cells().get_var(2, ta);
      ta *= 13.7;
    }
    else if (cell_model == "RiceTT2")
    {
      ephy->get_cells().get_monitored_values(0, ta);
      ta = 300.0 * ta;
    }
    else if (cell_model == "RiceORd")
    {
      ephy->get_cells().get_monitored_values(0, ta);
      ta = 10.0 * ta;
    }
    else if (cell_model == "ToRORdLand")
    {
      ephy->get_cells().get_var(49, ta);
      ta = ta/25.;
    }else if (cell_model == "Kerkoff2003")
    {
      ephy->get_cells().get_monitored_values(0, ta); 
    }

    int save_freq = static_cast<int>(dt_mech / dt); //TODO: improve this

    if (i % save_freq == 0) 
    {
      timer.enter("Elasticity");
      elas.set_active_stress(ta); 
      elas.solve();
      elas.reset();
      timer.leave(); 

      timer.enter("Writing");
      elas.get_displacements(u_field);
      elas.output_vtk(ii, vm, u_field);
      elas.storeStress(ii);
      ii += 1;
      timer.leave(); 
    }  

    // LUCAS:
    // if checkpoint:
    // salvar tudo
  }
  
  elas.timer.summary();
  ephy->timer.summary(); 
  timer.summary();
}
