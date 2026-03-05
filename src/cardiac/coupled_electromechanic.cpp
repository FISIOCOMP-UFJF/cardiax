#include <armadillo>
#include "util/command_line_args.h"
#include "coupled_electromechanic.hpp"

Electromechanic::Electromechanic(const std::string & epmodel) : dt_mech(0.1), timer()
{
  // TODO: change from EPHY to EP
  // create electrophysiology model
  //if (epmodel == "mono")
  //  ep = new MonodomainDeformation();
  //else
  //  ep = new BidomainDeformation();
}

Electromechanic::~Electromechanic()
{
  // freeing memory
  //delete ep;

  for (uint i = 0; i < vec_stress.size(); i++)
    delete vec_stress[i];
  for (uint i = 0; i < vec_fib.size(); i++)
    delete vec_fib[i];
  for (uint i = 0; i < vec_fib0.size(); i++)
    delete vec_fib0[i];
}

void Electromechanic::config(const string & basename)
{
  double T, dt;
  string monofile = basename;
  string parfile = basename + "";
  string mshfile = basename + "";
  string typfile = basename + "";
  string cellmodel, odesolver;

  cout << endl << "Starting coupled electromechanics problem" << endl;

#ifdef DEBUG
  cout << "*** DEBUG MODE ON ***" << endl;
#endif
  T = CommandLineArgs::read("-t", 10.0);
  dt = CommandLineArgs::read("-dt", 0.1);
  dt_mech = CommandLineArgs::read("-dt_mech", 0.1);
  cellmodel = CommandLineArgs::read("-c", "NP");
  odesolver = CommandLineArgs::read("-m", "ExplicitEuler");

  //double TaConstante = CommandLineArgs::read("-ta",50);

  double pr = 1000; // int(0.1/dt);// dt_mech / dt;
  double pa = 1000; // int(0.1/dt);// dt_mech / dt;

  ephy.setup(monofile, cellmodel, odesolver, dt, T, pr, pa);

  ephy.update_matrix(false);

  ephy.init();

  ephy.setup_types(typfile);
  elas.config(mshfile, parfile);

  elas.set_output_step(false);
  elas.init();
  elas.setup_data_writer(int(T/1)); // save for each 0.1 ms


  // Configure fibers for mechanical problem
  for (int i = 0; i < elas.get_mesh().get_n_elements(); i++)
  {
    arma::vec3 f0 = ephy.get_fiber(i);
    arma::vec3 s0 = ephy.get_trans(i);
    arma::vec3 n0 = ephy.get_normal(i);
    elas.get_mesh().set_element_fiber(i, f0);
    elas.get_mesh().set_element_trans(i, s0);
    elas.get_mesh().set_element_normal(i, n0);
  }

  // Configure fiber-sheet-normal directions
  for (int i = 0; i < ephy.get_mesh().get_n_elements(); i++)
  {
    arma::mat33 * R = new arma::mat33();
    arma::mat33 * R0 = new arma::mat33();

    R->col(0) = ephy.get_fiber(i);
    R->col(1) = ephy.get_trans(i);
    R->col(2) = ephy.get_normal(i);

    R0->col(0) = ephy.get_fiber(i);
    R0->col(1) = ephy.get_trans(i);
    R0->col(2) = ephy.get_normal(i);

    vec_fib.push_back(R);
    vec_fib0.push_back(R0);
  }

  // Check if mechanical and electrical meshes are the same
  assert(elas.get_mesh().get_n_points() == ephy.get_mesh().get_n_points());
  assert(elas.get_mesh().get_n_elements() == ephy.get_mesh().get_n_elements());
}

void Electromechanic::solve()
{
  cout << "Solving coupled electromechanical problem" << endl;

  int nstep;
  int cont  = 0;
  int size  = ephy.get_mesh().get_n_points();
  int neln  = ephy.get_mesh().get_nen();
  int nelem = elas.get_mesh().get_n_elements();

  arma::vec vm(size), ta(size), ta0(size), dta(size);
  arma::vec vm2(size), dta2(size);
  arma::vec u_field(3*size);
  arma::mat xmat;

  // name of the cell model
  double dt = CommandLineArgs::read("-dt", 0.1);
  string cmname = CommandLineArgs::read("-c", "NP");

  // configure number of steps until mechanical solve
  nstep = dt_mech / dt;
  cout << "Time step of mechanics: " << dt_mech << endl;
  cout << "Number of steps until mechanical solve: " << nstep << endl;

  // get reference to the vector of tensor F in elasticity problem
  const std::vector<arma::mat33*> & vec_ftens = elas.get_vec_F();

  TimeParameters tip(ephy.get_time_parameters());

  // configure vector of displacements / stresses
  u_field.zeros();
  for (int i = 0; i < nelem; i++)
    for (int j = 0; j < neln; j++)
      vec_stress.push_back(new arma::mat33());


  // initialization, initial conditions and Vm
  ephy.initial_conditions();
  ephy.get_cells().get_var(0, vm);
  //ephy.write_data(vm, u_field, "vm", cont);
  //ephy.write_data_text(vm, cont);
  //elas.output_vtk(0, cont, "vm", vm, vec_fib0);
  ta0.zeros();
  // Used to get the initial value of Ta
  //ephy.get_cells().get_monitored_values(0, ta0);

  elas.output_vtk(0, vm, u_field);

  // loop in time
  while (!tip.finished())
  {
    tip.increase_time();
    cout << "\nTime: " << tip.time() << endl;

    // ------------------------------------------------------------------------
    // Solve Monodomain
    // ------------------------------------------------------------------------

    timer.enter("Monodomain");
    ephy.advance(vec_ftens);
    timer.leave();

    // ------------------------------------------------------------------------
    // Solve Mechanics
    // ------------------------------------------------------------------------
    ephy.get_cells().get_var(0, vm);
    
    if ((tip.it() % nstep) == 0)
    {
      timer.enter("Mechanics");

      if (cmname == "NP")         ephy.get_cells().get_var(2, ta);
      else if (cmname == "MINI")   ephy.get_cells().get_var(5, ta);
      else if (cmname == "MMSilva")   ephy.get_cells().get_var(5, ta);
      else if (cmname == "MNP")   ephy.get_cells().get_var(2, ta);
      else if (cmname == "TT2Ta") ephy.get_cells().get_var(19, ta);
      else if (cmname == "MS")
      {
        ephy.get_cells().get_var(2, ta);
        ta *= 13.7;
      }
      else if (cmname == "RiceTT2")
      {
        ephy.get_cells().get_monitored_values(0, ta);
        ta = 300.0 * ta;
      }
      else if (cmname == "RiceORd")
      {
        ephy.get_cells().get_monitored_values(0, ta);
        ta = 10.0 * ta;
      }
      else if (cmname == "ToRORdLand")
      {
        ephy.get_cells().get_var(49, ta);
        ta = ta/25.;
      }
    
      // Assemble initial stress and solve
      dta = ta - ta0;
      ta0 = ta;

      /*
      // BIDOMAIN - RENUMBER DOFs
      for(uint k=0; k<ephy.bidomap.size(); k++)
        dta2(k) = dta(ephy.bidomap[k]);
      */

      // UL
      //elas.assemble_active(dta, vec_stress, vec_fib);

      // TL
      
      elas.assemble_active(dta, vec_stress, vec_fib0);
      elas.set_Ta(ta, dta);
      


      cout<<"antes do solve"<<endl;
      elas.solve(); //PROBLEMA!
      cout<<"depois do solve"<<endl; 


      cout<<"AQUI!"<<endl; 
      return; 
      
      // This should be used for the case where we update the coordinates
      // instead of using div( F^{-1} D F^{-T}) in the monodomain/bidomain eqs.
      elas.get_displacements(xmat);
      ephy.update_coords(xmat);

      elas.reset();

      // ----------------------------------------------------------------------
      // Output information
      // ----------------------------------------------------------------------

      
      elas.update_vectors(vec_fib0, vec_fib);
      
      //ephy.write_data(vm, u_field, "vm", cont);
      cout << "Writing: " << int(1/dt) << endl;
      if ((tip.it() % int(1/dt)) == 0) // save for each 0.1 ms
      {
        cont += 1;
        elas.get_displacements(u_field);
        elas.output_vtk(cont, vm, u_field);
      }
      

      /*
      // BIDOMAIN
      for(uint k=0; k<ephy.bidomap.size(); k++)
        vm2(k) = vm(ephy.bidomap[k]);
      elas.output_vtk(0, cont, "vm", vm2, vec_fib);
      */

      ephy.update_matrix(true);
      timer.leave();
    }
    //cont += 1;
    //ephy.write_data(vm, u_field, "vm", cont);
    
    
  }

  elas.timer.summary();
  ephy.timer.summary();
  timer.summary();
}
