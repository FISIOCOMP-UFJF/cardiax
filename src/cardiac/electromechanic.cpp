#include <armadillo>
#include "util/command_line_args.h"
#include "electromechanic.hpp"
#include "util/pugixml.hpp"

void Electromechanic::config(const string & basename)
{
  filename = basename;
  read_options();
  setup_electrophysiology();
  setup_mechanics();
  setup_output();
}

void Electromechanic::read_options()
{
  t_final    = CommandLineArgs::read("-t",       1.0);
  dt_ephy    = CommandLineArgs::read("-dt_ephy", 0.001);
  dt_mech    = CommandLineArgs::read("-dt_mech", 0.1);
  dt_out     = CommandLineArgs::read("-dt_out",  1.0);
  cell_model = CommandLineArgs::read("-c", "ToRORdLand");
  odesolver  = CommandLineArgs::read("-m", "ExplicitEuler");   // make odesolver a member

  // number of electrical steps per interval; must be a whole number >= 1
  auto steps = [](double interval, double dt) {
    const long n = std::lround(interval / dt);
    return (n >= 1 && std::abs(n * dt - interval) <= 1e-9 * interval) ? static_cast<int>(n) : 0;
  };
  steps_per_mech = steps(dt_mech, dt_ephy);
  steps_per_out  = steps(dt_out,  dt_ephy);

  if (steps_per_mech == 0 || steps_per_out == 0 || steps_per_out % steps_per_mech != 0)
  {
    cerr << "Error: need dt_ephy | dt_mech | dt_out (each a multiple of the previous); got "
         << dt_ephy << ", " << dt_mech << ", " << dt_out << endl;
    exit(1);
  }

  cout << "Time step of electrics: " << dt_ephy << endl
       << "Time step of mechanics: " << dt_mech << endl
       << "Time step of outputs  : " << dt_out  << endl
       << "Total simulation time : " << t_final << endl;
}

void Electromechanic::setup_electrophysiology()
{
  if (ep_type == "mono")
    ephy = new Monodomain();
  else if (ep_type == "monodef")
    ephy = new MonodomainDeformation();
  else
  {
    cerr << "Error: unknown electrophysiology model '" << ep_type
         << "' (use mono or monodef)" << endl;
    exit(1);
  }
  cout << "Electrophysiology model: " << ep_type << endl;

  const double screen_rate = 1.0;          // time between "Time:" prints
  ephy->setup(filename, cell_model, odesolver, dt_ephy, t_final, screen_rate, screen_rate);
  ephy->set_write_output(false);           // only elas writes
  ephy->update_matrix(false);
  ephy->init();

  const CellModel & cm = ephy->get_cells().get_model();
  if (!cm.has_active_tension())
  {
    cerr << "Error: cell model '" << cell_model
         << "' has no active tension; it can't drive the mechanics." << endl;
    exit(1);
  }
}

void Electromechanic::setup_mechanics()
{
  elas.config(filename, filename);
  elas.set_output_step(false);
  elas.save_checkpoint(false);
  elas.init();

  Mesh & mm = elas.get_mesh();
  const Mesh & me = ephy->get_mesh();
  if (mm.get_n_points() != me.get_n_points() || mm.get_n_elements() != me.get_n_elements())
  {
    cerr << "Error: mechanical and electrical meshes differ" << endl;
    exit(1);
  }

  // the mechanics uses the electrophysiology fibres
  for (int e = 0; e < mm.get_n_elements(); e++)
  {
    arma::vec3 f0 = ephy->get_fiber(e);
    arma::vec3 s0 = ephy->get_trans(e);
    arma::vec3 n0 = ephy->get_normal(e);
    mm.set_element_fiber (e, f0);
    mm.set_element_trans (e, s0);
    mm.set_element_normal(e, n0);
  }
}

void Electromechanic::setup_output()
{
  const int n_points = elas.get_mesh().get_n_points();
  vm.zeros(n_points);
  ta.zeros(n_points);
  u_field.zeros(3 * n_points);

  // frames at t = 0, dt_out, ..., counted with the loop's own step numbers
  const int n_steps  = ephy->get_time_parameters().get_nsteps();
  const int n_frames = n_steps / steps_per_out;

  // fields the mechanics doesn't register itself; before setup_data_writer(),
  // which opens the file
  WriterHDF5 & w = elas.get_writer();
  w.add_field("vm",          FieldLoc::Point, FieldKind::Transient);
  w.add_field("ta",          FieldLoc::Point, FieldKind::Transient);
  w.add_field("long_strain", FieldLoc::Cell,  FieldKind::Transient);   // storeStress()
  w.add_field("circ_strain", FieldLoc::Cell,  FieldKind::Transient);
  w.add_field("rad_strain",  FieldLoc::Cell,  FieldKind::Transient);
  w.add_field("fiber",       FieldLoc::Cell,  FieldKind::Static, 3);

  elas.setup_data_writer(n_frames);

  w.write_static_field(ephy->direction_field(0).data(), "fiber");
}

// void Electromechanic::config(const string &basename)
// {
//   string parfile = basename;
//   string mshfile = basename;
//   string monofile = basename;
//   string odesolver;

//   filename = basename;

//   t_final = CommandLineArgs::read("-t", 1.0);
//   dt_out  = CommandLineArgs::read("-dt_out", 1.0);
//   dt_ephy = CommandLineArgs::read("-dt_ephy", 0.001);
//   dt_mech = CommandLineArgs::read("-dt_mech", 0.1);
//   cell_model = CommandLineArgs::read("-c", "ToRORdLand"); // using ToRORdLand as a standard value
//   odesolver = CommandLineArgs::read("-m", "ExplicitEuler");

//   // number of electrical steps per interval; must be a whole number >= 1
//   auto steps = [](double interval, double dt) {
//     const long n = std::lround(interval / dt);
//     return (n >= 1 && std::abs(n * dt - interval) <= 1e-9 * interval) ? static_cast<int>(n) : 0;
//   };

//   steps_per_mech = steps(dt_mech, dt_ephy);
//   steps_per_out  = steps(dt_out,  dt_ephy);
  
//   if (steps_per_mech == 0 || steps_per_out == 0 || steps_per_out % steps_per_mech != 0)
//   {
//     cerr << "Error: need dt_ephy | dt_mech | dt_out (each a multiple of the previous); got "
//          << dt_ephy << ", " << dt_mech << ", " << dt_out << endl;
//     exit(1);
//   }

//   cout << "Time step of mechanics: " << dt_mech << endl;
//   cout << "Time step of electrics: " << dt_ephy << endl;
//   cout << "Time step of outputs  : " << dt_out << endl;
//   cout << "Total simulation time : " << t_final << endl;

//   pugi::xml_document doc;
//   pugi::xml_parse_result result = doc.load_file(mshfile.c_str());
//   if (!result)
//   {
//     std::cout << "XML parsed with errors\n";
//     std::cout << "Error description: " << result.description() << "\n";
//     std::cout << "Error offset: " << result.offset;
//     std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
//   }

//   if (ep_type == "mono") {
//     ephy = new Monodomain();
//     cout << "Using the Monodomain electrophysiology model." << endl;
//   } else if (ep_type == "monodef") {
//     ephy = new MonodomainDeformation();
//     cout << "Using the MonodomainDeformation electrophysiology model." << endl;

//   } else {
//     cout << "Error: Unknown electrophysiology model '" << ep_type << "'!" << endl;
//     cout << "Available options: mono (Monodomain) or monodef (MonodomainDeformation)"<<endl; 
//     exit(1);
//   }

//   const int n_steps  = ephy->get_time_parameters().get_nsteps();   // EP steps in the run
//   const int n_frames = n_steps / steps_per_out + 1;                 // +1: frame at t = 0 (item 5)

//   ephy->setup(monofile, cell_model, odesolver, dt_ephy, t_final, 1.0, 1.0);
//   ephy->set_write_output(false); // only elas writes; vm goes through elas
//   ephy->update_matrix(false);
//   ephy->init();

//   CellModel & cm = ephy->get_cells().get_model();
//   if (!cm.has_active_tension())
//   {
//     cerr << "Error: cell model '" << cell_model
//          << "' has no active tension; it can't drive the mechanics." << endl;
//     exit(1);
//   }

//   elas.config(mshfile, parfile);
//   elas.set_output_step(false);
//   elas.save_checkpoint(false);
//   elas.init();
  
//   // Old
//   // elas.setup_data_writer( static_cast<int> (T / dt_mech) + 1 );

//   // New
//   // Fields this coupled solver writes through elas that the mechanics
//   // doesn't register itself. Must come before setup_data_writer(),
//   // which calls open().
//   WriterHDF5 & w = elas.get_writer();
//   w.add_field("vm",          FieldLoc::Point, FieldKind::Transient);
//   w.add_field("ta",          FieldLoc::Point, FieldKind::Transient);  // output_vtk(ii, vm, ...)
//   w.add_field("long_strain", FieldLoc::Cell,  FieldKind::Transient);  // storeStress(ii)
//   w.add_field("circ_strain", FieldLoc::Cell,  FieldKind::Transient);
//   w.add_field("rad_strain",  FieldLoc::Cell,  FieldKind::Transient);
//   w.add_field("fiber",       FieldLoc::Cell,  FieldKind::Static, 3);

//   elas.setup_data_writer(n_frames);

//   w.write_static_field(ephy->direction_field(0).data(), "fiber");

//   // Configure fibers for mechanical problem
//   for (int i = 0; i < elas.get_mesh().get_n_elements(); i++)
//   {
//     arma::vec3 f0 = ephy->get_fiber(i);
//     arma::vec3 s0 = ephy->get_trans(i);
//     arma::vec3 n0 = ephy->get_normal(i);
//     elas.get_mesh().set_element_fiber(i, f0);
//     elas.get_mesh().set_element_trans(i, s0);
//     elas.get_mesh().set_element_normal(i, n0);
//   }

//   // Configure fiber-sheet-normal directions
//   for (int i = 0; i < ephy->get_mesh().get_n_elements(); i++)
//   {
//     arma::mat33 *R = new arma::mat33();
//     arma::mat33 *R0 = new arma::mat33();

//     R->col(0) = ephy->get_fiber(i);
//     R->col(1) = ephy->get_trans(i);
//     R->col(2) = ephy->get_normal(i);

//     R0->col(0) = ephy->get_fiber(i);
//     R0->col(1) = ephy->get_trans(i);
//     R0->col(2) = ephy->get_normal(i);

//     vec_fib.push_back(R);
//     vec_fib0.push_back(R0);
//   }
  
//   int neln = ephy->get_mesh().get_nen();
//   int nelem = elas.get_mesh().get_n_elements();
//   int nnodes  = ephy->get_mesh().get_n_points(); 
//   for (int i = 0; i < nelem; i++) 
//     for (int j = 0; j < neln; j++)
//       vec_stress.push_back(new arma::mat33());

//   // Check if mechanical and electrical meshes are the same
//   assert(elas.get_mesh().get_n_points() == ephy->get_mesh().get_n_points());
//   assert(elas.get_mesh().get_n_elements() == ephy->get_mesh().get_n_elements());

//   ta.zeros(nnodes);
// }

void Electromechanic::solve()
{
  cout << "Solving coupled electromechanical problem" << endl;

  const int n_points = elas.get_mesh().get_n_points();
  std::string restart_file = CommandLineArgs::read("-restore", "");

  const std::vector<arma::mat33*> & vec_ftens = elas.get_vec_F();

  arma::vec vm(n_points, arma::fill::zeros);
  arma::vec u_field(3 * n_points, arma::fill::zeros);

  TimeParameters tip(ephy->get_time_parameters());
  elas.pre_solve();
  ephy->initial_conditions();
  tip.reset();

  int frame = 0;  // next output frame

  if (restart_file.empty())
  {
    // (5) frame 0: state at t = 0
    ephy->get_cells().get_var(0, vm);
    elas.get_displacements(u_field);
    elas.store_point_field(frame, vm, "vm");
    elas.store_point_field(frame, ta, "ta");
    elas.get_writer().write_displ_step(frame, u_field.memptr());
    frame += 1;
  }
  else
  {
    int restored_step;
    double restored_time;

    double* state_vars = ephy->get_cells().get_state_vars();
    elas.read_ep_checkpoint(restart_file, vm.memptr(), state_vars, restored_step, restored_time);

    elas.restore_checkpoint(restart_file);
    elas.lc.reset();

    tip.restore_state(restored_step, restored_time);

    TimeParameters& internal_tip = const_cast<TimeParameters&>(ephy->get_time_parameters());
    internal_tip.restore_state(restored_step, restored_time);

    frame = restored_step / steps_per_out + 1;      // (5) +1 for the t = 0 frame
  }

  //
  // time loop
  //
  while (!tip.finished())
  {
    tip.increase_time();
    const int it = tip.it();                        // (6) the only step counter

    if (tip.time2print())
      cout << "Time: " << tip.time() << endl;       // (7)

    // electrophysiology: one step
    timer.enter("Monodomain");
    if (ep_type == "monodef")
      static_cast<MonodomainDeformation*>(ephy)->advance(vec_ftens);
    else
      ephy->advance();

    ephy->get_cells().get_var(0, vm);
    ephy->get_cells().get_active_tension(ta);
    timer.leave();

    // mechanics: every steps_per_mech electrical steps
    if (it % steps_per_mech == 0)
    {
      timer.enter("Elasticity");
      elas.set_active_stress(ta);
      elas.solve();
      timer.leave();

      // output: every steps_per_out electrical steps
      if (it % steps_per_out == 0)
      {
        timer.enter("Writing");
        elas.get_displacements(u_field);
        elas.store_point_field(frame, vm, "vm");
        elas.store_point_field(frame, ta, "ta");
        elas.get_writer().write_displ_step(frame, u_field.memptr());
        elas.storeStress(frame);
        frame += 1;
        timer.leave();
      }

      // (4) checkpoint: test > 0 first
      if (checkpoint_interval > 0 && it % checkpoint_interval == 0)
      {
        timer.enter("Save Checkpoint");
        cout << "Saving coupled checkpoint at t = " << tip.time() << endl;
        elas.save_coupled_checkpoint(it, tip.time(), vm.memptr(),
                                     ephy->get_cells().get_state_vars(),
                                     ephy->get_cells().get_ode_size());
        timer.leave();
      }

      elas.reset();
    }
  }

  elas.timer.summary();
  ephy->timer.summary();
  timer.summary();
}




// void Electromechanic::solve()
// {
//   cout << "Solving coupled electromechanical problem" << endl;

//   int ii = 0;
//   int size = elas.get_mesh().get_n_points();
//   std::string restart_file = CommandLineArgs::read("-restore", "");
 
//   const std::vector<arma::mat33*> & vec_ftens = elas.get_vec_F();

//   arma::vec vm(size), u_field(3 * size);
//   vm.zeros();


//   TimeParameters tip(ephy->get_time_parameters());

//   elas.pre_solve();
//   ephy->initial_conditions(); 
  
//   tip.reset();

//   const int steps_per_mech = static_cast<int>(std::lround(dt_mech / dt_ephy)); // EP steps between mechanics solves
//   const int steps_per_out  = static_cast<int>(std::lround(dt_out  / dt_ephy)); // EP steps between output frames

//   if (!restart_file.empty())
//   {
//       int restored_step;
//       double restored_time;

//       double* state_vars = ephy->get_cells().get_state_vars();
//       elas.read_ep_checkpoint(restart_file, vm.memptr(), state_vars, restored_step, restored_time);

//       elas.restore_checkpoint(restart_file);
//       elas.lc.reset(); 

//       tip.restore_state(restored_step, restored_time);
//       i = restored_step;

//       TimeParameters& internal_tip = const_cast<TimeParameters&>(ephy->get_time_parameters());
//       internal_tip.restore_state(restored_step, restored_time);
      
//       ii = tip.it() / steps_per_out;
      
//   }

//   //
//   // time loop for solution
//   //
//   while (!tip.finished())
//   {
//     tip.increase_time();
//     if (tip.time2print()) cout << "Time: " << tip.time();
    
//     timer.enter("Monodomain");
//     // Update the active stress value
//     if(ep_type == "monodef")
//     {
//       MonodomainDeformation* monodef = static_cast<MonodomainDeformation*>(ephy);
//       monodef->advance(vec_ftens);
//     }
//     else
//     {
//       ephy->advance();
//     }

//     ephy->get_cells().get_var(0, vm);
//     ephy->get_cells().get_active_tension(ta);
//     timer.leave(); 

//     if (tip.it() % steps_per_mech == 0) 
//     {
//       timer.enter("Elasticity");
//       elas.set_active_stress(ta); 
//       elas.solve();
//       timer.leave(); 

//       if (tip.it() % steps_per_out == 0)
//       {
//         timer.enter("Writing");
//         elas.get_displacements(u_field);
//         // elas.output_vtk(ii, vm, u_field);
//         // elas.output_vtk(ii, ta, u_field);
//         elas.store_point_field(ii, vm, "vm");
//         elas.store_point_field(ii, ta, "ta");
//         elas.get_writer().write_displ_step(ii, u_field.memptr());
//         elas.storeStress(ii);
//         ii += 1;
//         timer.leave(); 
//       } 

//       // DEBUG
//       // arma::vec ta_eff;
//       // elas.effective_active_tension(ta, ta_eff);
//       // cout << " Ta: cell max = " << ta.max() << "   effective max (x ta_scale) = " << ta_eff.max() << endl;
//       // cout << " max |u| = " << arma::abs(u_field).max() << " um" << endl;

//       // t =   65.00 | EP lin:  1.0 its, max rnorm 0.0e+00 | Newton:  4 its, rnorm 4.8e-08 | max Vm 0.99  Ta 13.28  |u| 3736.0 um

//       if(tip.it() % checkpoint_interval == 0 && checkpoint_interval > 0)
//         {
//           timer.enter("Save Checkpoint");
//           cout << "Saving Coupled Checkpoint at t = " << tip.time() << " ms" << endl;
//           int num_state_vars =  ephy->get_cells().get_ode_size();
//           const double* state_vars = ephy->get_cells().get_state_vars();  
//           elas.save_coupled_checkpoint(tip.it(), tip.time(), vm.memptr(), state_vars, num_state_vars);
//           timer.leave(); 
//         }

//       elas.reset();
//     }  

//   }
  
//   elas.timer.summary();
//   ephy->timer.summary(); 
//   timer.summary();
// }
