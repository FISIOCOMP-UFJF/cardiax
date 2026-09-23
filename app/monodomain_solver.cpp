#include <iostream>
#include <string>
#include <limits>
#include <filesystem>
#include "util/toml.hpp"
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "cardiac/monodomain.hpp"
#include "cardiac/monodomain_purkinje.hpp"
#include "cardiac/coupled_purkinje.hpp"

static char help[] = "Monodomain Cardiac Solver.\n\n";

using namespace std;
namespace fs = std::filesystem;

void usage()
{
  cout << endl;
  cout << " Usage: monodomain -f <config.toml> [OPTIONS]" << endl << endl;
  cout << "    -f       <config.toml>            TOML input file (required)" << endl;
  cout << endl;
  cout << "  Command-line options override the values from the TOML file:" << endl;
  cout << "    -dt      <dt>                     time step (ms)" << endl;
  cout << "    -t       <tend>                   total time of simulation (ms)" << endl;
  cout << "    -pr      <pr>                     print rate to output file" << endl;
  cout << "    -c       <cellmodel>              string that identifies the ionic model" << endl;
  cout << "    -m       <odesolver>              ODE solver (ExplicitEuler, Implicit, ...)" << endl;
  cout << "    -ep      <model>                  monodomain | purkinje | coupled_purkinje" << endl;
  cout << "    -fp      <pkmesh>                 Purkinje mesh" << endl;
  cout << "    -restore <restorefile>            Restore state file" << endl;
  cout << "    -save_state <interval>            Save state interval (ms)" << endl;
  cout << "    -num_threads <n>                  number of OpenMP threads" << endl;
  cout << endl;
  exit(0);
}

void print_banner()
{
  cout << "========================================" << endl;
  cout << "   M O N O D O M A I N   S O L V E R    " << endl;
  cout << "========================================" << endl;
}

namespace {

// Fetch a mandatory value from the TOML; abort with a clear message if absent.
template <typename T>
T require(const toml::table & cfg, std::string_view path)
{
  auto v = cfg.at_path(path).value<T>();
  if (!v)
  {
    cerr << " *** ERROR: missing mandatory config key: " << path << endl;
    exit(1);
  }
  return *v;
}

// CLI-override helpers. Each reads the flag with a sentinel default; if the
// returned value differs from the sentinel, the flag was passed on the command
// line and should overwrite the value already taken from the TOML file.
void override_str(std::string & value, const char * flag)
{
  static const std::string SENT = "\x01__unset__\x01";
  std::string got = CommandLineArgs::read(flag, SENT);
  if (got != SENT) value = got;
}

void override_double(double & value, const char * flag)
{
  const double SENT = std::numeric_limits<double>::lowest();
  double got = CommandLineArgs::read(flag, SENT);
  if (got != SENT) value = got;
}

void override_int(int & value, const char * flag)
{
  const int SENT = std::numeric_limits<int>::lowest();
  int got = CommandLineArgs::read(flag, SENT);
  if (got != SENT) value = got;
}

} // namespace

int main(int argc, const char *argv[])
{
  CommandLineArgs::init(argc, argv);

  // -f is mandatory: the TOML config file.
  std::string inputfile = CommandLineArgs::read("-f", "");
  if (inputfile.empty()) usage();

  // Parse the TOML file.
  toml::table cfg;
  try
  {
    cfg = toml::parse_file(inputfile);
  }
  catch (const toml::parse_error & err)
  {
    cerr << " *** ERROR parsing '" << inputfile << "':\n" << err << endl;
    return 1;
  }

  print_banner();

  // --------------------------------------------------------------
  //  Values from the TOML file (base configuration).
  //  Mandatory keys use require<>; optional keys use value_or.
  // --------------------------------------------------------------
  std::string mshname   = require<std::string>(cfg, "problem.mesh");
  std::string model     = cfg["problem"]["model"].value_or("monodomain");

  double dt = cfg["time"]["dt"].value_or(0.1);
  double T  = cfg["time"]["t_end"].value_or(10.0);
  double tp = cfg["time"]["dt_output"].value_or(1.0);

  std::string cellmodel = cfg["cell_model"]["name"].value_or("TT2");
  std::string odesolver = cfg["numerics"]["ode_solver"].value_or("ExplicitEuler");

  // Purkinje
  std::string pkmshname    = cfg["problem"]["purkinje_mesh"].value_or("emptymesh");
  
  // Restore
  std::string restfilename = cfg["restore"]["file"].value_or("");
  double checkpoint_interval = cfg["restore"]["save_interval"].value_or(-1.0);

  int num_threads = cfg["numerics"]["num_threads"].value_or(
                      omp_get_max_threads() > 0 ? omp_get_max_threads() : 1);

  
  // --------------------------------------------------------------
  //  Command-line overrides: any flag actually passed wins over the
  //  file. Same flags as the original solver.
  // --------------------------------------------------------------
  override_double(dt,                  "-dt");
  override_double(T,                   "-t");
  override_double(tp,                  "-pr");
  override_str   (mshname,             "-f_mesh");   // see note below
  override_str   (pkmshname,           "-fp");
  override_str   (cellmodel,           "-c");
  override_str   (odesolver,           "-m");
  override_str   (model,               "-ep");
  override_str   (restfilename,        "-restore");
  override_double(checkpoint_interval, "-save_state");
  override_int   (num_threads,         "-num_threads");

  std::string typefile = mshname + ".typ";

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

  if(model == "monodomain")
  {
    cout << "Monodomain solver\n";

    omp_set_num_threads(num_threads);
    std::cout << "Solving Monodomain problem with: " << num_threads << " OpenMP threads." << std::endl;
    {
      Monodomain monodomain;
      msg("Reading parameters file");
      monodomain.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
      monodomain.set_parameters(cfg);
      monodomain.init(restfilename != "");

      if(!file_exists(typefile))
        cout << "Cells: all cells are of the same type\n";
      else
        monodomain.setup_types(typefile);

      monodomain.initial_conditions();
      monodomain.set_checkpoint_interval(checkpoint_interval);
      if(restfilename != "")
      {
        monodomain.restore_checkpoint(restfilename);
      }

      monodomain.solve();
    }
  }
  else if(model == "purkinje")
  {
    cout << "Monodomain-Purkinje Solver\n";
    {
      MonodomainPurkinje mp;
      mp.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
      mp.init();
      mp.solve();
    }
  }
  else if(model == "coupled_purkinje")
  {
    cout << "Coupled tissue/Purkinje Monodomain Solver\n";
    {
      CoupledPurkinje cp;
      if(pkmshname == "emptymesh")
      {
        std::cerr << "Error: unknown -fp PurkinjeMesh" << endl;
        exit(1);
      }
      cp.setup(mshname, pkmshname, cellmodel, odesolver, dt, T, tp, tp);
      cp.solve();
    }
  }
  else
  {
    cout << "Cardiac PDE Model " << model << " does not exist." << endl;
  }

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}


// #include <iostream>
// #include <string>
// #include "util/util.hpp"
// #include "util/command_line_args.h"
// #include "cardiac/monodomain.hpp"
// #include "cardiac/monodomain_purkinje.hpp"
// #include "cardiac/coupled_purkinje.hpp"
// #include<filesystem>

// static char help[] = "Monodomain Cardiac Solver.\n\n";

// using namespace std;
// namespace fs = std::filesystem;

// void usage()
// {
//   cout << endl;
//   cout << " Usage: ";
//   cout << "monodomain [OPTIONS]" << endl << endl;
//   cout << "    -f     <meshbase>              prefix of the mesh file" << endl;
//   cout << "    -dt    <dt>                    time step (ms)" << endl;
//   cout << "    -t     <tend>                  total time of simulation (ms)" << endl;
//   cout << "    -pr    <pr>                    print rate to output file" << endl;
//   cout << "    -c     <cellmodel>             string that identifies the ionic model" << endl;
//   cout << "    -m     <odesolver>             ODE solver (ExplicitEuler, Implicit, ...)" << endl;
// 	cout << "    -ep    <model>                 monodomain or purkinje" << endl;
//   cout << "    -fp    <pkmesh>                Purkinje mesh" << endl;
//   cout << "  -restore <restorefile>           Restore state file" <<endl; 
//   cout << "-save_state <checkpoint_interval>  Save state interval (ms)"<<endl; 
//   cout << "-num_threads\t number of OpenMP threads " << endl;
//   cout << endl;
//   exit(0);
// }

// int main(int argc, const char *argv[])
// {
//   double dt, T, tp, checkpoint_interval;
//   string mshname, pkmshname, cellmodel, odesolver, typefile, restfilename;
// 	string model;
//   int num_threads = 1; 
//   if (argc < 7) usage();

//   // Parse command line options
//   CommandLineArgs::init(argc, argv);
//   dt = CommandLineArgs::read("-dt",0.1);
//   T  = CommandLineArgs::read("-t",10.0);
//   tp = CommandLineArgs::read("-pr",1.0);
//   mshname = CommandLineArgs::read("-f", "emptymesh");
//   pkmshname = CommandLineArgs::read("-fp", "emptymesh");
//   cellmodel = CommandLineArgs::read("-c","TT2"); 
//   odesolver = CommandLineArgs::read("-m","ExplicitEuler");
// 	model = CommandLineArgs::read("-ep","monodomain");
//   restfilename = CommandLineArgs::read("-restore", "");
//   checkpoint_interval = CommandLineArgs::read("-save_state", -1.0);
// 	typefile  = mshname + ".typ";
//   num_threads = CommandLineArgs::read("-num_threads", omp_get_max_threads() > 0 ? omp_get_max_threads() : 1);

//   // Start PETSc
//   PetscMPIInt rank;
//   PetscMPIInt size;
//   PetscErrorCode ierr;

//   ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
//   ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

// 	if(model == "monodomain")
// 	{
// 		// Start PDE solver for Monodomain model
// 		cout << "Monodomain solver\n";

//     omp_set_num_threads(num_threads); 
//     std::cout << "Solving Monodomain problem with: " << num_threads << " OpenMP threads." << std::endl;
// 		{
// 			Monodomain monodomain;
//       msg("Reading parameters file");
// 			monodomain.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
// 			monodomain.init(restfilename != "");
			
// 			if(!file_exists(typefile))
// 				cout << "Cells: all cells are of the same type\n";
// 			else    
// 				monodomain.setup_types(typefile);

// 			monodomain.initial_conditions();
//       monodomain.set_checkpoint_interval(checkpoint_interval);
//       if(restfilename != "")
//       {
//         monodomain.restore_checkpoint(restfilename); 
//       }

//       monodomain.solve();   
// 		}		
// 	}
// 	else if(model == "purkinje")
// 	{
// 		cout << "Monodomain-Purkinje Solver\n";
// 		{
// 			MonodomainPurkinje mp;
// 			mp.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
// 			mp.init();					
// 			mp.solve();   
// 		}		
// 	}
//   else if(model == "coupled_purkinje")
//   {
//     cout << "Coupled tissue/Purkinje Monodomain Solver\n";
//     {
//       CoupledPurkinje cp;
//       if(pkmshname == "emptymesh")
//       {
//         std::cerr << "Error: unknown -fp PurkinjeMesh" << endl;
//         exit(1);
//       }
//       cp.setup(mshname, pkmshname, cellmodel, odesolver, dt, T, tp, tp);
//       cp.solve();
//     }
//   }
// 	else
// 	{
// 		cout << "Cardiac PDE Model " << model << " does not exist." << endl;
// 	}
	
//   // End PDE solver
  
//   msg("Done.");
//   ierr = PetscFinalize(); CHKERRQ(ierr);
  
//   return 0;
// }