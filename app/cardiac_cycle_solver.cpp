#include "util/log.hpp"
#include "util/command_line_args.h"
#include <filesystem>
//#include "cardiac/coupled_electromechanic.hpp"
#include "cardiac/cardiac_cycle.hpp"
#include "../src/util/command_line_args.h"

using namespace std;
namespace fs = std::filesystem;


static char help[] = "coupled monodomain-mechanics solver";

void usage()
{
  cout << endl;
  cout << " Usage: cardiac_cycle -f <mesh.xml> [OPTIONS]" << endl << endl;
  cout << " Units everywhere: time ms, length um, pressure Pa, volume mL." << endl << endl;
  cout << " General" << endl;
  cout << "    -f    <mesh.xml>   mesh + material + EP setup (one file)" << endl;
  cout << "    -c    <cell>       cell model: Kerkoff2003 | ToRORdLand" << endl;
  cout << "    -m    <odesolver>  ODE solver (ExplicitEuler, ImplicitEuler)" << endl;
  cout << "    -ep   <EP model>   mono | bido" << endl;
  cout << "    -cond <condtype>   conductivity type (0 1 ... 5)" << endl;
  cout << "    -num_threads <n>   OpenMP threads" << endl;
  cout << "    -amgx <file>       AMGX config (./configs/CG_DILU.json)" << endl;
  cout << endl << " Circulation model (default: legacy pressure curve)" << endl;
  cout << "    -circ <0|1>        closed-loop 0D circulation (Regazzoni 2022)" << endl;
  cout << "    -pw   <0|1>        piecewise PV loop (Wang et al. 2025); exclusive with -circ" << endl;
  cout << "    -beats <n>         beats to run (-circ / -pw)" << endl;
  cout << "    -rr   <ms>         cardiac period (-circ / -pw; default 1000)" << endl;
  cout << endl << " Time" << endl;
  cout << "    -dt   <ms>         EP step (default: Kerkoff2003 1.0, ToRORdLand 0.02)" << endl;
  cout << "    -dt_mech <ms>      coupling step (default 1.0); the EP runs dt_mech/dt" << endl;
  cout << "                       steps per coupling step. Ignored by the legacy model." << endl;
  cout << "    -t    <ms>         with -circ/-pw: MAXIMUM time (reduces -beats);" << endl;
  cout << "                       legacy: duration of the run" << endl;
  cout << "    -bcl  <ms>         cell warm-up BCL (default: the cardiac period)" << endl;
  cout << endl << " Output and diagnostics" << endl;
  cout << "    -prc  <n>          save the mesh every n steps (legacy 1, others 8)" << endl;
  cout << "    -logc <n>          detailed log every n steps (25)" << endl;
  cout << "    -dumpflags <0|1>   print every flag read, [cli]/[default] (1)" << endl;
  cout << endl << " Per-model flags (-circ, -pw, cell model, mechano-electric feedback):" << endl;
  cout << " see the flag table, or the block \"flags desta rodada\" printed at the" << endl;
  cout << " end of the configuration." << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  string basename, meshname, ep_model;
  int condtype, num_threads;

  if (argc <= 1) usage();

  CommandLineArgs::init(argc, argv);
  if (CommandLineArgs::wants_help()) usage();
  basename = CommandLineArgs::read("-f","emptymesh");
  ep_model = CommandLineArgs::read("-ep","mono");
  condtype = CommandLineArgs::read("-cond",0);
  num_threads = CommandLineArgs::read("-num_threads", omp_get_max_threads() > 0 ? omp_get_max_threads() : 1);


  // check and clean output directory
 if ( fs::exists("output") ) 
  {
    fs::remove_all("output");
    assert(!fs::exists("output"));
  }
  
  fs::create_directory("output");
  fs::create_directories("output/vm_text/");

  // start solution
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  // start PDE solver
  omp_set_num_threads(num_threads); 
  std::cout << "Solving Mechanical problem with: " << num_threads << " OpenMP threads." << std::endl;

  {
    CardiacElectromechanic model(ep_model);
    model.config(basename);
    cout << "Config OK" << endl;
    model.ref().set_conductivity(condtype);
    // O modelo de circulacao (legado, -circ 1 ou -pw 1) e escolhido em config()
    cout << "Iniciando solve" << endl;
    model.solve();
  }
  // end of PDE solver

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}