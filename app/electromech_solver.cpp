#include "util/log.hpp"
#include "util/command_line_args.h"
#include "boost/filesystem.hpp"
//#include "cardiac/coupled_electromechanic.hpp"
#include "cardiac/cardiac_cycle.hpp"
#include "../src/util/command_line_args.h"

using namespace std;
using namespace boost::filesystem; 

static char help[] = "coupled monodomain-mechanics solver";

void usage()
{
  cout << endl;
  cout << " Usage: monomech [OPTIONS]" << endl << endl;
  cout << "    -f    <meshbase>   prefix of the mesh file" << endl;
  cout << "    -dt   <dt>         time step (ms)" << endl;
  cout << "    -t    <time>       total time of simulation (ms)" << endl;
  cout << "    -c    <cell>       name of the cell model (FHN,NP,LR1,TT2)" << endl; 
  cout << "    -m    <odesolver>  ODE solver (ExplicitEuler, ImplicitEuler)" << endl;
  cout << "    -cond <condtype>   conductivity type (0 1 ... 5)" << endl;
  cout << "    -ep   <EP model>   Monodomain or Bidomain (mono or bido)" << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  string basename, meshname, ep_model;
  int condtype;

  if (argc <= 1) usage();

  CommandLineArgs::init(argc, argv);
  basename = CommandLineArgs::read("-f","emptymesh");
  ep_model = CommandLineArgs::read("-ep","mono");
  condtype = CommandLineArgs::read("-cond",0);

  //meshname = basename + ".msh";

  //if (!file_exists(meshname)) error("mesh file not found");

  // check and clean output directory
  if ( exists("output") )
  {
    remove_all("output");
    assert(!exists("output"));
  }

  create_directory("output");
  create_directories("output/vm_text/");

  // start solution
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  // start PDE solver
  {
    CardiacElectromechanic model(ep_model);
    model.config(basename);
cout << "Config OK" << endl;
    model.ref().set_conductivity(condtype);
    cout << "Iniciando solve"<<endl;
    model.solve();
  }
  // end of PDE solver

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}

