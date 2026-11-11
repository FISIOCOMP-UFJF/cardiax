#include "util/log.hpp"
#include "util/command_line_args.h"
#include "cardiac/electromechanic.hpp"
#include "../src/util/command_line_args.h"
#include<filesystem>

using namespace std;
namespace fs = std::filesystem;

static char help[] = "coupled monodomain-mechanics solver";

void usage()
{
  cout << endl;
  cout << " Usage: electromech -f <mesh.xml> [OPTIONS] [PETSc options]" << endl << endl;

  cout << "  Input:" << endl;
  cout << "    -f           <file>      mesh + mechanics + stimulus (XML)" << endl;
  cout << "    -ep          <model>     electrophysiology: mono | monodef" << endl;
  cout << endl;

  cout << "  Time (each step a multiple of the previous: dt_ephy | dt_mech | dt_out):" << endl;
  cout << "    -t           <time>      total simulation time (ms)" << endl;
  cout << "    -dt_ephy     <step>      electrophysiology time step (ms)" << endl;
  cout << "    -dt_mech     <step>      mechanics: solve every dt_mech (ms)" << endl;
  cout << "    -dt_out      <step>      output: one frame every dt_out (ms)" << endl;
  cout << endl;

  cout << "  Cell model:" << endl;
  cout << "    -c           <model>     cell model with active tension (NP, MNP, MS, TT2Ta, ToRORdLand, ...)" << endl;
  cout << "    -m           <method>    ODE solver (ExplicitEuler, ImplicitEuler, RungeKutta4)" << endl;
  cout << endl;

  cout << "  Checkpoint / restart:" << endl;
  cout << "    -save_state  <interval>  save a coupled checkpoint every <interval> ms" << endl;
  cout << "    -restore     <file>      restart from a coupled checkpoint" << endl;
  cout << endl;

  cout << "  Parallel:" << endl;
  cout << "    -num_threads <n>         number of OpenMP threads" << endl;
  cout << endl;

  cout << "  Example:" << endl;
  cout << "    electromech -f em_cylinder_5x48x12.xml -ep mono -c NP -t 100 \\" << endl;
  cout << "                -dt_ephy 0.01 -dt_mech 0.1 -dt_out 1 -ksp_type preonly -pc_type lu" << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  string basename, meshname, ep_model;
  int condtype, checkpoint_rate, num_threads;

  if (argc <= 1) usage();

  CommandLineArgs::init(argc, argv);
  basename = CommandLineArgs::read("-f","emptymesh");
  ep_model = CommandLineArgs::read("-ep","mono");
  condtype = CommandLineArgs::read("-cond",0);
  checkpoint_rate = CommandLineArgs::read("-save_state", -1);
  num_threads = CommandLineArgs::read("-num_threads", omp_get_max_threads() > 0 ? omp_get_max_threads() : 1);

  // Old for GMSH
  //meshname = basename + ".msh";

  // New for XML
  meshname = basename;

  if (!file_exists(meshname)) show_error("mesh file not found");

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
    Electromechanic model(ep_model);
    model.config(basename);
    if(ep_model == "monodef")
    {
      MonodomainDeformation& monodef = static_cast<MonodomainDeformation&>(model.ref());
      monodef.set_conductivity(condtype);
    }
    model.set_checkpoint_interval(checkpoint_rate);
    model.solve();
  }
  // end of PDE solver

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}

