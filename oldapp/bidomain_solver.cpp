//
// FEM solver for Bidomain equations
// Bernardo M. Rocha
//

#include "cardiac/bidomain.hpp"
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;


static char help[] = "Bidomain Cardiac Solver.\n\n";

void usage()
{
  cout << endl;
  cout << " Usage: ";
  cout << "bidomain [OPTIONS]" << endl << endl;
  cout << "    -f   <meshbase>    prefix of the mesh file" << endl;
  cout << "    -dt  <dt>          time step (ms)" << endl;
  cout << "    -t   <tend>        total time of simulation (ms)" << endl;
  cout << "    -pr  <pr>          print rate to output file" << endl;
  cout << "    -pa  <pa>          print rate to output Vm to file for APD computation" << endl;
  cout << "    -c   <cellmodel>   string that identifies the ionic model" << endl;
  cout << "    -m   <odesolver>   ODE solver (ExplicitEuler, Implicit, ...)" << endl;
  cout << endl;
  exit(0);   
}

int main(int argc, const char *argv[])
{
  double dt, T, pr, pa;
  string filename, cellmodel, odesolver, typefile;

  if (argc < 7) usage();
 
  // parse command line options
  CommandLineArgs::init(argc, argv);
  dt = CommandLineArgs::read("-dt",0.1);
  T  = CommandLineArgs::read("-t",10.0);
  pr = CommandLineArgs::read("-pr",1.0);
  pa = CommandLineArgs::read("-pa",pr);
  filename  = CommandLineArgs::read("-f", "emptymesh");
  cellmodel = CommandLineArgs::read("-c", "TT2"); 
  odesolver = CommandLineArgs::read("-m", "ExplicitEuler");
  typefile  = filename + ".typ"; 
  
  // check and clean output directory
  if ( fs::exists("output") ) 
  {
    fs::remove_all("output");
    assert(!fs::exists("output"));
  }
  
  fs::create_directory("output");
  fs::create_directories("output/vm_text/");

  // start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);     
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  // start PDE solver
  {
      
    Bidomain bidomain;
    bidomain.setup(filename, cellmodel, odesolver, dt, T, pr, pa);
    bidomain.init();
    if(file_exists(typefile))
      bidomain.setup_types(typefile);
    bidomain.solve();

  } // <--- Bidomain's destructor will get called here

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}
