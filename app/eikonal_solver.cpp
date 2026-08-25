#include <iostream>
#include <string>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "cardiac/eikonal.hpp"
#include "cardiac/coupled_purkinje.hpp"
#include<filesystem>

static char help[] = "Eikonal Cardiac Solver.\n\n";

using namespace std;
namespace fs = std::filesystem;

void usage()
{
  cout << endl;
  cout << " Usage: ";
  cout << "Eikonal [OPTIONS]" << endl << endl;
  cout << "    -f     <meshbase>              prefix of the mesh file" << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char *argv[])
{
  double dt = 0.1, T = 10.0, tp = 1.0; 
  string mshname, pkmshname, cellmodel, odesolver, typefile, restfilename;

  // Parse command line options
  CommandLineArgs::init(argc, argv);
  T  = CommandLineArgs::read("-t",10.0);
  mshname = CommandLineArgs::read("-f", "emptymesh");
  cellmodel = CommandLineArgs::read("-c","Kerkoff2003"); 
  odesolver = CommandLineArgs::read("-m","ExplicitEuler");

  typefile  = mshname + ".typ";

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  cout << "Eikonal solver\n";
  Eikonal* eikonal = new Eikonal();

  msg("Reading parameters file");

  eikonal->setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
  eikonal->init();
  eikonal->solve(mshname);   

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);
  
  return 0;
}
