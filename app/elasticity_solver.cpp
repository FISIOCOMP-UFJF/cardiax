#include <iostream>
#include <iterator>
#include <fstream>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "pdes/elasticity.hpp"

static char help[] = "Linear Elasticity Solver.\n\n";

using namespace std;

void usage()
{
  std::cout << std::endl;
  std::cout << " Usage: elasticity <options>\n";
  std::cout << " Options: " << std::endl;
  std::cout << "    -f    \t mesh file" << std::endl;
  std::cout << "    -p    \t parameter file" << std::endl;
  std::cout << "    -o    \t output file name" << std::endl;
  std::cout << std::endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  CommandLineArgs::init(argc, argv);

  std::string filename, parfile, output, extension;
  filename  = CommandLineArgs::read("-f","null");
  parfile   = CommandLineArgs::read("-p","null");
  output    = CommandLineArgs::read("-o","output");
  extension = file_extension(filename);

  if(argc <= 1) usage();

  if(extension == "xml")
    parfile = filename;

  if(extension != "xml" && !file_exists(parfile))
    print_error("elasticity_solver.cpp", "open mesh file", "mesh file not found");

  // Start PETSC
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);  

  cout << "Linear Elasticity Solver" << endl;

  // Start PDE solver
  {
    Elasticity esolver;
    msg("Setting boundary conditions");
    esolver.config(parfile);
    msg("Solving problem");
    esolver.run(filename);
    msg("Writing data file");
    esolver.write_data(output);
  }

  msg("Done");

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}
