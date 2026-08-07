#include <iostream>
#include <fstream>
#include <armadillo>
#include "fem/projection.hpp"
#include "util/command_line_args.h"

static char help[] = "L2 Projection Solver.\n\n";

void usage()
{
  std::cout << std::endl;
  std::cout << " Usage: poisson <options>\n";
  std::cout << " Options: " << std::endl;
  std::cout << "    -m    \t mesh file" << std::endl;
  std::cout << "    -p    \t parameter file" << std::endl;
  std::cout << "    -o    \t output file name" << std::endl;
  std::cout << std::endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  CommandLineArgs::init(argc, argv);

  bool output_text;
  string filename, output;

  filename    = CommandLineArgs::read("-m","null");
  output      = CommandLineArgs::read("-o","output");
  output_text = CommandLineArgs::read("-text",0);

  if(filename == "null")
  {
    cout << "OI" << endl;
    usage();
  }

  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  // Start PETSC
  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  title("L2 Projection Solver");

  // Start PDE solver
  {
    Projection *psolver = new Projection();

    msg("Running");
    psolver->run(filename);
    msg("Writing data file");
    psolver->write_data(output, output_text);

  } // <--- L2 Projection's destructor will get called here

  msg("Done\n");

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}
