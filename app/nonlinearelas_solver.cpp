#include <iostream>
#include <string>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "pdes/total_lagrangian.hpp"
//#include "pdes/total_lagrangian_snes.hpp"
#include "pdes/updated_lagrangian.hpp"

static char help[] = "Non-Linear Elasticity Finite Element Solver.\n\n";

using namespace std;

void usage()
{
  cout << "\n Usage: nonlinear_elasticity_solver [OPTIONS]" << endl << endl;
  cout << "    -m   \t mesh file" << endl;
  cout << "    -p   \t parameter file" << endl;
  cout << "    -o   \t output directory" << endl;
  cout << "    -s   \t solution method: ul or tl" << endl;
  cout << "         \t    tl : Total Lagrangian" << endl;
  cout << "         \t    ul : Updated Lagrangian" << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  std::string smethod, mshfile, parfile, outdir, extension;

  // Parse command line options
  if (argc <= 1) usage();

  CommandLineArgs::init(argc, argv);
  outdir  = CommandLineArgs::read("-o","output/");
  mshfile = CommandLineArgs::read("-m","null");
  parfile = CommandLineArgs::read("-p","null");
  smethod = CommandLineArgs::read("-s","tl");

  extension = file_extension(mshfile);
  if(extension == "xml")
    parfile = mshfile;

  if(extension != "xml" && !file_exists(parfile))
    print_error("elasticity_solver.cpp", "open mesh file", "parameters file not found");

  if (!file_exists(mshfile))
    print_error("nonlinearelas_solver.cpp", "open file", "mesh file not found");

  // Start PETSC
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  cout << "Nonlinear Elasticity FEM Solver" << endl;

  // Start PDE solver
  {
    NonlinearElasticity * esolver;

    if(smethod=="ul")
    {
      msg("Updated Lagrangian Formulation");
      esolver = new UpdatedLagrangian();
    }
    else if(smethod=="tl")
    {
      msg("Total Lagrangian Formulation");
      esolver = new TotalLagrangian();
    }
      /*
      else if(smethod=="tlsnes")
      {
        msg("Total Lagrangian Formulation SNES");
        esolver = new TotalLagrangianSNES();
      }
      */
    else
    {
      throw std::runtime_error("Unknown formulation. Please use TL or UL.");
    }

    // msg("Setting material and elasticity type");
    // msg("Setting boundary conditions");
    // esolver->config(mshfile, parfile);
    // esolver->set_output_step(true);

    msg("Solving problem");
    esolver->run(mshfile, parfile);

    delete esolver;
  }

  msg("Done.");

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
