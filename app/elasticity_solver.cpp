#include <iostream>
#include <iterator>
#include <fstream>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "util/toml.hpp"
#include "pdes/elasticity.hpp"

static char help[] = "Linear Elasticity Solver.\n\n";

using namespace std;

void print_banner()
{
  std::cout <<
  "==============================================================\n"
  "   L I N E A R   E L A S T I C I T Y   S O L V E R\n"
  "==============================================================\n";
}

void usage()
{
  std::cout << std::endl;
  std::cout << " Usage: elasticity -f <config.toml> [options]\n";
  std::cout << " Options:\n";
  std::cout << "    -f    \t configuration file (.toml)\n";
  std::cout << "    -m    \t mesh file (.xml/.msh) — overrides 'mesh' in the TOML\n";
  std::cout << "    -o    \t output file name\n";
  std::cout << std::endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  CommandLineArgs::init(argc, argv);

  if(argc <= 1) usage();

  std::string cfgfile = CommandLineArgs::read("-f","null");
  std::string output  = CommandLineArgs::read("-o","output");

  if(cfgfile == "null" || !file_exists(cfgfile))
    print_error("elasticity_solver.cpp", "config", "TOML config file not found (pass -f <config.toml>)");

  // Parse the TOML config
  toml::table cfg;
  try {
    cfg = toml::parse_file(cfgfile);
  } catch (const toml::parse_error & err) {
    std::cerr << "Failed to parse '" << cfgfile << "':\n" << err << "\n";
    return 1;
  }

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

  print_banner();

  // Start PDE solver
  {
    Elasticity esolver;
    msg("Setting parameters");
    esolver.set_parameters(cfg);
    msg("Solving problem");
    esolver.run();           
    msg("Writing data file");
    esolver.write_data( cfg["output"]["filename"].value_or(std::string("output")) );
  }

  msg("Done");

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}

// #include <iostream>
// #include <iterator>
// #include <fstream>
// #include "util/util.hpp"
// #include "util/command_line_args.h"
// #include "pdes/elasticity.hpp"

// static char help[] = "Linear Elasticity Solver.\n\n";

// using namespace std;

// void usage()
// {
//   std::cout << std::endl;
//   std::cout << " Usage: elasticity <options>\n";
//   std::cout << " Options: " << std::endl;
//   std::cout << "    -f    \t mesh file" << std::endl;
//   std::cout << "    -p    \t parameter file" << std::endl;
//   std::cout << "    -o    \t output file name" << std::endl;
//   std::cout << std::endl;
//   exit(0);
// }

// int main(int argc, const char* argv[])
// {
//   CommandLineArgs::init(argc, argv);

//   std::string filename, parfile, output, extension;
//   filename  = CommandLineArgs::read("-f","null");
//   parfile   = CommandLineArgs::read("-p","null");
//   output    = CommandLineArgs::read("-o","output");
//   extension = file_extension(filename);

//   if(argc <= 1) usage();

//   if(extension == "xml")
//     parfile = filename;

//   if(extension != "xml" && !file_exists(parfile))
//     print_error("elasticity_solver.cpp", "open mesh file", "mesh file not found");

//   // Start PETSC
//   PetscMPIInt rank;
//   PetscMPIInt size;
//   PetscErrorCode ierr;

//   ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
//   ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);  

//   cout << "Linear Elasticity Solver" << endl;

//   // Start PDE solver
//   {
//     Elasticity esolver;
//     msg("Setting boundary conditions");
//     esolver.config(parfile);
//     msg("Solving problem");
//     esolver.run(filename);
//     msg("Writing data file");
//     esolver.write_data(output);
//   }

//   msg("Done");

//   ierr = PetscFinalize();
//   CHKERRQ(ierr);

//   return 0;
// }
