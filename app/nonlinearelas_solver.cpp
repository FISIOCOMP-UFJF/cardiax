#include <iostream>
#include <string>
#include <limits>
#include "util/toml.hpp"
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "pdes/total_lagrangian.hpp"
#include "pdes/updated_lagrangian.hpp"

static char help[] = "Non-Linear Elasticity Solver.\n\n";

void print_banner()
{
  cout << "============================================================" << endl;
  cout << "  N O N L I N E A R   E L A S T I C I T Y   S O L V E R     " << endl;
  cout << "============================================================" << endl;
}

void usage()
{
  cout << "\n Usage: nonlinear_elasticity -f <config.toml> [OPTIONS]" << endl << endl;
  cout << "    -f          <config.toml>   TOML input file (required)" << endl;
  cout << endl;
  cout << "  Command-line options override the values from the TOML file:" << endl;
  cout << "    -o          <outdir>        output directory" << endl;
  cout << "    -s          <ul|tl>         solution method (ul=Updated, tl=Total Lagrangian)" << endl;
  cout << "    -restore    <restorefile>   restore state file" << endl;
  cout << "    -save_state <0|1>           save checkpoint" << endl;
  cout << "    -num_threads <n>            number of OpenMP threads" << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char* argv[])
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

  // values from the TOML (base configuration) 
  auto mshfile_opt = cfg.at_path("problem.mesh").value<std::string>();
  if (!mshfile_opt)
  {
    cerr << " *** ERROR: missing mandatory config key: problem.mesh" << endl;
    return 1;
  }
  std::string mshfile = *mshfile_opt;
  std::string smethod = cfg["solver"]["method"].value_or("ul");      // ul | tl
  std::string outdir  = cfg["output"]["directory"].value_or("output/");
  std::string restorefile = cfg["restore"]["file"].value_or("");
  int save_flag = cfg["restore"]["save_checkpoint"].value_or(0);
  int num_threads = cfg["solver"]["num_threads"].value_or(
                      omp_get_max_threads() > 0 ? omp_get_max_threads() : 1);

  // command-line overrides: any flag passed wins over the file
  override_str(mshfile,     "-f_mesh");   // see note
  override_str(smethod,     "-s");
  override_str(outdir,      "-o");
  override_str(restorefile, "-restore");
  override_int(save_flag,   "-save_state");
  override_int(num_threads, "-num_threads");

  bool save_checkpoint = (save_flag != 0);

  if (!file_exists(mshfile))
    print_error("nonlinearelas_solver.cpp", "open file", "mesh file not found");

  //
  // Start PETSc
  //
  PetscMPIInt rank, size;
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = PetscOptionsClearValue(NULL, "-f"); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

  apply_petsc_options_from_toml(cfg);   // TOML [solver.petsc] -> PETSc options DB

  cout << "Nonlinear Elasticity FEM Solver" << endl;
  omp_set_num_threads(num_threads);
  std::cout << "Solving Mechanical problem with: " << num_threads << " OpenMP threads." << std::endl;

  {
    NonlinearElasticity * esolver;

    if(smethod == "ul")
    {
      msg("Updated Lagrangian Formulation");
      esolver = new UpdatedLagrangian();
      esolver->save_checkpoint(save_checkpoint);
      esolver->run(mshfile, cfg, restorefile); 
    }
    else if(smethod == "tl")
    {
      msg("Total Lagrangian Formulation");
      esolver = new TotalLagrangian();
      esolver->run(mshfile, cfg, "");
    }
    else
    {
      throw std::runtime_error("Unknown formulation. Please use TL or UL.");
    }

    delete esolver;
  }

  msg("Done.");
  
  ierr = PetscFinalize(); 
  CHKERRQ(ierr);

  return 0;
}
