#include <iostream>
#include <fstream>
#include <armadillo>
#include "pdes/poisson.hpp"
#include "pdes/laplace.hpp"
#include "util/command_line_args.h"

#pragma push_macro("error")
#undef error
#include "util/toml.hpp"
#pragma pop_macro("error")

#pragma push_macro("error")
#undef error
#include "exprtk.hpp"
#pragma pop_macro("error")
#include "util/expr_function.hpp"

static char help[] = "Poisson Solver.\n\n";

void usage()
{
  std::cout << std::endl;
  std::cout << " Usage: poisson -f config.toml\n";
  std::cout << std::endl;
  exit(0);
}

void print_banner()
{
  cout << "========================================" << endl;
  cout << "     P O I S S O N   S O L V E R        " << endl;
  cout << "========================================" << endl;
}

/** Some exact solutions for testing Poisson problem
       2d ->  u[i] = x*(x-1.)*y*(y-1.);
       3d ->  u(i) = x*(x-1.)*y*(y-1.)*z*(z-1.);
 */

// void calc_exact_solution_2D(const Mesh & msh, arma::vec & u)
// {
//   double x,y;
//   std::vector<arma::vec3> pts = msh.get_points();

//   u.resize((int)pts.size());
//   for(uint i=0; i<pts.size(); i++){
//     x = pts[i](0);
//     y = pts[i](1);
//     u[i] = (1.0/(2.0*M_PI*M_PI))*(sin(M_PI*x)*sin(M_PI*y));
//   }
// }

// void calc_exact_solution_3D(const Mesh & msh, arma::vec & u)
// {
//   double x,y,z;
//   std::vector<arma::vec3> pts = msh.get_points();

//   u.resize((int)pts.size());
//   for(uint i=0; i<pts.size(); i++){
//     x = pts[i](0);
//     y = pts[i](1);
//     z = pts[i](2);
//     u(i) = (1.0/(3.0*M_PI*M_PI))*(sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
//   }
// }

int main(int argc, const char* argv[])
{
  CommandLineArgs::init(argc, argv);

  if (argc < 2 || CommandLineArgs::read("-h", false) || CommandLineArgs::read("--help", false))
    usage();

  print_banner();

  bool is_poisson = true;
  string configfile, filename, output, problem;

  // Single entry point: -f is the TOML config file.
  configfile = CommandLineArgs::read("-f", "config.toml");

  // parse TOML file
  toml::table tbl;
  try
  {
    tbl = toml::parse_file(configfile);
  }
  catch (const toml::parse_error & err)
  {
    std::cerr << " *** ERROR parsing '" << configfile << "':\n" << err << std::endl;
    return 1;
  }

  // [problem]
  problem  = tbl["problem"]["type"].value_or("poisson");   // poisson | laplace
  filename = tbl["problem"]["mesh"].value_or("null");      // XML for now, .h5 later
  output   = tbl["output"]["filename"].value_or("output");

  if(filename == "null" || problem == "null")
    usage();

  if(problem != "poisson")
    is_poisson = false;

  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  // Start PETSC
  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

  cout << "Poisson-Laplace Solver" << endl;

  // Start PDE solver
  {
    Poisson *psolver;

    if (problem == "poisson")
      psolver = new Poisson();
    else
      psolver = new Laplace();

    msg(("Reading config file: " + configfile).c_str());
    psolver->config(tbl);

    msg("Running solver");
    psolver->run(filename);

    msg("Writing data file");
    psolver->write_data(output);

    if (auto exact_str = tbl["validation"]["exact"].value<std::string>())
    {
      ExprScalarFunction exact(*exact_str);
      double enorm = psolver->calc_l2_error(exact);
      cout << " L2 error norm: " << std::setprecision(6) << enorm << endl;
    }
    else
    {
      cout << " No exact solution was provided in config; skipping L2 error." << endl;
    }

    // if (is_poisson)
    // {
    //   cout << "Post-processing for Poisson problem" << endl;
    //   arma::sp_mat M;
    //   arma::vec ue, uh, e;
    //   if(psolver->get_mesh().get_n_dim() == 2 )
    //     calc_exact_solution_2D(psolver->get_mesh(), ue);
    //   else
    //     calc_exact_solution_3D(psolver->get_mesh(), ue);
    //   uh = psolver->get_solution();
    //   e  = ue - uh;
    //   psolver->calc_mass_matrix(M);
    //   double enorm = FETools::calc_L2_norm(M, e);
    //   cout << "L2 Error norm: " << std::setprecision(6) << enorm << endl;
    // }

  } // <--- Poisson's destructor will get called here

  cout << "Done" << endl;

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}
