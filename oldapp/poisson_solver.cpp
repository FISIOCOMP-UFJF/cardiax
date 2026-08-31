#include <iostream>
#include <fstream>
#include <armadillo>
#include "pdes/poisson.hpp"
#include "pdes/laplace.hpp"
#include "util/command_line_args.h"

static char help[] = "Poisson Solver.\n\n";

void usage()
{
  std::cout << std::endl;
  std::cout << " Usage: poisson <options>\n";
  std::cout << " Options: " << std::endl;
  std::cout << "    -f    \t mesh file" << std::endl;
  std::cout << "    -p    \t parameter file" << std::endl;
  std::cout << "    -o    \t output file name" << std::endl;
  std::cout << "    -text \t write to text (instead of VTK)" << std::endl;
  std::cout << "    -type \t poisson or laplace" << std::endl;
  std::cout << std::endl;
  exit(0);
}

/** Some exact solutions for testing Poisson problem
       2d ->  u[i] = x*(x-1.)*y*(y-1.);
       3d ->  u(i) = x*(x-1.)*y*(y-1.)*z*(z-1.);
 */

void calc_exact_solution_2D(const Mesh & msh, arma::vec & u)
{
  double x,y;
  std::vector<arma::vec3> pts = msh.get_points();

  u.resize((int)pts.size());
  for(uint i=0; i<pts.size(); i++){
    x = pts[i](0);
    y = pts[i](1);
    u[i] = (1.0/(2.0*M_PI*M_PI))*(sin(M_PI*x)*sin(M_PI*y));
  }
}

void calc_exact_solution_3D(const Mesh & msh, arma::vec & u)
{
  double x,y,z;
  std::vector<arma::vec3> pts = msh.get_points();

  u.resize((int)pts.size());
  for(uint i=0; i<pts.size(); i++){
    x = pts[i](0);
    y = pts[i](1);
    z = pts[i](2);
    u(i) = (1.0/(3.0*M_PI*M_PI))*(sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
  }
}

int main(int argc, const char* argv[])
{
  CommandLineArgs::init(argc, argv);

  bool is_poisson = true;
  string filename, parfile, output, problem;

  problem     = CommandLineArgs::read("-type","poisson");
  filename    = CommandLineArgs::read("-f","null");
  parfile     = CommandLineArgs::read("-p","null");
  output      = CommandLineArgs::read("-o","output");

  if(file_extension(filename) == "xml")
    parfile = filename;

  if(filename == "null" || parfile == "null" || problem == "null")
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

    msg("Reading parameters file");
    psolver->config(parfile);
    
    msg("Running solver");
    psolver->run(filename);

    msg("Writing data file");
    psolver->write_data(output);

    /*
    if (is_poisson)
    {
      cout << "Post-processing for Poisson problem" << endl;
      arma::sp_mat M;
      arma::vec ue, uh, e;
      
      if(psolver->get_mesh().get_n_dim() == 2 )
        calc_exact_solution_2D(psolver->get_mesh(), ue);
      else
				calc_exact_solution_3D(psolver->get_mesh(), ue);

      uh = psolver->get_solution();
      e  = ue - uh;

      psolver->calc_mass_matrix(M);

      double enorm = FETools::calc_L2_norm(M, e);
      cout << "L2 Error norm: " << std::setprecision(6) << enorm << endl;      
    }
    */

  } // <--- Poisson's destructor will get called here

  cout << "Done" << endl;

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  return 0;
}
