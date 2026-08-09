#include <iostream>
#include <string>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "cardiac/eikonal.hpp"
#include "cardiac/monodomain_purkinje.hpp"
#include "cardiac/monodomain.hpp"
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
  // cout << "    -dt    <dt>                    time step (ms)" << endl;
  cout << "    -t     <tend>                  total time of simulation (ms)" << endl;
  // cout << "    -pr    <pr>                    print rate to output file" << endl;
  // cout << "    -c     <cellmodel>             string that identifies the ionic model" << endl;
  // cout << "    -m     <odesolver>             ODE solver (ExplicitEuler, Implicit, ...)" << endl;
	// cout << "    -ep    <model>                 monodomain or purkinje" << endl;
  // cout << "    -fp    <pkmesh>                Purkinje mesh" << endl;
  // cout << "  -restore <restorefile>           Restore state file" <<endl; 
  // cout << "-save_state <checkpoint_interval>  Save state interval (ms)"<<endl; 
  cout << endl;
  exit(0);
}

int main(int argc, const char *argv[])
{
  double dt, T, tp, checkpoint_interval;
  string mshname, pkmshname, cellmodel, odesolver, typefile, restfilename;
	string model;

  // if (argc < 7) usage();

  // Parse command line options
  CommandLineArgs::init(argc, argv);
  // dt = CommandLineArgs::read("-dt",0.1);
  T  = CommandLineArgs::read("-t",10.0);
  // tp = CommandLineArgs::read("-pr",1.0);
  mshname = CommandLineArgs::read("-f", "emptymesh");
  // pkmshname = CommandLineArgs::read("-fp", "emptymesh");
  cellmodel = CommandLineArgs::read("-c","TT2"); 
  // odesolver = CommandLineArgs::read("-m","ExplicitEuler");
	// model = CommandLineArgs::read("-ep","monodomain");
  // restfilename = CommandLineArgs::read("-restore", "");
  // checkpoint_interval = CommandLineArgs::read("-save_state", -1.0);
	typefile  = mshname + ".typ";

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

		// Start PDE solver for Monodomain model
  cout << "Eikonal solver\n";

  Eikonal eikonal;
  msg("Reading parameters file");
  
  eikonal.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
  cout<<"depois do setup"<<endl;
  eikonal.init();
  cout<<"depois do init"<<endl; 
  
  // if(!file_exists(typefile))
  //   cout << "Cells: all cells are of the same type\n";
  // else    
  //   eikonal.setup_types(typefile);

  // eikonal.initial_conditions();
  cout<<"antes do solve"<<endl; 
  eikonal.solve(mshname);   
  cout<<"depois do solve"<<endl; 
	
	
  // End PDE solver
  
  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);
  
  return 0;
}
