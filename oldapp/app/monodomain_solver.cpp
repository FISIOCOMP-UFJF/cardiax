#include <iostream>
#include <string>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "cardiac/monodomain.hpp"
#include "cardiac/monodomain_purkinje.hpp"
#include "cardiac/coupled_purkinje.hpp"
#include<filesystem>

static char help[] = "Monodomain Cardiac Solver.\n\n";

using namespace std;
namespace fs = std::filesystem;

void usage()
{
  cout << endl;
  cout << " Usage: ";
  cout << "monodomain [OPTIONS]" << endl << endl;
  cout << "    -f     <meshbase>    prefix of the mesh file" << endl;
  cout << "    -dt    <dt>          time step (ms)" << endl;
  cout << "    -t     <tend>        total time of simulation (ms)" << endl;
  cout << "    -pr    <pr>          print rate to output file" << endl;
  cout << "    -c     <cellmodel>   string that identifies the ionic model" << endl;
  cout << "    -m     <odesolver>   ODE solver (ExplicitEuler, Implicit, ...)" << endl;
	cout << "    -ep    <model>       monodomain or purkinje" << endl;
  cout << "    -fp    <pkmesh>      Purkinje mesh" << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char *argv[])
{
  double dt, T, tp;
  string mshname, pkmshname, cellmodel, odesolver, typefile;
	string model;

  if (argc < 7) usage();

  // Parse command line options
  CommandLineArgs::init(argc, argv);
  dt = CommandLineArgs::read("-dt",0.1);
  T  = CommandLineArgs::read("-t",10.0);
  tp = CommandLineArgs::read("-pr",1.0);
  mshname = CommandLineArgs::read("-f", "emptymesh");
  pkmshname = CommandLineArgs::read("-fp", "emptymesh");
  cellmodel = CommandLineArgs::read("-c","TT2"); 
  odesolver = CommandLineArgs::read("-m","ExplicitEuler");
	model = CommandLineArgs::read("-ep","monodomain");

	typefile  = mshname + ".typ";

  if ( fs::exists("output") ) 
  {
    fs::remove_all("output");
    assert(!fs::exists("output"));
  }
  
  fs::create_directory("output");
  fs::create_directories("output/vm_vtu/");
  fs::create_directories("output/vm_text/");

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

	if(model == "monodomain")
	{
		// Start PDE solver for Monodomain model
		cout << "Monodomain solver\n";
		{
			Monodomain monodomain;
      msg("Reading parameters file");
//      monodomain.config(mshname);
			monodomain.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
			monodomain.init();
			
			if(!file_exists(typefile))
				cout << "Cells: all cells are of the same type\n";
			else    
				monodomain.setup_types(typefile);
			
			monodomain.solve();   
		}		
	}
	else if(model == "purkinje")
	{
		cout << "Monodomain-Purkinje Solver\n";
		{
			MonodomainPurkinje mp;
			mp.setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
			mp.init();					
			mp.solve();   
		}		
	}
  else if(model == "coupled_purkinje")
  {
    cout << "Coupled tissue/Purkinje Monodomain Solver\n";
    {
      CoupledPurkinje cp;
      if(pkmshname == "emptymesh")
      {
        std::cerr << "Error: unknown -fp PurkinjeMesh" << endl;
        exit(1);
      }
      cp.setup(mshname, pkmshname, cellmodel, odesolver, dt, T, tp, tp);
      cp.solve();
    }
  }
	else
	{
		cout << "Cardiac PDE Model " << model << " does not exist." << endl;
	}
	
  // End PDE solver
  
  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);
  
  return 0;
}
