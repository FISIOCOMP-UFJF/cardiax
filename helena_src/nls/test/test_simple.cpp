/* 
 * File:   test_simple.cpp
 * Author: rocha
 *
 * Created on January 31, 2014, 3:51 PM
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "../newton.hpp"
#include "simple.hpp"

using namespace std;

static char help[] = "NLSolver test\n\n";

int main(int argc, char** argv)
{
  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  { 
    cout << "Creating simple problem" << endl;
    SimpleNLDemo * s = new SimpleNLDemo();
    
    cout << "Creating Newton solver for simple problem" << endl;
    Newton * nls = new Newton(s);
    
    cout << "Initializing Newton" << endl;
    nls->init();
    
    cout << "Solving problem" << endl;
    nls->solve();  
    
    cout << "Freeing memory" << endl;
    delete nls;
  }
  
  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);
  
  return 0;
}

