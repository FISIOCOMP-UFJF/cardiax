#include <cstdlib>
#include "../pardiso_solver.hpp"

using namespace std;
static char help[] = "PARDISO test\n\n";

/*!
 * Simple program with a simple linear system taken from
 * PARDISO manual to test the implementation and wrapper.
 */
int main (int argc, char** argv)
{
   // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
 
  // System
  petsc::Matrix * mat = new petsc::Matrix();
  mat->create(8,8,10);
  mat->set_symmetric();
  
  mat->set(0,0,7);
  mat->set(0,2,1); //mat->set(2,0,1);
  mat->set(0,5,2); //mat->set(5,0,2);
  mat->set(0,6,7); //mat->set(6,0,7);
  
  mat->set(1,1,-4);
  mat->set(1,2,8); //mat->set(2,1,8);
  mat->set(1,4,2); //mat->set(4,1,2);
  
  mat->set(2,2,1); 
  mat->set(2,7,5); //mat->set(7,2,5);
  
  mat->set(3,3,7);
  mat->set(3,6,9); //mat->set(6,3,9);
  
  mat->set(4,4,5);
  mat->set(4,5,1); //mat->set(5,4,-1);
  mat->set(4,6,5); //mat->set(6,4,5);
 
  mat->set(5,5,0);
  mat->set(5,7,5); //mat->set(7,5,5);
  
  mat->set(6,6,11);
  
  mat->set(7,7,5);    
  
  mat->assemble();    
  mat->view();
  
  double b[8], x[8];
  for (int i = 0; i < 8; i++) { x[i] = 0; }
  for (int i = 0; i < 8; i++) { b[i] = i; } 
 
  // Call solver to solve the system
  PardisoSolver solver;    
  solver.solve(*mat, x, b);
  
  delete mat;
  
  // End
  
  ierr = PetscFinalize(); 
  CHKERRQ(ierr);
 
  return 0;
}

