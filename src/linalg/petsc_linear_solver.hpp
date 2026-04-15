#ifndef PETSC_LINEAR_SOLVER_HPP_
#define PETSC_LINEAR_SOLVER_HPP_

#include <iostream>
#include <map>
// #include <amgx_c.h>
#include "petscksp.h"
#include "petsc_matrix.hpp"
#include "petsc_vector.hpp"
#include "util/util.hpp"

namespace petsc
{

/** A simple wrapper around PETSc's KSP object for solving
    linear system of equations.
*/
class LinearSolver
{
public:

  //! Constructor
  LinearSolver()
  {}

  //! Destructor
  ~LinearSolver()
  {

    if(_ksp != NULL)
    {
      ierr = KSPDestroy(&_ksp);
      CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

  }

  //! Reason a Krylov method was said to have converged or diverged
  void converged_reason();


  //! Create KSP and call set from options
  void init();

  //! Check whether the KSP context is null
  bool is_null() const
  { return (_ksp == NULL); }

  //! Change solver type (cg,gmres,bicg,...)
  void set_solver_type();

  //! Change solver type (cg,gmres,bicg,...)
  void set_solver_type(std::string & type);

  //! Change solver type (cg,gmres,bicg,...)
  void set_solver_type(const char * type);

  //! Change preconditioner type (ilu,icc,jacobi,...)
  void set_preconditioner(std::string & type);

  //! Change preconditioner type (ilu,icc,jacobi,...)
  void set_preconditioner(const char * type);

  //! Change the ordering of the matrix
  void set_ordering(const char * otype);

  //! Uses UMFPACK LU as a direct solver
  void use_umfpack();

  //! Uses MUMPS LU as a direct solver
  void use_mumps();

  //! Uses HYPRE 
  void use_hypre();

  //! Uses PC FIELDSPLIT for displacement components
  void use_fieldsplit(int N);
  void use_fieldsplit_2D(int N);
  
  //! Print KSP object information
  void view();

  void MatlabView(petsc::Matrix & A);



  //! Solve the linear system
  std::pair<PetscInt, PetscReal> solve (petsc::Matrix & A,
                                        petsc::Vector & x,
                                        petsc::Vector & b,
                                        const double tol=1.0e-16);
  std::pair<PetscInt, PetscReal> solveFieldSplit (petsc::Matrix & A,
                                          petsc::Vector & x,
                                          petsc::Vector & b,
                                          const double tol=1.0e-16);
  std::pair<PetscInt, PetscReal> solve_3_FieldSplit (petsc::Matrix & A,
                                                    petsc::Vector & x,
                                                    petsc::Vector & b,
                                                    const double tol=1.0e-16);

private:

  //! The KSP context object
  KSP _ksp;

  //! The preconditioner object to be used for the solution
  PC _pc;

  //! PETSc int error code
  PetscErrorCode ierr;

  const char *getConfig(const std::string &configName);

};

}
#endif /* PETSC_LINEAR_SOLVER_HPP_ */
