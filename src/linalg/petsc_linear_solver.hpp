#ifndef PETSC_LINEAR_SOLVER_HPP_
#define PETSC_LINEAR_SOLVER_HPP_

#include <iostream>
#include <map>
#ifdef AMGX_SOLVER
  #include <amgx_c.h>
#endif

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

  #ifdef AMGX_SOLVER
    // Only tear down what this instance actually built. Destroying a null or
    // never-created handle is not a no-op in AMGX: it reads the mode out of
    // the handle and aborts with "Incorrect C API mode".
    if (_amgx_ready)
    {
      if (_amgx_solver) AMGX_SAFE_CALL(AMGX_solver_destroy(_amgx_solver));
      if (_amgx_x)      AMGX_SAFE_CALL(AMGX_vector_destroy(_amgx_x));
      if (_amgx_b)      AMGX_SAFE_CALL(AMGX_vector_destroy(_amgx_b));
      if (_amgx_A)      AMGX_SAFE_CALL(AMGX_matrix_destroy(_amgx_A));
      if (_amgx_rsrc)   AMGX_SAFE_CALL(AMGX_resources_destroy(_amgx_rsrc));
      if (_amgx_config) AMGX_SAFE_CALL(AMGX_config_destroy(_amgx_config));

      _amgx_solver = nullptr; _amgx_x = nullptr; _amgx_b = nullptr;
      _amgx_A = nullptr; _amgx_rsrc = nullptr; _amgx_config = nullptr;
      _amgx_ready = false;

      // AMGX_finalize() is global. Calling it from every destructor meant the
      // first LinearSolver to die shut the library down for all the others,
      // which then operated on a finalised AMGX. Finalise only when the last
      // holder goes away.
      if (--_amgx_live == 0)
        AMGX_SAFE_CALL(AMGX_finalize());
    }
  #endif
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

  //! The KSP context object.
  //! Initialised here, not only in init(): the destructor tests it against
  //! NULL, and several LinearSolver instances in the code base are
  //! constructed but never init()'d (Monodomain::solver and
  //! MonodomainDeformation::solver only call init() under `if(test)`).
  //! Without this initialiser those instances reach the destructor holding
  //! whatever was on the stack.
  KSP _ksp = NULL;

  #ifdef AMGX_SOLVER
      // AMGX handles. Same reasoning as _ksp, and the consequence is worse
      // here: AMGX decodes the mode (precision, host/device) from the handle
      // itself, so destroying an uninitialised one makes the library look up
      // a garbage mode and report "Mode not found" / "Incorrect C API mode".
      // That error surfaces at teardown, far from the code that caused it.
      AMGX_config_handle    _amgx_config  = nullptr;
      AMGX_resources_handle _amgx_rsrc    = nullptr;
      AMGX_matrix_handle    _amgx_A       = nullptr;
      AMGX_vector_handle    _amgx_b       = nullptr;
      AMGX_vector_handle    _amgx_x       = nullptr;
      AMGX_solver_handle    _amgx_solver  = nullptr;

      //! True once init() has built the AMGX objects of THIS instance.
      bool _amgx_ready = false;

      //! Number of instances currently holding AMGX objects. AMGX_initialize
      //! and AMGX_finalize are global, so they must be paired once for the
      //! whole process, not once per LinearSolver.
      static int _amgx_live;
  #endif

  //! The preconditioner object to be used for the solution
  PC _pc;

  //! PETSc int error code
  PetscErrorCode ierr;

  const char *getConfig(const std::string &configName);

};

}
#endif /* PETSC_LINEAR_SOLVER_HPP_ */
