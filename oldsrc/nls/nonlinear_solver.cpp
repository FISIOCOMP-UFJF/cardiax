#include "nonlinear_solver.hpp"

NonlinearSolver::NonlinearSolver(NonlinearProblem * p) : nlp(p)
{
  ls = new petsc::LinearSolver();
  ls->init();

  //ls->use_fieldsplit(p->size());

  //ls->use_mumps();
  //ls->set_ordering("rcm");

//  ls->use_umfpack();
//  ls->set_ordering("amd");
}

NonlinearSolver::~NonlinearSolver()
{
  delete ls;
}
