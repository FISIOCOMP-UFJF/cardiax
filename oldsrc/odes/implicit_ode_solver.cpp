#include "implicit_ode_solver.hpp"

ImplicitODESolver::ImplicitODESolver(CellModel * ode)
  : ODESolver(ode), maxits(10), jsize(ode->get_num_state_vars()),
    update_jacobian(true), atol(1e-5), rtol(1e-5)
{
  const int n = ode->get_num_state_vars();

  // Initialize matrices and vectors
  jac.resize(n,n);
  dy.resize(n);
}
