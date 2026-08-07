#include "ode_solver.hpp"
#include "cellmodel.hpp"
#include "runge_kutta4.hpp"
#include "explicit_euler.hpp"
#include "implicit_euler.hpp"

ODESolver * ODESolver::create(std::string method, CellModel *cell)
{
  ODESolver *ptr = NULL;

  if (method == "ExplicitEuler")
    ptr = new ExplicitEuler(cell);
  else if (method == "RungeKutta4")
    ptr = new RungeKutta4(cell);
  else if (method == "ImplicitEuler")
    ptr = new ImplicitEuler(cell);
  else
    throw std::invalid_argument("ode_solver::Unknown ODESolver.");

  return ptr;
}
