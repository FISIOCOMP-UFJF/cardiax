#include "explicit_euler.hpp"

ExplicitEuler::ExplicitEuler(CellModel *model) 
  : ODESolver(model), dydt(model->get_num_state_vars())
{
  cout << "ODE solver: Explicit Euler" << endl;
}

void ExplicitEuler::advance(double * y, double & t, double & dt)
{
  const int n = ode->get_num_state_vars();

  // Evaluate ODE rhs: dydt = derivatives for all vars;
  // equation() also fills ode->rl_inf[i] / ode->rl_tau[i] for RL vars.
  ode->equation(t, y, dydt.memptr());

  const double* __restrict dy = dydt.memptr();

  // Advance
  for (int i = 0; i < n; i++)
  {
    if (ode->is_rl[i])
      y[i] = ode->rl_inf[i] + (y[i] - ode->rl_inf[i]) * std::exp(-dt / ode->rl_tau[i]);
    else
      y[i] += dt * dy[i];
  }
}
