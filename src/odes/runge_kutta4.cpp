#include "runge_kutta4.hpp"

RungeKutta4::RungeKutta4(CellModel *model)
  : ODESolver(model)
{
  cout << "\n ODE solver:"
       << "               "
       << "Runge Kutta 4" << endl;
}

void RungeKutta4::advance(double *sv, double &t, double & dt, double istim) 
{
  int n = ode->get_num_state_vars();
  double rY[n], k1[n], k2[n], k3[n], k4[n];
  //double htime;

  // evaluate k1
  for(int i=0; i<n; i++)
    rY[i] = sv[i];
  ode->equation(t,rY,k1, istim);
  
  // evaluate k2
  //htime = t + dt/2.;
  for (int i=0; i<n; i++)
    rY[i] = sv[i] + (dt/2.)*k1[i];
  ode->equation(t,rY,k2, istim);
  
  // evaluate k3
  for (int i = 0; i<n; i++)
    rY[i] = sv[i] + (dt/2.)*k2[i];
  ode->equation(t,rY,k3, istim);
  
  // evaluate k4
  //htime = t + dt;
  for (int i=0; i<n; i++)
    rY[i] = sv[i] + dt*k3[i];
  ode->equation(t,rY,k4, istim);
 
  for (int i=0; i<n; i++)
    sv[i] = sv[i] + (dt/6.)*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);
}

