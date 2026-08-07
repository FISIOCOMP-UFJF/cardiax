#include "implicit_euler.hpp"

ImplicitEuler::ImplicitEuler(CellModel * ode)
  : ImplicitODESolver(ode)
{
  const int n = ode->get_num_state_vars();

  // Initialize vectors and matrices
  A.resize(n,n);
  b.resize(n);
  f.resize(n);
  yo.resize(n);
  yn.resize(n);
  yc.resize(n);
  
  // Initialize identity matrix only once
  I = arma::eye(n,n);

  // Message
  cout << " ODE solver: Implicit Euler" << endl;
}

void ImplicitEuler::advance(double * y, double & t, double & dt)
{
  int it, num_states = ode->get_num_state_vars();
  bool done;
  double dnorm, ynorm;
  double tn1 = t + dt;

  // Old y value (=y_n)
  for(int i=0; i < num_states; i++)
    yo(i) = y[i];
  
  yn = yo;      // Initial Newton iterate
  it = 0;       // Newton iteration counter
  done = false; // Iteration control

  // Main Newton iteration loop 
  while(!done)
  {
    yc = yn;

    // Evaluate ODE using computed solution
    ode->equation(tn1, yc.memptr(), f.memptr());
    
    // Build LHS of the linear system
    if (update_jacobian)
    {
      ode->compute_jacobian(yc.memptr(), tn1, jac);     
      A = I - dt*jac;
    }
     
    // Build RHS of the linear system
    b = yo - yc + dt*f;
    
    // Solve linear system
    dy = arma::solve(A,b);   
    yn = yc + dy;
    it = it + 1;
    
    dnorm = arma::norm(dy,2);
    ynorm = arma::norm(yn,2);

    if(dnorm < rtol*ynorm || it == maxits)
      done = true;
  }

  if (it==maxits)
    error("implicit_euler::newton did not converge");

  // Update
  for(int i=0; i < num_states; i++)
    y[i] = yn(i);
  
}
