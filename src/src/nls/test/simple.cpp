#include "simple.hpp"

SimpleNLDemo::SimpleNLDemo()
{
  // set initial condition on x
  x[0] = 2.0;
  x[1] = 0.5;
}

SimpleNLDemo::~SimpleNLDemo()
{
  // nothing 
}

void SimpleNLDemo::evaluate(petsc::Vector & r)
{
  double f[2];
  f[0] = x[0]*x[0] + x[1]*x[1] - 2.0;
  f[1] = exp(x[0]-1.0) + x[1]*x[1] - 2.0; 
  r.set(0, f[0]);
  r.set(1, f[1]);
  r.assemble();
}

void SimpleNLDemo::jacobian(petsc::Matrix & K)
{
  double a,b,c,d,x0,x1;
  x0 = x[0]; x1 = x[1];
  a = 2*x0;  b = 2*x1;
  c = exp(x0-1.0); d = 2*x1;
  K.set(0,0,a); K.set(0,1,b);
  K.set(1,0,c); K.set(1,1,d);
  K.assemble();
}

void SimpleNLDemo::update(petsc::Vector & u, double s)
{
  double dx[2];
  u.get_data(dx);
  x[0] = x[0] - dx[0];
  x[1] = x[1] - dx[1];
}

bool SimpleNLDemo::converged(petsc::Vector & du)
{
  double dx[2];
  double f0, f1, x0, x1, norm;
  
  du.get_data(dx);
   
  x0 = x[0];
  x1 = x[1];
  f0 = x0*x0 + x1*x1 - 2.0;
  f1 = exp(x0-1.0) + x1*x1 - 2.0;
  
  norm = sqrt(f0*f0 + f1*f1);
  
  std::cout << " F norm " << norm << std::endl;
    
  if(fabs(norm) < 1e-6) 
    return true;
  else
    return false;
}
