#include "lafortune.hpp"

double LaFortune::strain_energy(MaterialData * md, const arma::mat & E) const
{
  // Undeformed fiber and directions
  arma::vec3 & f0 = md->f;
  arma::vec3 & s0 = md->s;

  // Debug
#ifdef DEBUG
  f0.zeros(); s0.zeros();
  f0(0)=1.;
  s0(1)=1.;
#endif

  // Right Cauchy-Green deformation tensor and its trace
  arma::mat C = 2*E + I;

  // Jacobian
  double J = sqrt(arma::det(C));

  // Isochoric right Cauchy-Green deformation tensor and its trace
  arma::mat Cbar = pow(J,(-2.0/3.0)) * C; 
  double Icbar = arma::trace(Cbar);

  // Invariants
  const arma::vec3 cs = C*s0;
  double I4f  = dot(f0, C*f0);
  double I4s  = dot(s0, cs);
  double I8fs = dot(f0, cs);

  const double K = 100;
  double term1 = (a/(2.0*b))*(exp(b*(Icbar-3.0))-1.0);
  double term2 = -(a/2)*(Icbar-3);
  double term3;
  if(I4f > 1.0)
    term3 = (af/(2*bf))*(exp(bf*(I4f-1.0)*(I4f-1.0))-1.0);
  else
    term3 = ((6.5*af)/(2*bf))*(exp(bf*(I4f-1.0)*(I4f-1.0))-1.0);
  double term4 = (K/2.0)*((J-1)*(J-1)); 
  return term1 + term2 + term3 + term4;
}

void LaFortune::cauchy_stress(int ndim, double detf, const arma::mat & B,
                              arma::mat & sigma) const
{
  throw runtime_error("LaFortune: not implemented");
}

void LaFortune::piola2_stress(double detf, const arma::mat & C, arma::mat & S) const
{
  throw runtime_error("LaFortune: not implemented");
}

void LaFortune::sp_elastensor (ElasticityType elastype, double detf,
                               arma::mat & D) const
{
  throw runtime_error("LaFortune: not implemented");
}

void LaFortune::mt_elastensor (ElasticityType elastype,
                               double detf, const arma::mat & C,
                               arma::mat & D) const
{
  throw runtime_error("LaFortune: not implemented");
}

