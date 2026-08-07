#include "exponential.hpp"

double Exponential::strain_energy(MaterialData * md, const arma::mat & E) const
{
  // Right Cauchy-Green deformation tensor and its trace
  arma::mat C = 2*E + arma::eye<arma::mat>(3,3); 
 
  // Jacobian
  double J = sqrt(arma::det(C));

  // Isochoric right Cauchy-Green deformation tensor and its trace
  arma::mat Cbar = pow(J,(-2.0/3.0)) * C;
  double Icbar = trace(Cbar);
  
  // Compute strain energy function terms
  double term1 = (a/(2.0*b))*(exp(b*(Icbar-3.0))-1.0);
  double term2 = (3333.0/2.0)*((J-1)*(J-1));  // 333 or 3333
  return term1 + term2;
}

double Exponential::active_strain_energy(MaterialData * md,
                                      const arma::mat &) const
{
  return 0;
}

void Exponential::cauchy_stress (MaterialData * md, arma::mat & sigma) const
{
  throw std::runtime_error("Exponential -> not implemented");
}

void Exponential::piola2_stress(MaterialData * md, arma::mat & S) const
{
  throw std::runtime_error("Exponential -> not implemented");
}

void Exponential::sp_elastensor(MaterialData * md, arma::mat & D) const
{
  throw std::runtime_error("Exponential -> not implemented");
}

void Exponential::mt_elastensor(MaterialData * md, arma::mat & D) const
{
  throw std::runtime_error("Exponential -> not implemented");
}

