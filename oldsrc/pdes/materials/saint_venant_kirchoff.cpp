#include "saint_venant_kirchoff.hpp"

double SaintVenantKirchoff::strain_energy(MaterialData * md, const arma::mat & E) const
{
  return 0.5*lambda*pow(arma::trace(E),2) + mu*contract(E,E);
}

double SaintVenantKirchoff::active_strain_energy(MaterialData * md,
                                             const arma::mat &) const
{
  return 0;
}

void SaintVenantKirchoff::cauchy_stress (MaterialData * md, arma::mat & sigma) const
{
  throw std::runtime_error("StVK: not implemented");
}

void SaintVenantKirchoff::piola2_stress(MaterialData * md, arma::mat & S) const
{
  arma::mat C = md->right_cauchy_green();
  arma::mat E = 0.5*(C - arma::eye<arma::mat>(3,3));
  S = lambda*arma::trace(E)*arma::eye<arma::mat>(3,3) + 2*mu*E;
}

void SaintVenantKirchoff::sp_elastensor(MaterialData * md, arma::mat & D) const
{
  throw runtime_error("StVK: not implemented");
}

void SaintVenantKirchoff::mt_elastensor(MaterialData * md, arma::mat & D) const
{
  D.zeros();

  if (ndim == 2)
  {
    error("2D not implemented yet!");
  }
  else if (ndim == 3)
  {
    D(0,0) = lambda + 2.0*mu;
    D(0,1) = lambda;
    D(0,2) = lambda;

    D(1,0) = lambda;
    D(1,1) = lambda + 2.0*mu;
    D(1,2) = lambda;

    D(2,0) = lambda;
    D(2,1) = lambda;
    D(2,2) = lambda + 2.0*mu;
    
    D(3,3) = mu;
    D(4,4) = mu;
    D(5,5) = mu;
  }
  else
  {
    throw runtime_error("elasticity type not found");
  }
}
