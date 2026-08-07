#include "plane_stress_neo_hookean.hpp"

void PlaneStressNeoHookean::cauchy_stress(int ndim, double detf,
						 const arma::mat & btens,
						 arma::mat & sigma) const
{
  sigma.zeros();

  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      sigma(i,j) = mu * btens(i,j);
    }
    sigma(i,i) = sigma(i,i) - mu/(detf*detf);
  }
}

inline void PlaneStressNeoHookean::pk2_stress(double detf,
                                              const arma::mat & C,
                                              arma::mat & pk2) const
{
  throw std::runtime_error("Plane stress NH not implemented");
}

inline void PlaneStressNeoHookean::elas_tensor(ElasticityType elastype,
					       double detf, arma::mat & D) const
{
  double xmu = mu;
  double xme = xmu/(detf*detf);
  double xlambda = 2.0 * xme;
  double l2m = xlambda + 2.0 * xme;

  D.zeros();

  if (elastype == PLANE_STRESS)
  {
    D(0,0) = l2m;
    D(0,1) = xlambda;
    D(1,0) = xlambda;
    D(1,1) = l2m;
    D(2,2) = xme;
  }
  else
  {
    throw std::runtime_error("elasticity type not found");
  }
}

void PlaneStressNeoHookean::mat_elas_tensor (ElasticityType elastype,
					     double detf, const arma::mat & C,
					     arma::mat & D) const
{
  throw std::runtime_error("not implemented");
}



