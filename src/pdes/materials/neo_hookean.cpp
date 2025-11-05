#include "neo_hookean.hpp"

double NeoHookean::strain_energy(MaterialData * md, const arma::mat & E) const
{
  // Right Cauchy-Green deformation tensor
  arma::mat C = 2*E + arma::eye<arma::mat>(3,3); 

  // Jacobian
  double J = sqrt(arma::det(C));

  // Isochoric right Cauchy-Green deformation tensor
  arma::mat Cbar = pow(J,(-2.0/3.0)) * C;
  
  // Trace of isochoric right Cauchy-Green tensor
  double Icbar = trace(Cbar);

  //return (mu/2.0)*(Icbar - 3.0) + (kappa/4.0)*((J-1)*(J-1)+(log(J)*log(J)));
  return (mu/2.0)*(Icbar - 3.0 - 2*log(J));
    //+ (kappa/2.0)*((log(J)*log(J)));
}

double NeoHookean::active_strain_energy(MaterialData * md,
                                         const arma::mat &) const
{
  return 0;
}

void NeoHookean::deviatoric_stress(MaterialData * md, arma::mat & sigma) const
{
  double J = md->get_J();
  arma::mat33 btens = md->left_cauchy_green();
  double xm = mu/(pow(J,(5./3.)));
  double trace = arma::trace(btens);

  sigma.zeros();

  for(int i=0; i<ndim; i++)
  {
    for(int j=0; j<ndim; j++)
      sigma(i,j) = xm * btens(i,j);
    sigma(i,i) = sigma(i,i) - xm * trace/3.;
  }
}

void NeoHookean::deviatoric_elastensor(MaterialData * md, Tensor4 & A) const
{
  double J, xm, tr;
  arma::mat b = md->left_cauchy_green();
  static arma::mat33 delta = arma::eye(3,3);

  J  = md->get_J();
  xm = mu/(3.*pow(J,(5./3.)));
  tr = arma::trace(b); 

  A.zero();  

  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
        for(int l=0; l<ndim; l++)
          A(i,j,k,l) = (2.*xm*tr/3.0) * delta(i,j) * delta(k,l)
                     + xm * tr * delta(i,k) * delta(j,l)
	                 + xm * tr * delta(i,l) * delta(j,k)
	                 - 2. * xm * b(i,j) * delta(k,l)
	                 - 2. * xm * delta(i,j) * b(k, l);
}

/*
void NeoHookean::cisotp(double detf, Tensor4 & A)
{  
  double xlamb = lambda/detf;
  double xmu = mu/detf;

  arma::mat delta = arma::eye<arma::mat>(ndim,ndim);

  A.zero();

  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
	for(int l=0; l<ndim; l++)
	  A(i,j,k,l) = xlamb * delta(i,j) * delta(k,l) 
                     + xmu   * delta(i,k) * delta(j,l)
                     + xmu   * delta(i,l) * delta(j,k);

}
*/


