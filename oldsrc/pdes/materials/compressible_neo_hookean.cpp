#include "compressible_neo_hookean.hpp"

double CompressibleNeoHookean::strain_energy(MaterialData * md, const arma::mat & E) const
{
  arma::mat C = 2*E + arma::eye<arma::mat>(3,3);
  double Ic = trace(C);
  double J = sqrt(arma::det(C));
  return (mu/2.)*(Ic-3.) - mu*log(J) + (lambda/2.)*(log(J)*log(J));
}

double CompressibleNeoHookean::active_strain_energy(MaterialData * md,
                                                 const arma::mat &) const
{
  return 0;
}

inline void CompressibleNeoHookean::cauchy_stress(MaterialData * md,
                                                  arma::mat & sigma) const
{
  double J = md->get_J();
  double xme = (mu - lambda * log(J))/J;
  arma::mat b = md->left_cauchy_green();

  sigma.zeros();

  for(int i=0; i<ndim; i++)
  {
    for(int j=0; j<ndim; j++)
      sigma(i,j) = (mu/J) * b(i,j);
    sigma(i,i) = sigma(i,i) - xme;
  }

  if (ndim == 2)
    sigma(2,2) = (mu/J)*(b(2,2) - 1.0) + (lambda/J)*log(J);
}

inline void CompressibleNeoHookean::piola2_stress(MaterialData * md, arma::mat & S) const
{
  double J = md->get_J();
  arma::mat C = md->left_cauchy_green();
  double c11 = C(0,0);
  double c12 = C(0,1);
  double c13 = C(0,2);
  double c21 = C(1,0);
  double c22 = C(1,1);
  double c23 = C(1,2);
  double c31 = C(2,0);
  double c32 = C(2,1);
  double c33 = C(2,2);

  S.zeros();

  S(0,0) =  1.0*mu + 2.0*lambda*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(c22*c33/2 - c23*c32/2)/pow(J, 2);
  S(0,1) =  2.0*lambda*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(-c21*c33/2 + c23*c31/2)/pow(J, 2);
  S(0,2) =  2.0*lambda*(c21*c32/2 - c22*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(c21*c32/2 - c22*c31/2)/pow(J, 2);
  S(1,0) =  2.0*lambda*(-c12*c33/2 + c13*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(-c12*c33/2 + c13*c32/2)/pow(J, 2);
  S(1,1) =  1.0*mu + 2.0*lambda*(c11*c33/2 - c13*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(c11*c33/2 - c13*c31/2)/pow(J, 2);
  S(1,2) =  2.0*lambda*(-c11*c32/2 + c12*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(-c11*c32/2 + c12*c31/2)/pow(J, 2);
  S(2,0) =  2.0*lambda*(c12*c23/2 - c13*c22/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(c12*c23/2 - c13*c22/2)/pow(J, 2);
  S(2,1) =  2.0*lambda*(-c11*c23/2 + c13*c21/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(-c11*c23/2 + c13*c21/2)/pow(J, 2);
  S(2,2) =  1.0*mu + 2.0*lambda*(c11*c22/2 - c12*c21/2)*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*mu*(c11*c22/2 - c12*c21/2)/pow(J, 2);
}


inline void CompressibleNeoHookean::sp_elastensor(MaterialData * md, 
						  arma::mat & D) const
{
  double J = md->get_J();
  double xlambda = lambda/J;
  double xmu = mu/J;
  double xme = (xmu - xlambda * log(J))/J;
  double l2m = xlambda + 2.0*xme;

  D.zeros();

  if (ndim == 2)
  {
    D(0,0) = l2m;
    D(0,1) = xlambda;
    D(1,0) = xlambda;
    D(1,1) = l2m;
    D(2,2) = xme;
  }
  else if (ndim == 3)
  {
    D(0,0) = l2m;
    D(0,1) = xlambda;
    D(0,2) = xlambda;
    D(1,0) = xlambda;
    D(1,1) = l2m;
    D(1,2) = xlambda;
    D(2,0) = xlambda;
    D(2,1) = xlambda;
    D(2,2) = l2m;
    D(3,3) = xme;
    D(4,4) = xme;
    D(5,5) = xme;
  }
}

void CompressibleNeoHookean::mt_elastensor(MaterialData * md, 
					   arma::mat & D) const
{
  arma::mat33 C = md->left_cauchy_green();
  double J = md->get_J();
  double c11 = C(0,0);
  double c12 = C(0,1);
  double c13 = C(0,2);
  double c21 = C(1,0);
  double c22 = C(1,1);
  double c23 = C(1,2);
  double c31 = C(2,0);
  double c32 = C(2,1);
  double c33 = C(2,2);

  D.zeros();

  if(ndim == 2)
  {
    D(0,0)  = 4.0*lambda*(-c22*c33 + c23*c32)*(c22*c33/2 - c23*c32/2)*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) + 4.0*lambda*pow(c22*c33/2 - c23*c32/2, 2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) - 4.0*mu*(-c22*c33 + c23*c32)*(c22*c33/2 - c23*c32/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2);
    D(0,1)  = 2.0*c33*lambda*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31) - 2.0*c33*mu/(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31) + 4.0*lambda*(-c11*c33 + c13*c31)*(c22*c33/2 - c23*c32/2)*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(c22*c33/2 - c23*c32/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) - 4.0*mu*(-c11*c33 + c13*c31)*(c22*c33/2 - c23*c32/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2);
    D(0,2)  = 4.0*lambda*(-c21*c33/2 + c23*c31/2)*(c22*c33/2 - c23*c32/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) + 4.0*lambda*(c21*c33 - c23*c31)*(c22*c33/2 - c23*c32/2)*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) - 4.0*mu*(c21*c33 - c23*c31)*(c22*c33/2 - c23*c32/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2);
    D(0,2) /= 2.0;

    D(1,0)  = D(0,1); // symmetric
    D(1,1)  = 4.0*lambda*(-c11*c33 + c13*c31)*(c11*c33/2 - c13*c31/2)*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) + 4.0*lambda*pow(c11*c33/2 - c13*c31/2, 2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) - 4.0*mu*(-c11*c33 + c13*c31)*(c11*c33/2 - c13*c31/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2);
    D(1,2)  = 4.0*lambda*(c11*c33/2 - c13*c31/2)*(-c21*c33/2 + c23*c31/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(c21*c33 - c23*c31)*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) - 4.0*mu*(c11*c33/2 - c13*c31/2)*(c21*c33 - c23*c31)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2);
    D(1,2) /= 2.0;

    D(2,0)  = D(0,2); // symmetric
    D(2,1)  = D(1,2); // symmetric
    D(2,2)  = -2.0*c33*lambda*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31) + 2.0*c33*mu/(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31) + 4.0*lambda*(-c12*c33/2 + c13*c32/2)*(-c21*c33/2 + c23*c31/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) + 4.0*lambda*(c12*c33 - c13*c32)*(-c21*c33/2 + c23*c31/2)*log(sqrt(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31))/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2) - 4.0*mu*(c12*c33 - c13*c32)*(-c21*c33/2 + c23*c31/2)/pow(c11*c22*c33 - c11*c23*c32 - c12*c21*c33 + c12*c23*c31 + c13*c21*c32 - c13*c22*c31, 2);
    D(2,2) /= 2.0;
  }
  else if (ndim == 3)
  {
    //cout << "Novo D -> Bonet" << endl;

    D(0,0) = 4.0*lambda*(-c22*c33 + c23*c32)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*pow(c22*c33/2 - c23*c32/2, 2)/pow(J, 4) - 4.0*mu*(-c22*c33 + c23*c32)*(c22*c33/2 - c23*c32/2)/pow(J, 4);
    D(0,1) = 2.0*c33*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*c33*mu/pow(J, 2) + 4.0*lambda*(-c11*c33 + c13*c31)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(-c11*c33 + c13*c31)*(c22*c33/2 - c23*c32/2)/pow(J, 4);
    D(0,2) = 2.0*c22*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*c22*mu/pow(J, 2) + 4.0*lambda*(-c11*c22 + c12*c21)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(-c11*c22 + c12*c21)*(c22*c33/2 - c23*c32/2)/pow(J, 4);
    D(0,3) = 2.0*lambda*(-c12*c33/2 + c13*c32/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) + 2.0*lambda*(c12*c33 - c13*c32)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c21*c33/2 + c23*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) + 2.0*lambda*(c21*c33 - c23*c31)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c12*c33 - c13*c32)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 2.0*mu*(c21*c33 - c23*c31)*(c22*c33/2 - c23*c32/2)/pow(J, 4);
    D(0,4) = -1.0*c23*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c23*mu/pow(J, 2) - 1.0*c32*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c32*mu/pow(J, 2) + 2.0*lambda*(-c11*c23/2 + c13*c21/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) + 2.0*lambda*(c11*c23 - c13*c21)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) + 2.0*lambda*(c11*c32 - c12*c31)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c23 - c13*c21)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 2.0*mu*(c11*c32 - c12*c31)*(c22*c33/2 - c23*c32/2)/pow(J, 4);
    D(0,5) = 2.0*lambda*(-c12*c23 + c13*c22)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) + 2.0*lambda*(-c21*c32 + c22*c31)*(c22*c33/2 - c23*c32/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c21*c32/2 - c22*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 2.0*mu*(-c12*c23 + c13*c22)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 2.0*mu*(-c21*c32 + c22*c31)*(c22*c33/2 - c23*c32/2)/pow(J, 4);

    D(1,0) = 2.0*c33*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*c33*mu/pow(J, 2) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(-c22*c33 + c23*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(c11*c33/2 - c13*c31/2)*(-c22*c33 + c23*c32)/pow(J, 4);
    D(1,1) = 4.0*lambda*(-c11*c33 + c13*c31)*(c11*c33/2 - c13*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*pow(c11*c33/2 - c13*c31/2, 2)/pow(J, 4) - 4.0*mu*(-c11*c33 + c13*c31)*(c11*c33/2 - c13*c31/2)/pow(J, 4);
    D(1,2) = 2.0*c11*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*c11*mu/pow(J, 2) + 4.0*lambda*(-c11*c22 + c12*c21)*(c11*c33/2 - c13*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(c11*c33/2 - c13*c31/2)/pow(J, 4) - 4.0*mu*(-c11*c22 + c12*c21)*(c11*c33/2 - c13*c31/2)/pow(J, 4);
    D(1,3) = 2.0*lambda*(c11*c33/2 - c13*c31/2)*(-c12*c33/2 + c13*c32/2)/pow(J, 4) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(c12*c33 - c13*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(c21*c33 - c23*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c33/2 - c13*c31/2)*(c12*c33 - c13*c32)/pow(J, 4) - 2.0*mu*(c11*c33/2 - c13*c31/2)*(c21*c33 - c23*c31)/pow(J, 4);
    D(1,4) = 2.0*lambda*(-c11*c23/2 + c13*c21/2)*(c11*c33/2 - c13*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c23 - c13*c21)*(c11*c33/2 - c13*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c11*c33/2 - c13*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c32 - c12*c31)*(c11*c33/2 - c13*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c23 - c13*c21)*(c11*c33/2 - c13*c31/2)/pow(J, 4) - 2.0*mu*(c11*c32 - c12*c31)*(c11*c33/2 - c13*c31/2)/pow(J, 4);
    D(1,5) = -1.0*c13*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c13*mu/pow(J, 2) - 1.0*c31*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c31*mu/pow(J, 2) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(-c12*c23 + c13*c22)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(-c21*c32 + c22*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c33/2 - c13*c31/2)*(c21*c32/2 - c22*c31/2)/pow(J, 4) - 2.0*mu*(c11*c33/2 - c13*c31/2)*(-c12*c23 + c13*c22)/pow(J, 4) - 2.0*mu*(c11*c33/2 - c13*c31/2)*(-c21*c32 + c22*c31)/pow(J, 4);

    D(2,0) = 2.0*c22*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*c22*mu/pow(J, 2) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(-c22*c33 + c23*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(c11*c22/2 - c12*c21/2)*(-c22*c33 + c23*c32)/pow(J, 4);
    D(2,1) = 2.0*c11*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 2.0*c11*mu/pow(J, 2) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(-c11*c33 + c13*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(c11*c33/2 - c13*c31/2)/pow(J, 4) - 4.0*mu*(c11*c22/2 - c12*c21/2)*(-c11*c33 + c13*c31)/pow(J, 4);
    D(2,2) = 4.0*lambda*(-c11*c22 + c12*c21)*(c11*c22/2 - c12*c21/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*pow(c11*c22/2 - c12*c21/2, 2)/pow(J, 4) - 4.0*mu*(-c11*c22 + c12*c21)*(c11*c22/2 - c12*c21/2)/pow(J, 4);
    D(2,3) = -1.0*c12*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c12*mu/pow(J, 2) - 1.0*c21*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c21*mu/pow(J, 2) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(-c12*c33/2 + c13*c32/2)/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(c12*c33 - c13*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(c21*c33 - c23*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c22/2 - c12*c21/2)*(c12*c33 - c13*c32)/pow(J, 4) - 2.0*mu*(c11*c22/2 - c12*c21/2)*(c21*c33 - c23*c31)/pow(J, 4);
    D(2,4) = 2.0*lambda*(c11*c22/2 - c12*c21/2)*(-c11*c23/2 + c13*c21/2)/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(c11*c23 - c13*c21)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(-c11*c32/2 + c12*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(c11*c32 - c12*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c22/2 - c12*c21/2)*(c11*c23 - c13*c21)/pow(J, 4) - 2.0*mu*(c11*c22/2 - c12*c21/2)*(c11*c32 - c12*c31)/pow(J, 4);
    D(2,5) = 2.0*lambda*(c11*c22/2 - c12*c21/2)*(-c12*c23 + c13*c22)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(-c21*c32 + c22*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c11*c22/2 - c12*c21/2)*(c21*c32/2 - c22*c31/2)/pow(J, 4) - 2.0*mu*(c11*c22/2 - c12*c21/2)*(-c12*c23 + c13*c22)/pow(J, 4) - 2.0*mu*(c11*c22/2 - c12*c21/2)*(-c21*c32 + c22*c31)/pow(J, 4);

    D(3,0) = 4.0*lambda*(-c21*c33/2 + c23*c31/2)*(-c22*c33 + c23*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(-c21*c33/2 + c23*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(-c21*c33/2 + c23*c31/2)*(-c22*c33 + c23*c32)/pow(J, 4);
    D(3,1) = 4.0*lambda*(-c11*c33 + c13*c31)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) - 4.0*mu*(-c11*c33 + c13*c31)*(-c21*c33/2 + c23*c31/2)/pow(J, 4);
    D(3,2) = -2.0*c21*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 2.0*c21*mu/pow(J, 2) + 4.0*lambda*(-c11*c22 + c12*c21)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) - 4.0*mu*(-c11*c22 + c12*c21)*(-c21*c33/2 + c23*c31/2)/pow(J, 4);
    D(3,3) = -1.0*c33*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c33*mu/pow(J, 2) + 2.0*lambda*(-c12*c33/2 + c13*c32/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(c12*c33 - c13*c32)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*pow(-c21*c33/2 + c23*c31/2, 2)/pow(J, 4) + 2.0*lambda*(-c21*c33/2 + c23*c31/2)*(c21*c33 - c23*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c12*c33 - c13*c32)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) - 2.0*mu*(-c21*c33/2 + c23*c31/2)*(c21*c33 - c23*c31)/pow(J, 4);
    D(3,4) = 1.0*c31*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 1.0*c31*mu/pow(J, 2) + 2.0*lambda*(-c11*c23/2 + c13*c21/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c23 - c13*c21)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c32 - c12*c31)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c23 - c13*c21)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) - 2.0*mu*(c11*c32 - c12*c31)*(-c21*c33/2 + c23*c31/2)/pow(J, 4);
    D(3,5) = 1.0*c23*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 1.0*c23*mu/pow(J, 2) + 2.0*lambda*(-c12*c23 + c13*c22)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(-c21*c32 + c22*c31)*(-c21*c33/2 + c23*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c21*c32/2 - c22*c31/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) - 2.0*mu*(-c12*c23 + c13*c22)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) - 2.0*mu*(-c21*c32 + c22*c31)*(-c21*c33/2 + c23*c31/2)/pow(J, 4);

    D(4,0) = -2.0*c32*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 2.0*c32*mu/pow(J, 2) + 4.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c22*c33 + c23*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(-c11*c32/2 + c12*c31/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(-c11*c32/2 + c12*c31/2)*(-c22*c33 + c23*c32)/pow(J, 4);
    D(4,1) = 4.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c11*c33 + c13*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(-c11*c32/2 + c12*c31/2)*(c11*c33/2 - c13*c31/2)/pow(J, 4) - 4.0*mu*(-c11*c32/2 + c12*c31/2)*(-c11*c33 + c13*c31)/pow(J, 4);
    D(4,2) = 4.0*lambda*(-c11*c22 + c12*c21)*(-c11*c32/2 + c12*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(-c11*c32/2 + c12*c31/2)/pow(J, 4) - 4.0*mu*(-c11*c22 + c12*c21)*(-c11*c32/2 + c12*c31/2)/pow(J, 4);
    D(4,3) = 1.0*c31*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 1.0*c31*mu/pow(J, 2) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c12*c33/2 + c13*c32/2)/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c12*c33 - c13*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c21*c33 - c23*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(-c11*c32/2 + c12*c31/2)*(c12*c33 - c13*c32)/pow(J, 4) - 2.0*mu*(-c11*c32/2 + c12*c31/2)*(c21*c33 - c23*c31)/pow(J, 4);
    D(4,4) = -1.0*c11*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c11*mu/pow(J, 2) + 2.0*lambda*(-c11*c23/2 + c13*c21/2)*(-c11*c32/2 + c12*c31/2)/pow(J, 4) + 2.0*lambda*(c11*c23 - c13*c21)*(-c11*c32/2 + c12*c31/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*pow(-c11*c32/2 + c12*c31/2, 2)/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c11*c32 - c12*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c23 - c13*c21)*(-c11*c32/2 + c12*c31/2)/pow(J, 4) - 2.0*mu*(-c11*c32/2 + c12*c31/2)*(c11*c32 - c12*c31)/pow(J, 4);
    D(4,5) = 1.0*c12*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 1.0*c12*mu/pow(J, 2) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c12*c23 + c13*c22)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(-c21*c32 + c22*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c21*c32/2 - c22*c31/2)/pow(J, 4) - 2.0*mu*(-c11*c32/2 + c12*c31/2)*(-c12*c23 + c13*c22)/pow(J, 4) - 2.0*mu*(-c11*c32/2 + c12*c31/2)*(-c21*c32 + c22*c31)/pow(J, 4);

    D(5,0) = 4.0*lambda*(c12*c23/2 - c13*c22/2)*(-c22*c33 + c23*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c12*c23/2 - c13*c22/2)*(c22*c33/2 - c23*c32/2)/pow(J, 4) - 4.0*mu*(c12*c23/2 - c13*c22/2)*(-c22*c33 + c23*c32)/pow(J, 4);
    D(5,1) = -2.0*c13*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 2.0*c13*mu/pow(J, 2) + 4.0*lambda*(-c11*c33 + c13*c31)*(c12*c23/2 - c13*c22/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c33/2 - c13*c31/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) - 4.0*mu*(-c11*c33 + c13*c31)*(c12*c23/2 - c13*c22/2)/pow(J, 4);
    D(5,2) = 4.0*lambda*(-c11*c22 + c12*c21)*(c12*c23/2 - c13*c22/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 4.0*lambda*(c11*c22/2 - c12*c21/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) - 4.0*mu*(-c11*c22 + c12*c21)*(c12*c23/2 - c13*c22/2)/pow(J, 4);
    D(5,3) = 1.0*c23*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 1.0*c23*mu/pow(J, 2) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(-c12*c33/2 + c13*c32/2)/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(c12*c33 - c13*c32)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(-c21*c33/2 + c23*c31/2)/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(c21*c33 - c23*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c12*c23/2 - c13*c22/2)*(c12*c33 - c13*c32)/pow(J, 4) - 2.0*mu*(c12*c23/2 - c13*c22/2)*(c21*c33 - c23*c31)/pow(J, 4);
    D(5,4) = 1.0*c12*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) - 1.0*c12*mu/pow(J, 2) + 2.0*lambda*(-c11*c23/2 + c13*c21/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) + 2.0*lambda*(c11*c23 - c13*c21)*(c12*c23/2 - c13*c22/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(-c11*c32/2 + c12*c31/2)*(c12*c23/2 - c13*c22/2)/pow(J, 4) + 2.0*lambda*(c11*c32 - c12*c31)*(c12*c23/2 - c13*c22/2)*log(sqrt(pow(J, 2)))/pow(J, 4) - 2.0*mu*(c11*c23 - c13*c21)*(c12*c23/2 - c13*c22/2)/pow(J, 4) - 2.0*mu*(c11*c32 - c12*c31)*(c12*c23/2 - c13*c22/2)/pow(J, 4);
    D(5,5) = -1.0*c22*lambda*log(sqrt(pow(J, 2)))/pow(J, 2) + 1.0*c22*mu/pow(J, 2) + 2.0*lambda*(-c12*c23 + c13*c22)*(c12*c23/2 - c13*c22/2)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*pow(c12*c23/2 - c13*c22/2, 2)/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(-c21*c32 + c22*c31)*log(sqrt(pow(J, 2)))/pow(J, 4) + 2.0*lambda*(c12*c23/2 - c13*c22/2)*(c21*c32/2 - c22*c31/2)/pow(J, 4) - 2.0*mu*(-c12*c23 + c13*c22)*(c12*c23/2 - c13*c22/2)/pow(J, 4) - 2.0*mu*(c12*c23/2 - c13*c22/2)*(-c21*c32 + c22*c31)/pow(J, 4);
  }
  else
  {
    throw runtime_error("elasticity type not found");
  }
}

