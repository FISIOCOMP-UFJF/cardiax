#include "material_data.hpp"

MaterialData::MaterialData(const Element & el, const arma::mat33 & ftens)
  : F(ftens)
{
  J = arma::det(ftens);

  f = el.get_fiber();
  s = el.get_trans();
  n = el.get_normal();
  marker = el.get_index();
}

arma::mat33 MaterialData::left_cauchy_green()
{
  arma::mat33 b;
  b(0,0) = F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(0,2)*F(0,2);
  b(0,1) = F(0,0)*F(1,0) + F(0,1)*F(1,1) + F(0,2)*F(1,2);
  b(0,2) = F(0,0)*F(2,0) + F(0,1)*F(2,1) + F(0,2)*F(2,2);
  b(1,0) = b(0,1);
  b(1,1) = F(1,0)*F(1,0) + F(1,1)*F(1,1) + F(1,2)*F(1,2);
  b(1,2) = F(1,0)*F(2,0) + F(1,1)*F(2,1) + F(1,2)*F(2,2);
  b(2,0) = b(0,2);
  b(2,1) = b(1,2);
  b(2,2) = F(2,0)*F(2,0) + F(2,1)*F(2,1) + F(2,2)*F(2,2);
  return b;
}

arma::mat33 MaterialData::right_cauchy_green()
{
  arma::mat33 C;
  C(0,0) = F(0,0)*F(0,0) + F(1,0)*F(1,0) + F(2,0)*F(2,0);
  C(0,1) = F(0,0)*F(0,1) + F(1,0)*F(1,1) + F(2,0)*F(2,1);
  C(0,2) = F(0,0)*F(0,2) + F(1,0)*F(1,2) + F(2,0)*F(2,2);
  C(1,0) = C(0,1);
  C(1,1) = F(0,1)*F(0,1) + F(1,1)*F(1,1) + F(2,1)*F(2,1);
  C(1,2) = F(0,1)*F(0,2) + F(1,1)*F(1,2) + F(2,1)*F(2,2);
  C(2,0) = C(0,2);
  C(2,1) = C(1,2);
  C(2,2) = F(0,2)*F(0,2) + F(1,2)*F(1,2) + F(2,2)*F(2,2);
  return C;
}

arma::mat33 MaterialData::isochoric_def_grad()
{  
  return pow(J, (-1./3.)) * F;
}

arma::mat33 MaterialData::isochoric_lcg()
{  
  arma::mat33 b = left_cauchy_green();
  return pow(J, -(2.0/3.0)) * b;
}

arma::mat33 MaterialData::isochoric_rcg()
{
  arma::mat33 C = right_cauchy_green();
  return pow(J, -(2.0/3.0)) * C;
}

arma::mat33 MaterialData::lagrangian_strain()
{
  return 0.5*((F.t()*F) - arma::eye(3,3));
}

