#include <cmath>
#include "fe.hpp"

// ----------------------------------------------------------------------------
// Finite Element: FESegm
// ----------------------------------------------------------------------------

void FESegm::calc_shape(const arma::vec3 & x, arma::vec & shape) const
{
  shape.resize(ndof);

  // barycentric coordinates
  //double lami[2] = {1-x[0],x[0]};

  const double xi = x(0);

  // this is valid for order = 1
  shape(0)=0.5*(1.0-xi);
  shape(1)=0.5*(1.0+xi);
}

void FESegm::calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const
{
  dshape.resize(ndofs(),1);

  dshape(0,0) = -0.5;
  dshape(1,0) =  0.5;
}

// ----------------------------------------------------------------------------
// Finite Element: FECubicHermite
// ----------------------------------------------------------------------------

void FECubicHermite::calc_shape(const arma::vec3 & x, arma::vec & shape) const
{
  const double xi  = x(0);
  const double xi2 = xi*xi;
  const double xi3 = xi*xi*xi;
  shape.resize(ndofs());  
  shape(0) = 0.5 - 0.75*xi + 0.25*xi3;
  shape(1) = 0.5 + 0.75*xi - 0.25*xi3;
  shape(2) =  0.25 - 0.25*xi - 0.25*xi2 + 0.25*xi3;
  shape(3) = -0.25 - 0.25*xi + 0.25*xi2 + 0.25*xi3;
}

void FECubicHermite::calc_deriv_shape(const arma::vec3 & x,
																			arma::mat & dshape) const
{
  const double xi  = x(0);
  const double xi2 = xi*xi;
  dshape.resize(ndofs(),1);
  dshape(0) = - 0.75 + 0.75*xi2;
  dshape(1) = + 0.75 - 0.75*xi2;
  dshape(2) = - 0.25 - 0.5*xi + 0.75*xi2;
  dshape(3) = - 0.25 + 0.5*xi + 0.75*xi2;
}

// ----------------------------------------------------------------------------
// Finite Element: FETrig
// ----------------------------------------------------------------------------

void FETrig::calc_shape(const arma::vec3 & x, arma::vec & shape) const
{
  shape.resize(ndof);
  
  // barycentric coordinates
  double lami[3] = {1-x(0)-x(1),x(0),x(1)};
  
  shape(0)=lami[0];
  shape(1)=lami[1];
  shape(2)=lami[2];
}

// dshape[i][j] ... d_phi_i/d_xj 
void FETrig::calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const
{
  dshape.resize(ndofs(),2);
  dshape.zeros();
  
  TinyVec<double,2> dlami[3]; 

  dlami[0](0) = -1.; dlami[0](1) = -1.; // 1-x-y
  dlami[1](0) =  1.; dlami[1](1) =  0;  // x
  dlami[2](0) =  0.; dlami[2](1) =  1.; // y
  
  for(int n=0;n<3;n++) 
    for(int l=0;l<2;l++) 
      dshape(n,l) = dlami[n](l);

} 

// ----------------------------------------------------------------------------
// Finite Element: FETetra
// ----------------------------------------------------------------------------

void FETetra::calc_shape(const arma::vec3 & x, arma::vec & shape) const
{
  shape.resize(ndof);

  double lami[4] = {1.0-x(0)-x(1)-x(2), x(0), x(1), x(2)};
  for(int i=0; i<ndof; i++)
    shape(i) = lami[i];
}

void FETetra::calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const
{
  dshape.resize(ndof,ndim);

  TinyVec<double,3> dlami[4];
  
  dlami[0](0) = -1.; dlami[0](1) = -1.; dlami[0](2) = -1.; // 1-x-y-z
  dlami[1](0) =  1.; dlami[1](1) =  0.; dlami[1](2) =  0.; // x
  dlami[2](0) =  0.; dlami[2](1) =  1.; dlami[2](2) =  0.; // y
  dlami[3](0) =  0.; dlami[3](1) =  0.; dlami[3](2) =  1.; // z

  for(int n=0; n<ndof; n++)
    for(int l=0; l<ndim; l++)
      dshape(n,l) = dlami[n](l);
}

// ----------------------------------------------------------------------------
// Finite Element: FEQuad
// ----------------------------------------------------------------------------

void FEQuad::calc_shape(const arma::vec3 & x, arma::vec & shape) const
{
  const double xi = x(0);
  const double eta = x(1);

  shape.resize(ndofs());

  shape(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
  shape(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
  shape(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
  shape(3) = 0.25 * (1.0 - xi) * (1.0 + eta);
}

void FEQuad::calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const
{
  const double xi = x(0);
  const double eta = x(1);

  dshape.resize(ndofs(),2);

  dshape(0,0) = -0.25 * (1.0 - eta);
  dshape(1,0) =  0.25 * (1.0 - eta);
  dshape(2,0) =  0.25 * (1.0 + eta);
  dshape(3,0) = -0.25 * (1.0 + eta);

  dshape(0,1) = -0.25 * (1.0 - xi);
  dshape(1,1) = -0.25 * (1.0 + xi);
  dshape(2,1) =  0.25 * (1.0 + xi);
  dshape(3,1) =  0.25 * (1.0 - xi);
} 

// ----------------------------------------------------------------------------
// Finite Element: FEHexa
// ----------------------------------------------------------------------------

void FEHexa::calc_shape(const arma::vec3 & x, arma::vec & shape) const
{
  const double xi  = x(0);
  const double eta = x(1);
  const double psi = x(2);

  shape.resize(ndof);

  shape(0) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - psi);
  shape(1) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - psi);
  shape(2) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - psi);
  shape(3) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - psi);
  shape(4) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + psi);
  shape(5) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + psi);
  shape(6) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + psi);
  shape(7) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + psi);
}

void FEHexa::calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const
{
  const double xi  = x(0);
  const double eta = x(1);
  const double psi = x(2);

  dshape.resize(ndof,ndim);
  dshape.zeros();

  dshape(0,0) = -0.125 * (1.0 - eta) * (1.0 - psi);
  dshape(1,0) =  0.125 * (1.0 - eta) * (1.0 - psi);
  dshape(2,0) =  0.125 * (1.0 + eta) * (1.0 - psi);
  dshape(3,0) = -0.125 * (1.0 + eta) * (1.0 - psi);
  dshape(4,0) = -0.125 * (1.0 - eta) * (1.0 + psi);
  dshape(5,0) =  0.125 * (1.0 - eta) * (1.0 + psi);
  dshape(6,0) =  0.125 * (1.0 + eta) * (1.0 + psi);
  dshape(7,0) = -0.125 * (1.0 + eta) * (1.0 + psi);

  dshape(0,1) = -0.125 * (1.0 - xi) * (1.0 - psi);
  dshape(1,1) = -0.125 * (1.0 + xi) * (1.0 - psi);
  dshape(2,1) =  0.125 * (1.0 + xi) * (1.0 - psi);
  dshape(3,1) =  0.125 * (1.0 - xi) * (1.0 - psi);
  dshape(4,1) = -0.125 * (1.0 - xi) * (1.0 + psi);
  dshape(5,1) = -0.125 * (1.0 + xi) * (1.0 + psi);
  dshape(6,1) =  0.125 * (1.0 + xi) * (1.0 + psi);
  dshape(7,1) =  0.125 * (1.0 - xi) * (1.0 + psi);

  dshape(0,2) = -0.125 * (1.0 - xi) * (1.0 - eta);
  dshape(1,2) = -0.125 * (1.0 + xi) * (1.0 - eta);
  dshape(2,2) = -0.125 * (1.0 + xi) * (1.0 + eta);
  dshape(3,2) = -0.125 * (1.0 - xi) * (1.0 + eta);
  dshape(4,2) =  0.125 * (1.0 - xi) * (1.0 - eta);
  dshape(5,2) =  0.125 * (1.0 + xi) * (1.0 - eta);
  dshape(6,2) =  0.125 * (1.0 + xi) * (1.0 + eta);
  dshape(7,2) =  0.125 * (1.0 - xi) * (1.0 + eta);
}
