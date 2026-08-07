#include "mixed_fe.hpp"

// -------------------------------------------------------------------
// Finite Element: VecH1FESegm
// -------------------------------------------------------------------

void VecH1FESegm::calc_shape_u(const arma::vec3 & x, arma::vec & shape) const
{
  shape.resize(2*ndim_); // CHECK
  const double xi = x(0);
  for(int i=0; i<ndim_; i++)
  {
    shape(2*i)   = 0.5*(1.0-xi);
    shape(2*i+1) = 0.5*(1.0+xi);
  }
}

void VecH1FESegm::calc_deriv_shape_u(const arma::vec3 & x,
                                     arma::mat & dshape) const
{
  dshape.resize(2*ndim_,1); // CHECK
  for(int i=0; i<ndim_; i++)
  {
    dshape(0,0)   = -0.5;
    dshape(1,0) =  0.5;
  }
}

// -------------------------------------------------------------------
// Finite Element: VecH1FETrig
// -------------------------------------------------------------------

void VecH1FETrig::calc_shape_u(const arma::vec3 & x,
			       arma::vec & shape) const
{
  shape.resize(ndofs_u_);
  shape.zeros();

  // barycentric coordinates
  double lami[3] = {1-x(0)-x(1),x(0),x(1)};

  for(int i=0; i<ndim_; i++)
  {
    shape(3*i)   = lami[0];
    shape(3*i+1) = lami[1];
    shape(3*i+2) = lami[2];
  }
}

void VecH1FETrig::calc_deriv_shape_u(const arma::vec3 & x,
				     arma::mat & dshape) const
{
  dshape.resize(ndofs_u_,2);
  dshape.zeros();

  arma::mat dlami(3,2);

  dlami(0,0) = -1.; dlami(0,1) = -1.; // 1-x-y
  dlami(1,0) =  1.; dlami(1,1) =  0;  // x
  dlami(2,0) =  0.; dlami(2,1) =  1.; // y

  for(int k=0; k<ndim_; k++)
    for(int n=0; n<3; n++)
      for(int l=0; l<2; l++)
        dshape(n+(k*3),l) = dlami(n,l);
}


// -------------------------------------------------------------------
// Finite Element: VecH1FEQuad
// -------------------------------------------------------------------

void VecH1FEQuad::calc_shape_u(const arma::vec3 & x,
			       arma::vec & shape) const
{
  const double xi  = x(0);
  const double eta = x(1);

  shape.resize(ndofs_u_);
  //shape.zeros();

  shape(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
  shape(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
  shape(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
  shape(3) = 0.25 * (1.0 - xi) * (1.0 + eta);

  shape(4) = 0.25 * (1.0 - xi) * (1.0 - eta);
  shape(5) = 0.25 * (1.0 + xi) * (1.0 - eta);
  shape(6) = 0.25 * (1.0 + xi) * (1.0 + eta);
  shape(7) = 0.25 * (1.0 - xi) * (1.0 + eta);
}

void VecH1FEQuad::calc_deriv_shape_u(const arma::vec3 & x,
				     arma::mat & dshape) const
{
  const double xi = x(0);
  const double eta = x(1);

  dshape.resize(ndofs_u_,2);
  //dshape.zeros();

  // first
  dshape(0,0) = -0.25 * (1.0 - eta);
  dshape(1,0) =  0.25 * (1.0 - eta);
  dshape(2,0) =  0.25 * (1.0 + eta);
  dshape(3,0) = -0.25 * (1.0 + eta);

  dshape(4,0) = -0.25 * (1.0 - eta);
  dshape(5,0) =  0.25 * (1.0 - eta);
  dshape(6,0) =  0.25 * (1.0 + eta);
  dshape(7,0) = -0.25 * (1.0 + eta);

  // second
  dshape(0,1) = -0.25 * (1.0 - xi);
  dshape(1,1) = -0.25 * (1.0 + xi);
  dshape(2,1) =  0.25 * (1.0 + xi);
  dshape(3,1) =  0.25 * (1.0 - xi);

  dshape(4,1) = -0.25 * (1.0 - xi);
  dshape(5,1) = -0.25 * (1.0 + xi);
  dshape(6,1) =  0.25 * (1.0 + xi);
  dshape(7,1) =  0.25 * (1.0 - xi);

}

// -------------------------------------------------------------------
// Finite Element: VecH1FETetra
// -------------------------------------------------------------------

void VecH1FETetra::calc_shape_u(const arma::vec3 & x,
				arma::vec & shape) const
{
  shape.resize(ndofs_u_);
  //shape.zeros();

  double lami[4] = {1.0-x(0)-x(1)-x(2), x(0), x(1), x(2)};

  for(int i=0; i<ndofs_u_; i++)
    shape(i) = lami[i%4];

}

void VecH1FETetra::calc_deriv_shape_u(const arma::vec3 & x,
				      arma::mat & dshape) const
{
  dshape.resize(ndofs_u_,ndim_);
  dshape.zeros();

  arma::mat dlami(4,3);

  dlami(0,0) = -1.; dlami(0,1) = -1.; dlami(0,2) = -1.; // 1-x-y-z
  dlami(1,0) =  1.; dlami(1,1) =  0.; dlami(1,2) =  0.; // x
  dlami(2,0) =  0.; dlami(2,1) =  1.; dlami(2,2) =  0.; // y
  dlami(3,0) =  0.; dlami(3,1) =  0.; dlami(3,2) =  1.; // z

  for(int n=0; n<ndofs_u_; n++)
    for(int l=0; l<ndim_; l++)
      dshape(n,l) = dlami(n%4,l);

}

// -------------------------------------------------------------------
// Finite Element: VecH1FEHexa
// -------------------------------------------------------------------

void VecH1FEHexa::calc_shape_u(const arma::vec3 & x, 
			       arma::vec & shape) const
{
  const double xi  = x(0);
  const double eta = x(1);
  const double psi = x(2);

  shape.resize(ndofs_u_);
  shape.zeros();

  shape(0) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - psi);
  shape(1) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - psi);
  shape(2) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - psi);
  shape(3) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - psi);
  shape(4) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + psi);
  shape(5) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + psi);
  shape(6) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + psi);
  shape(7) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + psi);

  shape( 8) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - psi);
  shape( 9) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - psi);
  shape(10) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - psi);
  shape(11) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - psi);
  shape(12) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + psi);
  shape(13) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + psi);
  shape(14) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + psi);
  shape(15) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + psi);

  shape(16) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - psi);
  shape(17) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - psi);
  shape(18) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - psi);
  shape(19) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - psi);
  shape(20) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + psi);
  shape(21) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + psi);
  shape(22) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + psi);
  shape(23) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + psi);
}

void VecH1FEHexa::calc_deriv_shape_u(const arma::vec3 & x,
				     arma::mat & dshape) const
{
  const double xi = x(0);
  const double eta = x(1);
  const double psi = x(2);

  dshape.resize(ndofs_u_, ndim_);
  dshape.zeros();

  // 1
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

  // 2
  dshape( 8,0) = -0.125 * (1.0 - eta) * (1.0 - psi);
  dshape( 9,0) =  0.125 * (1.0 - eta) * (1.0 - psi);
  dshape(10,0) =  0.125 * (1.0 + eta) * (1.0 - psi);
  dshape(11,0) = -0.125 * (1.0 + eta) * (1.0 - psi);
  dshape(12,0) = -0.125 * (1.0 - eta) * (1.0 + psi);
  dshape(13,0) =  0.125 * (1.0 - eta) * (1.0 + psi);
  dshape(14,0) =  0.125 * (1.0 + eta) * (1.0 + psi);
  dshape(15,0) = -0.125 * (1.0 + eta) * (1.0 + psi);

  dshape( 8,1) = -0.125 * (1.0 - xi) * (1.0 - psi);
  dshape( 9,1) = -0.125 * (1.0 + xi) * (1.0 - psi);
  dshape(10,1) =  0.125 * (1.0 + xi) * (1.0 - psi);
  dshape(11,1) =  0.125 * (1.0 - xi) * (1.0 - psi);
  dshape(12,1) = -0.125 * (1.0 - xi) * (1.0 + psi);
  dshape(13,1) = -0.125 * (1.0 + xi) * (1.0 + psi);
  dshape(14,1) =  0.125 * (1.0 + xi) * (1.0 + psi);
  dshape(15,1) =  0.125 * (1.0 - xi) * (1.0 + psi);

  dshape( 8,2) = -0.125 * (1.0 - xi) * (1.0 - eta);
  dshape( 9,2) = -0.125 * (1.0 + xi) * (1.0 - eta);
  dshape(10,2) = -0.125 * (1.0 + xi) * (1.0 + eta);
  dshape(11,2) = -0.125 * (1.0 - xi) * (1.0 + eta);
  dshape(12,2) =  0.125 * (1.0 - xi) * (1.0 - eta);
  dshape(13,2) =  0.125 * (1.0 + xi) * (1.0 - eta);
  dshape(14,2) =  0.125 * (1.0 + xi) * (1.0 + eta);
  dshape(15,2) =  0.125 * (1.0 - xi) * (1.0 + eta);

  // 3
  dshape(16,0) = -0.125 * (1.0 - eta) * (1.0 - psi);
  dshape(17,0) =  0.125 * (1.0 - eta) * (1.0 - psi);
  dshape(18,0) =  0.125 * (1.0 + eta) * (1.0 - psi);
  dshape(19,0) = -0.125 * (1.0 + eta) * (1.0 - psi);
  dshape(20,0) = -0.125 * (1.0 - eta) * (1.0 + psi);
  dshape(21,0) =  0.125 * (1.0 - eta) * (1.0 + psi);
  dshape(22,0) =  0.125 * (1.0 + eta) * (1.0 + psi);
  dshape(23,0) = -0.125 * (1.0 + eta) * (1.0 + psi);

  dshape(16,1) = -0.125 * (1.0 - xi) * (1.0 - psi);
  dshape(17,1) = -0.125 * (1.0 + xi) * (1.0 - psi);
  dshape(18,1) =  0.125 * (1.0 + xi) * (1.0 - psi);
  dshape(19,1) =  0.125 * (1.0 - xi) * (1.0 - psi);
  dshape(20,1) = -0.125 * (1.0 - xi) * (1.0 + psi);
  dshape(21,1) = -0.125 * (1.0 + xi) * (1.0 + psi);
  dshape(22,1) =  0.125 * (1.0 + xi) * (1.0 + psi);
  dshape(23,1) =  0.125 * (1.0 - xi) * (1.0 + psi);

  dshape(16,2) = -0.125 * (1.0 - xi) * (1.0 - eta);
  dshape(17,2) = -0.125 * (1.0 + xi) * (1.0 - eta);
  dshape(18,2) = -0.125 * (1.0 + xi) * (1.0 + eta);
  dshape(19,2) = -0.125 * (1.0 - xi) * (1.0 + eta);
  dshape(20,2) =  0.125 * (1.0 - xi) * (1.0 - eta);
  dshape(21,2) =  0.125 * (1.0 + xi) * (1.0 - eta);
  dshape(22,2) =  0.125 * (1.0 + xi) * (1.0 + eta);
  dshape(23,2) =  0.125 * (1.0 - xi) * (1.0 + eta);

}
