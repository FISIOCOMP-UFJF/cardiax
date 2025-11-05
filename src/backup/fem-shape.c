#include "shape.h"

/*!
 Calculate the shape functions for isoparametric bilinear element
 @param[out] shl   Local shape functions and derivatives
 @param[out] w     Numerical integration weights
 */
void fem_iso_bilinear_quad (Real ***shl, Real *w)
{
  unsigned int i,l;
  unsigned int nint = 4;
  Real et,xi; 
  const Real qrtr = 1.0/4.0;
  const Real gaus = 0.577350269189626;
  const Real node[4][2] = {
    {-1.0, -1.0},
    { 1.0, -1.0},
    { 1.0,  1.0},
    {-1.0,  1.0}};

  for(i=0; i<4; i++){
    // weights
    w[i] = 1.0;
    for(l=0; l<nint; l++){
      xi = node[l][0] * gaus;
      et = node[l][1] * gaus;
      // shape function
      shl[0][i][l] = qrtr * (1 + node[i][0]*xi) * (1 + node[i][1]*et);
      // derivative w.r.t. xi
      shl[1][i][l] = node[i][0] * qrtr * (1 + node[i][1]*et);
      // derivative w.r.t. eta
      shl[2][i][l] = node[i][1] * qrtr * (1 + node[i][0]*xi);
    }
  }
}

/*!
 Calculate the shape functions for isoparametric linear triangle element
 @param[out] shl   Local shape functions and derivatives
 @param[out] w     Numerical integration weights
 */
void fem_iso_linear_tri (Real ***shl, Real *w)
{
  unsigned int l;
  unsigned const int nint = 3;
  // loop over integration points
  for(l=0; l<nint; l++){
    w[l] = 1.0/2.0;
    // shape function
    shl[0][0][l] = 1.0/3.0;
    shl[0][1][l] = 1.0/3.0;
    shl[0][2][l] = 1.0/3.0;
    // derivative w.r.t. r
    shl[1][0][l] = -1.0;
    shl[1][1][l] =  1.0;
    shl[1][2][l] =  0.0;
    // derivative w.r.t. s
    shl[2][0][l] = -1.0;
    shl[2][1][l] =  0.0;
    shl[2][2][l] =  1.0;
  }
}

/*!
 Calculate the shape functions for trilinear isoparametric hexahedron element
 @param[out] shl   Local shape functions and derivatives
 @param[out] w     Numerical integration weights
 */
void fem_iso_trilinear_hex (Real ***shl, Real *w)
{
  int i,l,nint;
  Real xi,et,ze;
  const Real oito = 1.0/8.0;
  const Real gaus = 0.577350269189626;
  static const Real node[8][3] = {
    {-1.0, -1.0, -1.0},
    { 1.0, -1.0, -1.0},
    { 1.0,  1.0, -1.0},
    {-1.0,  1.0, -1.0},
    {-1.0, -1.0,  1.0},
    { 1.0, -1.0,  1.0},
    { 1.0,  1.0,  1.0},
    {-1.0,  1.0,  1.0}};

  nint = 8;
 
  for(i=0; i<8; i++){
    // weights
    w[i] = 1.0;
    for(l=0; l<nint; l++){
      xi = node[l][0] * gaus;
      et = node[l][1] * gaus;
      ze = node[l][2] * gaus;
      // shape function
      shl[0][i][l] = oito*(1+node[i][0]*xi)*(1+node[i][1]*et)*(1+node[i][2]*ze);
      // derivative w.r.t. xi
      shl[1][i][l] = node[i][0]*oito*(1+node[i][1]*et)*(1+node[i][2]*ze);
      // derivative w.r.t. eta
      shl[2][i][l] = node[i][1]*oito*(1+node[i][0]*xi)*(1+node[i][2]*ze);
      // derivative w.r.t zeta
      shl[3][i][l] = node[i][2]*oito*(1+node[i][0]*xi)*(1+node[i][1]*et);
    }
  }
}
