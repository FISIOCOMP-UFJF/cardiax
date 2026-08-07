#include "mapping.hpp"

// ----------------------------------------------------------------------------
// Mapping implementations
// ----------------------------------------------------------------------------

Mapping::Mapping(std::vector<arma::vec3> pts, int ndim, int indx)
        : dim(ndim), points(pts), J(ndim, ndim), Jinv(ndim, ndim)
{
  // nothing
}

void Mapping::calc_jacobian(arma::mat &dshape, int nbf)
{
  int ndim = dshape.n_cols;
  int ndof = nbf ? nbf : dshape.n_rows;

  J.zeros();

  //
  // Isoparametric mapping for 2D and 3D elements
  //
  if (ndim > 1)
  {
    for (int i = 0; i < ndim; i++) {
      for (int j = 0; j < ndim; j++) {
        for (int k = 0; k < ndof; k++) {
          J(i, j) += dshape(k, j) * points[k][i];
        }
      }
    }

    // calculates inverse of the jacobian matrix
    Jinv = arma::inv(J);

    // calculates jacobian determinant
    determinant = arma::det(J);

    // check for error with negative Jacobian determinant
    //if(determinant < 0) {  cout << "NEG" << determinant << endl;  }
    //if(determinant < 0) determinant *= -1.0;
    assert(determinant > 0);
  }
  else
  {
    // TODO: this else was moved to calc_segment_jacobian !

    J.resize(1, 1);
    Jinv.resize(1, 1);

    double xl = points[0](0);
    double xr = points[1](0);
    J(0, 0) = (xr - xl) / 2.0;

    determinant = sqrt((xr - xl) * (xr - xl)); //+ (yr-yl)*(yr-yl) ) ; //+ (zr-zl)*(zr-zl) );
    double d2 = determinant * determinant;

    Jinv(0, 0) = (points[1](0) - points[0](0)) * (2.0 / d2);
    //Jinv(1,0) = ( points[1](1) - points[0](1) ) * (2.0/d2);
    //Jinv(2,0) = ( points[1](2) - points[0](2) ) * (2.0/d2);

    //Jinv = arma::inv(J);
    //determinant = arma::det(J);
  }
}

void Mapping::calc_segment_jacobian(arma::mat &dshape, int nbf)
{
  //
  // Affine mapping for 1D elements (Lagrange or Hermite Cubic)
  //
  int ndim = dshape.n_cols;
  int ndof = nbf ? nbf : dshape.n_rows;

  J.zeros();
  J.resize(3, 1);
  Jinv.resize(3, 1);

  // *** OK for 1D domain only ***
  const double xl = points[0](0);
  const double xr = points[1](0);

  const double yl = points[0](1);
  const double yr = points[1](1);

  const double zl = points[0](2);
  const double zr = points[1](2);

  J(0, 0) = (xr - xl) / 2.0;
  J(1, 0) = (yr - yl) / 2.0;
  J(2, 0) = (zr - zl) / 2.0;

  determinant = sqrt((xr - xl) * (xr - xl) + (yr - yl) * (yr - yl) + (zr - zl) * (zr - zl));

  const double d2 = determinant * determinant;
  Jinv(0, 0) = (points[1](0) - points[0](0)) * (2.0 / d2);
  Jinv(1, 0) = (points[1](1) - points[0](1)) * (2.0 / d2);
  Jinv(2, 0) = (points[1](2) - points[0](2)) * (2.0 / d2);

  //Jinv = arma::inv(J);
  //determinant = arma::det(J);
}


arma::vec3 Mapping::map_point(const arma::vec3 &x) const
{
  return arma::vec3();
}

arma::vec3 Mapping::map_point(const arma::vec3 &x,
                              const arma::vec &shape) const
{
  arma::vec3 pt;
  pt.zeros();
  for (uint i = 0; i < points.size(); i++) {
    for (int j = 0; j < dim; j++) {
      pt(j) = pt(j) + shape[i] * points[i](j);
    }
  }
  return pt;
}

void Mapping::affine(std::vector<arma::vec3> pts, int indx)
{
  // no tinymatrix 
  J.zeros();
  J(0, 0) = points[1](0) - points[0](0);
  J(1, 0) = points[1](1) - points[0](1);
  J(0, 1) = points[2](0) - points[0](0);
  J(1, 1) = points[2](1) - points[0](1);

  Jinv = arma::inv(J);
  determinant = arma::det(J);
  assert(determinant > 0);

  //calc_det_jacobian();  
  //calc_inv_jacobian();  
}

// ----------------------------------------------------------------------------
// SurfaceMapping implementations
// ----------------------------------------------------------------------------

SurfaceMapping::SurfaceMapping(std::vector<arma::vec3> pts, int ndim, int idx)
        : index(idx), dim(ndim), jmat(ndim + 1, ndim), points(pts)
{
  // unidimensional case
  if (ndim == 1) {
    surface_jac(0) = (pts[1](0) - pts[0](0));
    surface_jac(1) = (pts[1](1) - pts[0](1));

    det = sqrt(surface_jac(0) * surface_jac(0) + surface_jac(1) * surface_jac(1));

    // s = ||X-X0|| / ||X1-X0|| = (X1-X)^T (X-X) / h^2
    double hh_inv = 1.0 / (det * det);

    // jacinv_ = (points_[1]-points_[0]);
    surface_jacinv = (points[1] - points[0]) * hh_inv;

    // Normal vector
    nv(0) = -surface_jac(1);
    nv(1) = surface_jac(0);
  }
}

void SurfaceMapping::calc_jacobian(arma::mat &dshape, int nubf)
{
  // Computing surface integrals
  // Reference: O.C. Zienkiewicz
  //  "The Finite Element Method", Ch. 5, Page 148

  int ndof = nubf ? nubf : dshape.n_rows; //height();
  arma::vec3 a, b, c;

  jmat.zeros();

  for (int i = 0; i < dim + 1; i++)
    for (int j = 0; j < dim; j++)
      for (int k = 0; k < ndof; k++)
        jmat(i, j) += dshape(k, j) * points[k][i];

  if (dim == 2) {
    for (int i = 0; i < 3; i++) {
      a(i) = jmat(i, 0);
      b(i) = jmat(i, 1);
    }

    // cross-product
    c = arma::cross(a, b);
    det = norm(c, 2);
  }
  else {
    a.zeros();

    // Old
    //a(0) = jmat(0,0);
    //a(1) = jmat(1,0);
    //det = norm(a,2);

    a(0) = points[1][0] - points[0][0];
    a(1) = points[1][1] - points[0][1];
    det = norm(a, 2);

  }
}

arma::vec3 SurfaceMapping::map_point(double x) const
{
  arma::vec3 tmp = surface_jac;
  tmp *= x;
  return (points[0] + tmp);
}

// ----------------------------------------------------------------------------
// SegmentMapping implementations
// ----------------------------------------------------------------------------

/*
void SegmentMapping::calc_jacobian(arma::mat & dshape, int nbf)
{
  int ndim = dshape.n_cols;
  int ndof = nbf ? nbf : dshape.n_rows;

  J.zeros();

	//
	// Affine mapping for 1D elements (Lagrange or Hermite Cubic)
	//
	
	J.resize(2,1);
	Jinv.resize(2,1);
	
	// *** OK for 1D domain only ***
	double xl = points[0](0);
	double xr = points[1](0);
	J(0,0) = (xr-xl)/2.0;
	
	double yl = points[0](1);
	double yr = points[1](1);
	J(1,0) = (xr-xl)/2.0;
	
	//double zl = points[0](2);
	//double zr = points[1](2);
	//J(2,0) = (xr-xl)/2.0;
	
	determinant = sqrt( (xr-xl)*(xr-xl) + (yr-yl)*(yr-yl) ) ; //+ (zr-zl)*(zr-zl) );
	
	double d2 = determinant*determinant;
	Jinv(0,0) = ( points[1](0) - points[0](0) ) * (2.0/d2);
	Jinv(1,0) = ( points[1](1) - points[0](1) ) * (2.0/d2);
	//Jinv(2,0) = ( points[1](2) - points[0](2) ) * (2.0/d2);
	
	//Jinv = arma::inv(J);
	//determinant = arma::det(J);  
}

*/
