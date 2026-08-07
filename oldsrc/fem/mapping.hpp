#ifndef MAPPING_HPP
#define MAPPING_HPP

#include <iostream>
#include <vector>
#include <armadillo>
#include "linalg/linalg.h"

// Forward declarations
class FiniteElement;

/** 
    Isoparametric element mapping class
*/
class Mapping
{
public:

  //! Constructor
  Mapping(std::vector<arma::vec3> points, int dim, int index = 1);

  //! Constructor
  Mapping(std::vector<arma::vec3> points, int index = 1);

  //! Constructor based on Finite Element
  Mapping(FiniteElement& fe, std::vector<arma::vec3> points, int index);

  //! Compute affine mapping
  void affine(std::vector<arma::vec3> pts, int index);

  //! Compute isoparametric mapping
  void calc_jacobian(arma::mat & dshape, int nubf=0);

	//! Compute isoparametric mapping for Segment
	void calc_segment_jacobian(arma::mat & dshape, int nubf=0);

  //! Subdomain index of Element 
  int get_index() const { return index; }

  //! Global Points Coordinates
  const std::vector<arma::vec3> & get_points() const { return points; };

  //! Jacobian dx/dX 
  const arma::mat& get_jacobian() const { return J; }

  //! InverseJacobian dx/dX 
  const arma::mat& get_inv_jacobian() const { return Jinv; }

  //! Deterimant of Jacobian
  double get_det_jacobian() const { return determinant; }

  //! Map local point x to physical point X
  arma::vec3 map_point(const arma::vec3& x) const;

  //! Map local point x to physical point X
  arma::vec3 map_point(const arma::vec3& x, const arma::vec& shape) const;
   
protected:

  //! Dimension
  int dim;

  //! Subdomain Index
  int index;

  //! Jacobian determinant
  double determinant;
  
  //! Global points of element
  std::vector<arma::vec3> points;

  //! Jacobian matrix
  arma::mat J;
  
  //! Inverse Jacobian matrix
  arma::mat Jinv;

};

// ----------------------------------------------------------------------------

/** Surface element mapping class
    See: "Spectral/hp element methods for CFD"
          George Em Karniadakis, Spencer Sherwin
          Section 4.1.4 "Discrete evaluation of the surface Jacobian"
*/
class SurfaceMapping
{
private:
  int index;
  int dim;
  double det;
  
  arma::mat jmat;
  arma::vec3 surface_jac;
  arma::vec3 surface_jacinv;
  arma::vec3 nv;
  std::vector<arma::vec3> points;

public:

  SurfaceMapping(std::vector<arma::vec3> points, int ndim, int index = 1);

  void calc_jacobian(arma::mat & dshape, int nubf=0);

  int get_index() const
  {
    return index;
  }

  const std::vector<arma::vec3>& get_points() const
  {
    return points;
  }

  const arma::vec3 & get_jacobian(const Point1d& p) const
  {
    return surface_jac;
  }

  const arma::vec3 & get_inv_jacobian(const arma::vec3& p) const
  {
    return surface_jacinv;
  }

  double get_det_jacobian() const
  {
    return det;
  }

  arma::vec3 map_point(double x) const;
};


// ----------------------------------------------------------------------------

/*
class SegmentMapping : public Mapping
{
 public:

  //! Overrides the computation of the isoparametric mapping
  void calc_jacobian(arma::mat & dshape, int nubf=0);
	
};
*/

#endif
