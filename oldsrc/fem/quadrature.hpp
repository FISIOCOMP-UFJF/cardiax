#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>
#include "mesh/mesh.hpp"

/**
   Base class for quadrature rules .
   TO DO: I think a singleton pattern would be nice to apply here.
   Need to test this.
*/
class Quadrature 
{
protected:

  //! Integration points 
  std::vector<arma::vec3> ipoints;

  //! Weights associated with integration points 
  std::vector<double> weights;

public:
  
  //! Default constructor
  Quadrature(){} 

  //! Constructor with integration (polynomial) order
  Quadrature(int order){}

  //! Destructor
  virtual ~Quadrature(){}
  
  //! Get quadrature point
  const arma::vec3 & get_point(int i) const { return ipoints[i]; }

  //! Get quadrature weight
  const double & get_weight(int i) const { return weights[i]; }

  //! Get the number of integration of points associated with the quadrature
  int get_num_ipoints() const { return ipoints.size(); }

  //! Create a Quadrature object based on element type and order 
  static Quadrature * create(int order, ElementType etype);
  
};

// ----------------------------------------------------------------------------

//! 1D Gauss quadrature rule for intervall (0,1). Exact up to order 2n-1
class Quadrature1d : public Quadrature
{
public:
  Quadrature1d(int order);
  ~Quadrature1d(){}
};

//! 2D quadrature rule for triangles 
class QuadratureTri : public Quadrature
{
public:
  QuadratureTri(int order);
  ~QuadratureTri(){}
};


//! 2D Gauss quadrature rules for quadrilaterals 
class QuadratureQuad : public Quadrature
{
public:
  QuadratureQuad(int order);
  ~QuadratureQuad(){}
};

//! 3D Gauss quadrature rules for tetrahedra 
class QuadratureTetra : public Quadrature
{
public:
  QuadratureTetra(int order);
  ~QuadratureTetra(){}
};

//! 3D Gauss quadrature rule for hexahedra 
class QuadratureHexa : public Quadrature
{
public:
  QuadratureHexa(int order);
  ~QuadratureHexa(){}
};

#endif
