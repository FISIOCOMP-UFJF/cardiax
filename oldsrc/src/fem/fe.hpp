#ifndef FE_HPP
#define FE_HPP

#include <vector>
#include "fem/element.hpp"
#include "mesh/mesh.hpp"

// Forward declaration of Mesh
//class Mesh;

/**  This is an abstract base class for a finite element and provides 
     interface to compute shape functions and its derivatives.
*/
class FiniteElement
{
protected:
  
  //! Number of degrees of freedom 
  int ndof;

  //! Number of dimensions 
  int ndim;

  //! Geometric element type information 
  ElementType etype;

  //! Ref to mesh in order to get local tensor of conductivities and more 
  const Mesh * mesh;
  
public:
  
  //! Constructor 
  FiniteElement(int ndofs, int dim, ElementType etype, const Mesh * msh)
    : ndof(ndofs), ndim(dim), etype(etype), mesh(msh)
  {
    // do nothing
  }

  //! Virtual destructor 
  virtual ~FiniteElement(){}

  //! Get number of degrees of freedom 
  int ndofs() const { return ndof; }

  //! Get number of degrees of freedom 
  int get_ndof() const { return ndof; }

  //! Get number of spatial dimensions 
  int get_ndim() const { return ndim; }

  //! The element type information 
  ElementType type() const { return etype; } 

  //! Map the current element to the reference element 
  Mapping get_mapping(int i) const
  {
    return mesh->get_element_mapping(i);
  }

  //! Map the current element to the reference element 
  SurfaceMapping get_boundary_mapping(int i) const
  {
   return mesh->get_boundary_element_mapping(i);
  }
		
  //! Compute shape function in point x on the reference element 
  virtual void calc_shape(const arma::vec3 & x, arma::vec & shape) const=0;

  /*! Evaluates the shape function derivative at x on
      the reference element. Returns a matrix.
      Matrix entries: dshape[i][j] = d_phi_i / d_xj (x)
  */
  virtual void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const=0;

  //! Polynomial order of the element 
  virtual int order() const=0;

};

// -----------------------------------------------------------------------------

//! One-dimensional line segment element. 
class FESegm : public FiniteElement
{
public:

  FESegm(const Mesh * msh): FiniteElement(2, 1, ELEM_SEGM, msh) {}
  void calc_shape(const arma::vec3 & x, arma::vec & shape) const;
  //void calc_deriv_shape(const arma::vec & x, arma::vec & dshape) const;
  void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const;
  int order() const {return 1;}
};

// -----------------------------------------------------------------------------

//! One-dimensional cubic Hermite element
class FECubicHermite : public FiniteElement
{
public:

  FECubicHermite(const Mesh * msh): FiniteElement(4, 1, ELEM_SEGM, msh) {}
  void calc_shape(const arma::vec3 & x, arma::vec & shape) const;
  //void calc_deriv_shape(const arma::vec & x, arma::vec & dshape) const;
  void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const;
  int order() const { return 3; }
};

// -----------------------------------------------------------------------------

//! Two-dimensional triangular element. 
class FETrig : public FiniteElement
{   
public:

  FETrig(const Mesh * msh): FiniteElement(3, 2, ELEM_TRIG, msh) {}
  void calc_shape(const arma::vec3 & x, arma::vec & shape) const;
  void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const;  
  int order() const {return 1;}
};

//! Three-dimensional tetrahedral element 
class FETetra : public FiniteElement
{
public:
  FETetra(const Mesh * msh): FiniteElement(4, 3, ELEM_TETRA, msh) {}
  void calc_shape(const arma::vec3 & x, arma::vec & shape) const;
  void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const; 
  int order() const {return 1;}
};

// -----------------------------------------------------------------------------

//! Two-dimensional quadrilateral element 
class FEQuad : public FiniteElement
{   
public: 
  FEQuad(const Mesh * msh): FiniteElement(4, 2, ELEM_QUAD, msh) {}
  void calc_shape(const arma::vec3 & x, arma::vec & shape) const;
  void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const; 
  int order() const {return 1;}
};

// -----------------------------------------------------------------------------

//! Three-dimensional hexahedal element 
class FEHexa : public FiniteElement
{
public:
  FEHexa(const Mesh * msh): FiniteElement(8, 3, ELEM_HEXA, msh) {}
  void calc_shape(const arma::vec3 & x, arma::vec & shape) const;
  void calc_deriv_shape(const arma::vec3 & x, arma::mat & dshape) const; 
  int order() const {return 1;}
};

#endif
