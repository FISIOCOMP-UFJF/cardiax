#ifndef MIXED_FE_HPP_
#define MIXED_FE_HPP_

#include <vector>
#include <cmath>
#include "mesh/mesh.hpp"

/** This is an abstract base class for Mixed finite elements,
    which can also be used for vector-value problems, such as
    linear elasticity and nonlinear elasticity.
*/

class MixedFiniteElement
{
protected:

  //! Number of spatial dimensions of the element 
  int ndim_;

  //! Number of nodes
  int nnode_;

  //! Number of dofs of variable U 
  int ndofs_u_;

  //! Number of dofs of variable P 
  int ndofs_p_;

  //! Polynomial order of U 
  int order_u_;

  //! Polynomial order of P 
  int order_p_;

  //! Geometric element type 
  ElementType etype;

  //! Constant reference to mesh 
  const Mesh & mesh;

public:

  //! Default constructor 
  MixedFiniteElement(int dim, ElementType type, const Mesh & msh) 
    : ndim_(dim), ndofs_u_(0), ndofs_p_(0), 
      order_u_(-1), order_p_(-1), etype(type), mesh(msh)
    {};

  //! Virtual destructor 
  virtual ~MixedFiniteElement(){};

  //! Get number of spatial dimensions 
  int get_ndim() const { return ndim_; }

  //! Get number of nodes
  int get_nnode() const { return nnode_; }

  //! Get number of total dofs 
  int get_ndofs() const {return ndofs_u_ + ndofs_p_;};

  //! Get number of U variable dofs 
  int get_ndofs_u() const {return ndofs_u_;};

  //! Get number of P variable dofs 
  int get_ndofs_p() const {return ndofs_p_;};

  //! The polynomial order of U 
  int get_order_u() const { return order_u_; }

  //! The polynomial order of P 
  int get_order_p() const { return order_p_; }

  /** The type of the geometric element 
      (triangle, quadrilateral, tetrahedron, hexahedron).
  */
  ElementType get_type() const { return etype; }

  //! Get a mapping to the reference element of the element i 
  Mapping get_mapping(int i) const
  {
    return mesh.get_element_mapping(i);
  }

  //! Get a mapping to the reference element of the element i 
  Mapping get_mapping(int i, std::vector<arma::vec3> & x) const
  {
    return mesh.get_element_mapping(i, x);
  }

  Mapping * get_mapping_ptr(int i, std::vector<arma::vec3> & x) const
  {
    return mesh.element_mapping(i,x);
  }

  //! Get a mapping to the reference element of the boundary element i
  SurfaceMapping get_boundary_mapping(int i) const
  {
    return mesh.get_boundary_element_mapping(i);
  }

  //! Get a mapping to the reference element of the boundary element i.
  SurfaceMapping get_boundary_mapping(int i, std::vector<arma::vec3> & x) const
  {
    return mesh.get_boundary_element_mapping(i,x);
  }

  /** Interface to a member function that computes the value
      of the shape functions for U at point x. Implement this
      in a derived class.
   */
  virtual void calc_shape_u(const arma::vec3 & x, arma::vec & shape) const = 0;

  /** Interface to a member function that computes the value
      of the shape functions for P at point x. Implement this
      in a derived class.
   */
  virtual void calc_shape_p(const arma::vec3 & x, arma::vec & shape) const = 0;

  /** Interface to a member function that computes the value of the
      derivatives of the shape functions for variabel U
      and store them in a matrix:
      
      dshape(nb, ndim)
      
      | dN1/dxi , dN1/deta, dN1/dzeta |
      | dN2/dxi , dN2/deta, dN2/dzeta |
      |  .....  , .....   , .....     |
      
      Implement this in a derived class.
   */
  virtual void calc_deriv_shape_u(const arma::vec3 & x, arma::mat & dshape) const = 0;

};

// ----------------------------------------------------------------------------

/** Concrete finite element implementation providing vector-valued linear
    U = (u1,u2) on line elements segments.
*/
class VecH1FESegm : public MixedFiniteElement
{
public:
  VecH1FESegm(int dim, const Mesh & msh) 
    : MixedFiniteElement(dim, ELEM_SEGM, msh)
    { 
      nnode_   = 2;
      ndofs_u_ = 2*dim;
      order_u_ = 1;
      ndofs_p_ = 0;
      order_p_ = -1;
  }

  virtual void calc_shape_u(const arma::vec3 & x, arma::vec & shape) const;
  virtual void calc_shape_p(const arma::vec3 & x, arma::vec & shape) const {};
  virtual void calc_deriv_shape_u(const arma::vec3 & x, arma::mat & dshape) const;
};


/** Concrete Finite Element providing vector-valued  linear 
    U = (u1,u2) on triangles
*/
class VecH1FETrig : public MixedFiniteElement
{
public:
  VecH1FETrig(int dim, const Mesh & msh) 
    : MixedFiniteElement(dim, ELEM_TRIG, msh)
  { 
    nnode_   = 3;
    ndofs_u_ = 3*dim;
    order_u_ = 1;
    ndofs_p_ = 0;
    order_p_ = -1;
  }

  virtual void calc_shape_u(const arma::vec3 & x, arma::vec & shape) const;
  virtual void calc_shape_p(const arma::vec3 & x, arma::vec & shape) const {};
  virtual void calc_deriv_shape_u(const arma::vec3 & x, arma::mat & dshape) const;
};

/** Concrete Finite Element providing vector-valued linear 
    U = (u1,u2) on Quadrilaterals
*/
class VecH1FEQuad : public MixedFiniteElement
{
public:
  VecH1FEQuad(int dim, const Mesh & msh) 
    : MixedFiniteElement(dim, ELEM_QUAD, msh)
  { 
    nnode_   = 4;
    ndofs_u_ = 4*dim;
    order_u_ = 1;
    ndofs_p_ = 0;
    order_p_ = -1;
  }

  virtual void calc_shape_u(const arma::vec3 & x, arma::vec & shape) const;
  virtual void calc_shape_p(const arma::vec3 & x, arma::vec & shape) const {};
  virtual void calc_deriv_shape_u(const arma::vec3 & x, arma::mat & dshape) const;
};


/** Concrete Finite Element providing vector-valued linear 
    U = (u1,u2) on Tetrahedral
*/
class VecH1FETetra : public MixedFiniteElement
{
public:
  VecH1FETetra(int dim, const Mesh & msh) 
    : MixedFiniteElement(dim, ELEM_TETRA, msh)
  { 
    nnode_   = 4;
    ndofs_u_ = 4*dim;
    order_u_ = 1;
    ndofs_p_ = 0;
    order_p_ = -1;
  }

  virtual void calc_shape_u(const arma::vec3 & x, arma::vec & shape) const;
  virtual void calc_shape_p(const arma::vec3 & x, arma::vec & shape) const {};
  virtual void calc_deriv_shape_u(const arma::vec3 & x, arma::mat & dshape) const;
};

/** Concrete Finite Element providing vector-valued linear 
    U = (u1,u2) on Hexahedron
*/
class VecH1FEHexa : public MixedFiniteElement
{
public:
  VecH1FEHexa(int dim, const Mesh & msh) 
    : MixedFiniteElement(dim, ELEM_HEXA, msh)
  { 
    nnode_   = 8;
    ndofs_u_ = 8*dim;
    order_u_ = 1;
    ndofs_p_ = 0; 
    order_p_ = -1;
  }

  virtual void calc_shape_u(const arma::vec3 & x, arma::vec & shape) const;
  virtual void calc_shape_p(const arma::vec3 & x, arma::vec & shape) const {};
  virtual void calc_deriv_shape_u(const arma::vec3 & x, arma::mat & dshape) const;
};


// Some typedefs
typedef MixedFiniteElement MxFE;

#endif
