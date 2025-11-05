#ifndef MIXED_FESPACE_HPP
#define MIXED_FESPACE_HPP

#include <iostream>
#include <vector>
#include "mixed_fe.hpp"
#include "mesh/mesh.hpp"

/** 
    Mixed FE-Space with two variables u,p for mixed problem
*/
class MixedFESpace
{
public:

  //! Constructor
  MixedFESpace();

  //! Constructor with mesh as parameter
  MixedFESpace(const Mesh * mesh);

  //! Destructor
  virtual ~MixedFESpace(){};

  //! Get total number of dofs
  int get_ndofs() const{ return ndofs_; }

  //! Get number of dofs for variable u
  int get_ndofs_u() const{ return ndofs_u_; }

  //! Get number of dofs for variable p
  int get_ndofs_p() const{ return ndofs_p_; }

  //! Get global dof numbers of u-variable of element i
  virtual int get_element_dofs_u(int i, std::vector<int> & dnums) const
  {
    cout << " MixedFESpace::get_element_dofs_u() is purely virtual\n";
    exit(0);
  }

  //! Get global dof numbers of p-variable of element i
  virtual int get_element_dofs_p(int i, std::vector<int> & dnums) const
  {
    cout << " MixedFESpace::get_element_dofs_p() is purely virtual\n";
    exit(0);
  }

  //! Get global boundary dof numbers of u-variable of element i
  virtual int get_boundary_element_dofs_u(int i, std::vector<int> & dnums) const
  {
    cout << " MixedFESpace::get_boundary_element_dofs_u() is purely virtual\n";
    exit(0);
  }

  //! Get global boundary dof numbers of p-variable of element i
  virtual int get_boundary_element_dofs_p(int i, std::vector<int> & dnums) const
  {
    cout << " MixedFESpace::get_boundary_element_dofs_p() is purely virtual\n";
    exit(0);
  }

  //! Interface to create a mixed finite element object
  virtual MixedFiniteElement * createFE() const = 0;

  //! Interface to create a mixed finite element object.
  virtual MixedFiniteElement * create_boundary_FE() const = 0;

  //! Default dof distribution
  virtual void setup_dofs() = 0;

protected:

  //! Computational mesh
  const Mesh * mesh_;

  //! Number of spatial dimensions of the problem
  int ndim_;

  //! Number of total degrees of freedom
  int ndofs_;

  //! Number of degrees of freedom for U
  int ndofs_u_;

  //! Number of degrees of freedom for P
  int ndofs_p_;

  //! The dofs stored in an nelem x nnode matrix
  arma::umat dofs;

};

/** This class implements a dummy space having vector-valued H^1
    for U and nothing for P. Continuous, piecewise linear for U.
    (Old) Pre-state for Mini-element.
 */
class VecH1FESpace: public MixedFESpace
{
protected:

  // Global dof number of first dof of second component -> np 
  int first_u2_;

  // Global dof number of first dof of third component -> 2*np
  int first_u3_;

  // TODO: use a std::vector<int> firsts
  // first u1 -> 0
  // first u2 -> ndof
  // first u3 -> 2*ndof

public:

  //! Default constructor
  VecH1FESpace(): MixedFESpace() {}

  //! Constructor with mesh as parameter
  VecH1FESpace(const Mesh * mesh): MixedFESpace(mesh) {}

  //! Destructor
  ~VecH1FESpace(){};

  /*! Set mesh member variable to msh and setup some properties
      (number of spatial dimensions, ndofs, first_u2 and so on)
      according to the mesh variable.
   */
  void set_mesh(const Mesh * msh)
  {
    mesh_  = msh;

    ndim_     = msh->get_n_dim();
    ndofs_    = msh->get_n_points();
    ndofs_u_  = ndim_ * msh->get_n_points(); // Vector-Valued p.w. linears
    first_u2_ = msh->get_n_points();
    first_u3_ = 2 * msh->get_n_points();
    ndofs_p_  = 0;
    ndofs_    = ndofs_u_ + ndofs_p_;
  }

  //! Get element dof numbers for U
  virtual int get_element_dofs_u(int i, std::vector<int> & dnums) const;

  //! Get element boundary dof numbers for U
  virtual int get_boundary_element_dofs_u(int i, std::vector<int> & dnums) const;

  //! Concrete implementation of createFE
  virtual MixedFiniteElement * createFE() const;

  //! Concrete implementation of create_boundary_FE
  virtual MixedFiniteElement * create_boundary_FE() const;

  //! Concrete implementation of setup_dofs
  virtual void setup_dofs();

};

#endif
