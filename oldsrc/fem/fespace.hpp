#ifndef FESPACE_HPP
#define FESPACE_HPP

#include <iostream>
#include <vector>
#include "fe.hpp"
#include "mesh/mesh.hpp"

/*! 
 * Abstract base class for Finite Element Spaces
 * (actually this class is used for mapping D.O.F.s)
 *
 */
class FESpace
{
public:

  //! Default constructor
  FESpace(): mesh(), ndofs() {}
  
  //! Constructor
  FESpace(const Mesh * msh): mesh(msh) {}

  //! Destructor
  virtual ~FESpace(){}
 
  //! Number of Global Degrees of Freedom
  int get_num_dofs() const { return ndofs; }

  //! Set mesh
  void set_mesh(const Mesh * msh) { mesh = msh; }

  //! Configure DoFs
  void config();

  // *** Interface ***
  
  //! Get Element DoF numbers
  virtual int get_element_dofs (int i, std::vector<int> & dnums) const = 0;
    
  //! Get Surface DoF numbers
  virtual int get_boundary_element_dofs(int i, std::vector<int> & dnums) const = 0; 
  
  //! Create Finite Element
  virtual FiniteElement & createFE(int i) const = 0;

  //! Create Surface Finite Element
  virtual FiniteElement & create_boundary_FE(int i) const = 0;

    
protected:

  //! Mesh associated with the FE space
  const Mesh * mesh;

  //! Total number of DOFs
  int ndofs;
  
};

/*!
 *  A first simple class for H^1 finite element space of order 1.
 *  Provides access to / creation of Local Finite elements and
 *  global-local DOF-mapping.
 *
 *  H^1-conforming Finite Element Space of uniform polynomial degree 1.
*/

class H1FESpace : public FESpace
{
private:

  //! Mesh associated with the H1 FE space
  //const Mesh * mesh;

  //! Number of DOFs
  //int ndofs;
    
public:

  //! Default constructor
  H1FESpace() : FESpace() {}

  //! Constructor for linear H1-FESpace on mesh
  H1FESpace(const Mesh * msh): FESpace(msh)
  {
    ndofs = mesh->get_n_points(); //, //ndofs(mesh->get_n_points()){}
  }			       
    
  //! Destructor
  ~H1FESpace(){};

  //! Configure DoFs
  void config();
    
  //! Get Element Dof Nrs
  int get_element_dofs (int i, std::vector<int> & dnums) const;
    
  //! Get Surface DofNrs
  int get_boundary_element_dofs(int i, std::vector<int> & dnums) const; 
  
  //! Create Finite Elment
  FiniteElement & createFE(int i) const; 

  //! Create Surface Finite Element
  FiniteElement & create_boundary_FE(int i) const;

};


#endif
