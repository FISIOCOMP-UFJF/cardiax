#ifndef HERMITE_FESPACE_HPP
#define HERMITE_FESPACE_HPP

#include "fespace.hpp"

class CubicHermiteFESpace : public FESpace
{
public:

  //! Default constructor
  CubicHermiteFESpace() : FESpace() {}

  //! Constructor for cubic Hermite FESpace on mesh
  CubicHermiteFESpace(const Mesh * msh): FESpace(msh)
  {
    ndofs = 2 * mesh->get_n_points();    
  }			       
    
  //! Destructor
  ~CubicHermiteFESpace(){};

  //! Configure
  void config();
  
  //! Get element DoFs
  int get_element_dofs (int i, std::vector<int> & dnums) const;
    
  //! Get surface element DoFs
  int get_boundary_element_dofs(int i, std::vector<int> & dnums) const;

  //! Create Finite Element
  FiniteElement & createFE(int i) const; 

  //! Create Surface Finite Element
  FiniteElement & create_boundary_FE(int i) const;

protected:

  //! Global DoF number of first derivative value DoF
  int first_u2;
    
};

#endif
