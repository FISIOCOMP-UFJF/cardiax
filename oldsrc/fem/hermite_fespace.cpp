#include "hermite_fespace.hpp"

void CubicHermiteFESpace::config()
{
  if(mesh)
  {
    ndofs = 2 * mesh->get_n_points();
    first_u2 = mesh->get_n_points();
  }
}

int CubicHermiteFESpace::get_element_dofs(int i, std::vector<int> & dnums) const
{
  // TODO: extend this for 2D and 3D elements
  
  // So far it works for 1D elements only
  // First component of u the DoFs are the function values
  // Second component of u the DoFs are the derivative values

  // Get global point/vertices/node numbers
  std::vector<int> pnums;
  mesh->get_element_pt_nums(i, pnums);

  // number of points/vertices of an element
  int npel = pnums.size();
  
  // resize dof map vector
  dnums.resize(2*npel);

  for(int k=0; k<npel; k++)
    dnums[k] = pnums[k];

  for(int k=0; k<npel; k++)
    dnums[k+npel] = pnums[k] + first_u2;

  return(dnums.size());
  
}

int CubicHermiteFESpace::get_boundary_element_dofs(int i, std::vector<int> & dnums) const
{
  // implementar
  return 0;
}

FiniteElement & CubicHermiteFESpace::createFE(int i) const
{
  Element elem = mesh->get_element(i);
  ElementType etype = elem.get_type();

  if(etype == ELEM_SEGM)
  {
    FECubicHermite * fe = new FECubicHermite(mesh);
    return *fe;
  }
  else
  {
    throw std::runtime_error("CubicHermite FE Space: unknown element type");
  }

  FECubicHermite * fe = NULL;
  return (fe[0]); 
}


FiniteElement & CubicHermiteFESpace::create_boundary_FE(int i) const
{
  // TODO: need to implement
  FECubicHermite * fe = NULL;
  return (fe[0]); 
}
