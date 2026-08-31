#include "fespace.hpp"

// -------------------------------------------------------------------
//   LINEAR H1-Space
// -------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const H1FESpace & fespace)
{
  out << "H1FESpace " << std::endl;
  out << "NDofs " << fespace.get_num_dofs() << std::endl;
  return(out);
}

void H1FESpace::config()
{
  ndofs = mesh->get_n_points();
}

int H1FESpace::get_element_dofs(int i, std::vector<int> & dnums) const
{
  mesh->get_element_pt_nums(i,dnums);
  return(dnums.size());
}

int H1FESpace::get_boundary_element_dofs(int i, std::vector<int> & dnums) const
{
  mesh->get_boundary_element_pt_nums(i,dnums);
  return(dnums.size());
}

FiniteElement & H1FESpace::createFE(int i) const
{
  Element elem = mesh->get_element(i);
  ElementType etype = elem.get_type();

  if (etype == ELEM_SEGM)
  {
    FESegm * fe = new FESegm(mesh);
    return *fe;
  }
  else if (etype == ELEM_TRIG)
  {
    FETrig * fe = new FETrig(mesh);
    return *fe;
  }
  else if (etype == ELEM_QUAD)
  {
    FEQuad * fe = new FEQuad(mesh);
    return *fe;
  }
  else if (etype == ELEM_TETRA)
  {
    FETetra * fe = new FETetra(mesh);
    return *fe;
  }  
  else if (etype == ELEM_HEXA)
  {
    FEHexa * fe = new FEHexa(mesh);
    return *fe;
  }
  else
  {
    cerr << " Error: unknown element geometry in FESpace.createFE().\n";
    exit(1);
  }

}

FiniteElement & H1FESpace::create_boundary_FE(int i) const
{
  Element elem = mesh->get_boundary_element(i);
  ElementType etype = elem.get_type();

  if (etype == ELEM_SEGM)
  {
    FESegm * fe = new FESegm(mesh);
    return (fe[0]);
  }
  else if (etype == ELEM_TRIG)
  {
    FETrig * fe = new FETrig(mesh);
    return (fe[0]);
  }
  else if(etype == ELEM_QUAD)
  {
    FEQuad * fe = new FEQuad(mesh);
    return (fe[0]);
  }

  // IMPROVE THIS PIECE OF SHIT!
  std::cerr << " Error in H1FESpace::create_boundary_FE:"
	    << " boundary element not found.";
  
  FESegm * fe = NULL;
  return (fe[0]);

}

