#include "mixed_fespace.hpp"

// -------------------------------------------------------------------
// Base class MixedFESpace
// -------------------------------------------------------------------

MixedFESpace::MixedFESpace()
  : mesh_(), ndofs_()
{
  // do nothing
}

MixedFESpace::MixedFESpace(const Mesh * mesh) 
  : mesh_(mesh), ndofs_(mesh->get_n_points())
{
  // do nothing
}

// -------------------------------------------------------------------
//   LINEAR Vec H1-Space
// -------------------------------------------------------------------

MixedFiniteElement * VecH1FESpace::createFE() const
{
  // TODO: Pure meshes only i=0

  Element elem = mesh_->get_element(0);
  ElementType etype = elem.get_type();

  if (etype == ELEM_SEGM)
  {
    VecH1FESegm * fe = new VecH1FESegm(ndim_, *mesh_);
    return fe;
  }
  else if (etype == ELEM_TRIG)
  {
    VecH1FETrig * fe = new VecH1FETrig(ndim_, *mesh_);
    return fe;
  }
  else if (etype == ELEM_QUAD)
  {
    VecH1FEQuad * fe = new VecH1FEQuad(ndim_, *mesh_);
    return fe;
  }
  else if (etype == ELEM_TETRA)
  {
    VecH1FETetra * fe = new VecH1FETetra(ndim_, *mesh_);
    return fe;
  }
  else if (etype == ELEM_HEXA)
  {
    VecH1FEHexa * fe = new VecH1FEHexa(ndim_, *mesh_);
    return fe;
  }
  else
  {
    cerr << " Error: unknown element geometry in MixedFESpace.createFE().\n";
    exit(1);
    return NULL;
  }
}

MixedFiniteElement * VecH1FESpace::create_boundary_FE() const
{
  ElementType etype;

  if (mesh_->get_n_boundary_elements() > 0)
  {
    Element elem = mesh_->get_boundary_element(0);
    etype = elem.get_type();
  }
  else
  {
    return NULL;
  }

  // Create Boundary element
  if (etype == ELEM_SEGM)
  {
    VecH1FESegm * fe = new VecH1FESegm(ndim_, *mesh_);
    return fe;
  }
  else if (etype == ELEM_TRIG)
  {
    VecH1FETrig * fe = new VecH1FETrig(ndim_, *mesh_);
    return fe;
  }
  else if (etype == ELEM_QUAD)
  {
    VecH1FEQuad * fe = new VecH1FEQuad(ndim_, *mesh_);
    return fe;
  }
  else
  {
    throw std::runtime_error("Unknown element geometry in MixedFESpace.createFE()");
  }

  return NULL;
}

void VecH1FESpace::setup_dofs()
{
  // TO DO: implement!
}
/*
int VecH1FESpace::get_element_dofs_u(int i, std::vector<int> & dnums) const
{
  //
  // For the first component of u the dofs are the vertex numbers
  //     the 2nd component of u -> vertex number + number of points in the mesh
  //     the 3rd component of u -> vertex number + 2 * number of points in the mesh
  //
  vector<int> pnums;
  mesh_->get_element_pt_nums(i, pnums);

  int np_el = pnums.size();  // number of local vertices

  dnums.resize(ndim_*pnums.size());

  for(int l=0;l<np_el;l++) dnums[l] = pnums[l];

  if(ndim_ >= 2)
    for(int l=0;l<np_el;l++) dnums[l+np_el] = pnums[l]+first_u2_;

  if(ndim_ >= 3)
    for(int l=0;l<np_el;l++) dnums[l+2*np_el] = pnums[l]+first_u3_;

  return(dnums.size());
}

int VecH1FESpace::get_boundary_element_dofs_u(int i, std::vector<int> & dnums) const
{
  //
  // For the first component of u the dofs are the vertex numbers
  //     the 2nd component of u -> vertex number + number of points in the mesh
  //     the 3rd component of u -> vertex number + 2*number of points in the mesh
  //
  vector<int> pnums;
  mesh_->get_boundary_element_pt_nums(i,pnums);

  // number of element nodes
  int neln = pnums.size();  

  dnums.resize(ndim_*pnums.size());

  for(int l=0; l<neln; l++) dnums[l] = pnums[l];

  if(ndim_ >= 2)
    for(int l=0; l<neln; l++) dnums[l+neln] = pnums[l]+first_u2_;

  if(ndim_ >= 3)
    for(int l=0; l<neln; l++) dnums[l+2*neln] = pnums[l]+first_u3_;
  
    return(dnums.size());
}
*/

int VecH1FESpace::get_element_dofs_u(int i, std::vector<int> & dnums) const
{
  //
  // For the first component of u the dofs are the vertex numbers
  //     the 2nd component of u -> vertex number + number of points in the mesh
  //     the 3rd component of u -> vertex number + 2 * number of points in the mesh
  //
  vector<int> pnums;
  mesh_->get_element_pt_nums(i, pnums);

  int np_el = pnums.size();  // number of local vertices

  dnums.resize(ndim_*pnums.size());

  for(int l=0;l<np_el;l++) dnums[l] = ndim_*pnums[l];

  if(ndim_ >= 2)
    for(int l=0;l<np_el;l++) dnums[l+np_el] = ndim_*pnums[l]+1;

  if(ndim_ >= 3)
    for(int l=0;l<np_el;l++) dnums[l+2*np_el] = ndim_*pnums[l]+2;


  return(dnums.size());
}


int VecH1FESpace::get_boundary_element_dofs_u(int i, std::vector<int> & dnums) const
{
  //
  // For the first component of u the dofs are the vertex numbers
  //     the 2nd component of u -> vertex number + number of points in the mesh
  //     the 3rd component of u -> vertex number + 2*number of points in the mesh
  //
  vector<int> pnums;
  mesh_->get_boundary_element_pt_nums(i,pnums);

  // number of element nodes
  int neln = pnums.size();

  dnums.resize(ndim_*pnums.size());

  for(int l=0; l<neln; l++) dnums[l] = ndim_*pnums[l];

  if(ndim_ >= 2)
    for(int l=0; l<neln; l++) dnums[l+neln] = ndim_*pnums[l]+1;

  if(ndim_ >= 3)
    for(int l=0; l<neln; l++) dnums[l+2*neln] = ndim_*pnums[l]+2;

/*
    cout<<"\n";
    for(int i=0; i<pnums.size(); i++)
        cout<<pnums[i]<<" ";

    cout<<"\n";

    for(int i=0; i<dnums.size(); i++)
        cout<<dnums[i]<<" ";
    exit(0);
*/

  return(dnums.size());
}

