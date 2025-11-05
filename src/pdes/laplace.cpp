#include "laplace.hpp"
#include "petscksp.h"

void Laplace::assemble_system()
{
  int nnode = mesh.get_nen();
  int ndofs = fespace.get_num_dofs();

  arma::mat elmat(nnode,nnode);
  arma::vec elvec(nnode);
  vector<int> dnums;

  FiniteElement & fe = fespace.createFE(0);

  ZeroFunction<double> zerofunc;

  // PETSC matrix and vector creation
  K.create(ndofs,ndofs,30);
  u.create(ndofs);
  f.create(ndofs);

  cout << "Assembling PETSc Matrix for Laplace problem\n";

  // Assemble matrix
  for(int i=0; i < mesh.get_n_elements(); i++)
  {
    calc_elmat_poisson(i, fe, elmat);
    fespace.get_element_dofs(i,dnums);

    for(int k=0;k<nnode;k++)
      for(int l=0;l<nnode;l++)
        K.add(dnums[k],dnums[l],elmat(k,l));
  }

  // Apply Dirichlet boundary conditions u = 0 on \gamma_D
  FiniteElement & bfe = fespace.create_boundary_FE(0);
  int m = bfe.get_ndof();
  arma::mat belmat(m,m);
  arma::vec belvec(m);
  vector<int> bdnums;

  // Assemble boundary bilinear form to impose boundary conditions
  for(int i=0; i < mesh.get_n_boundary_elements(); i++)
  {
    calc_robin_elmat (i, bfe, belmat);
    calc_robin_elvec (i, bfe, zerofunc, zerofunc, belvec);

    fespace.get_boundary_element_dofs(i,bdnums);

    for(int k=0; k<m; k++)
      for(int l=0; l<m; l++)
        K.add(bdnums[k], bdnums[l], belmat(k,l));

    for(int k=0; k<m; k++)
      f.add(bdnums[k],belvec(k));
  }

  // Fix boundary nodes
  if (fixed_nodes_map.size())
  {
    cout << fixed_nodes_map.size() << " constrained nodes." << endl;
    FETools::apply_boundary_values(fixed_nodes_map, K, f);
  }

  K.assemble();

  delete &fe;
  delete &bfe;
}
