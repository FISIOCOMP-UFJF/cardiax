//
// Created by rocha on 1/11/16.
//

#include "projection.hpp"
#include "../util/timer.hpp"

void Projection::assemble_system()
{
  FiniteElement & fe = fespace.createFE(0);

  int n = fe.get_ndof(); //int n = msh.get_nen();
  int ndofs = fespace.get_num_dofs();

  arma::mat elmat(n,n);
  arma::vec elvec(n);
  std::vector<int> dnums;

  // Function to compute L2 projection
  RHSExample1<double> rhsfunc;

  // PETSC matrix and vector creation
  K.create(ndofs,ndofs,30);
  u.create(ndofs);
  f.create(ndofs);

  msg("Assembling PETSc Matrix");

  // Assemble matrix
  for(int i=0; i<msh.get_n_elements(); i++)
  {
    std::cout << "Element " << i << std::endl;

    calc_elmat_mass(i, fe, elmat);
    calc_elvec_source (i, fe, rhsfunc, elvec);

    fespace.get_element_dofs(i,dnums);

    for(int k=0; k<n; k++)
      for(int l=0;l<n;l++)
      {
        int I = dnums[k];
        int J = dnums[l];
        K.add(I, J, elmat(k,l));
      }

    for(int k=0; k<n; k++)
      f.add(dnums[k], elvec(k));
  }

  K.assemble();
  f.assemble();

  delete &fe;
}

void Projection::config(std::string optfile)
{
  //
  // TODO: need to start using the XML reader here
  //
  /*
  InputFile ifile(optfile);
  ifile.read_section("neumann", neumann_map);
  ifile.read_section("dirichlet", dirichlet_map);
  ifile.read_section("fix_node", fixed_nodes_map);
  ifile.close();
  */
}

void Projection::calc_elmat_mass(const int iel, const FiniteElement & fe,
                                arma::mat & elmat)
{
  double detJxW;
  arma::vec shape;
  arma::mat dshape;

  Mapping em = fe.get_mapping(iel);
  Quadrature * qd = Quadrature::create(2*fe.order(), fe.type());

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q), shape);
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    for(uint i=0; i<shape.size(); i++)
      for(uint j=0; j<shape.size(); j++)
        elmat(i,j) += detJxW * (shape(i) * shape(j));
  }

  delete qd;
}

void Projection::calc_elmat_mixed_mass(const int ielF,
                                       const int ielC,
                                       const FiniteElement & feF,
                                       const FiniteElement & feC,
                                       arma::mat & elmat)
{
  double detJxW;
  arma::vec shapeF, shapeC;
  arma::mat dshapeF, dshapeC;

  // Fine
  Mapping emF = feF.get_mapping(ielF);
  Quadrature * qd = Quadrature::create(2*feF.order(), feF.type());

  // Coarse
  Mapping emC = feC.get_mapping(ielC);

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    feF.calc_shape(qd->get_point(q), shapeF);
    feF.calc_deriv_shape(qd->get_point(q),dshapeF);
    emF.calc_jacobian(dshapeF);

    feC.calc_shape(qd->get_point(q), shapeC);
    feC.calc_deriv_shape(qd->get_point(q),dshapeC);

    detJxW = qd->get_weight(q) * emF.get_det_jacobian();

    for(uint i=0; i<shapeF.size(); i++)
      for (uint j = 0; j < shapeC.size(); j++)
        elmat(i, j) += detJxW * ( shapeF(i) * shapeC(j) );
  }

  delete qd;
}

void Projection::calc_elvec_source (const int iel, const FiniteElement & fe,
                                 const ScalarFunction<double> & rhs,
                                 arma::vec & elvec)
{
  double detJxW, f;
  arma::vec shape;
  arma::mat dshape;

  //Quadrature * qd = Quadrature::create(7fe.order(), fe.type());
  Quadrature * qd = Quadrature::create(2, fe.type());
  Mapping em = fe.get_mapping(iel);

  elvec.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q),shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    em.calc_jacobian(dshape);

    f = rhs( em.map_point(qd->get_point(q), shape) );
    detJxW = qd->get_weight(q) * em.get_det_jacobian();

    elvec += (f * detJxW) * shape;
  }

  delete qd;
}

void Projection::init()
{
  if(!file_exists(filename)) error("mesh file does not exist");

  std::string fext = file_extension(filename);
  if(fext == "msh")
  {
    GmshIO gmshreader(msh);
    gmshreader.read(filename);
  }
  else if(fext == "xml")
    msh.read_xml(filename);

  // define the mesh for the FE space
  fespace.set_mesh(&msh);
  fespace.config();

  std::cout << msh;
}

void Projection::write_data(const string & filename, bool as_text)
{
  std::string vtuname = filename + ".vtu";

  if(!as_text)
  {
    // write output data
    writer.open(filename, 1, 1);
    writer.write_vm_step(0, solution.memptr());
    writer.close();
  }
  //else
  //vtkout.write(solution, vtuname, false);
}

void Projection::run(const string & name)
{
  filename = name;
  init();
  assemble_system();
  solve();
}

void Projection::solve()
{
  Timer timer;
  int ndofs = fespace.get_num_dofs();
  std::pair<PetscInt,PetscReal> ir;

  // prepare linear system solver
  petsc::LinearSolver ls;
  ls.init();

  // solve
  msg("Solving linear system");
  timer.start();
  ir = ls.solve(K,u,f,1.e-8);
  timer.stop();

  cout << " Number of iterations: " << ir.first << endl;
  cout << " Residual norm: " << ir.second << endl;
  cout << " Total time: " << timer.get_total_time() << endl;

  // Copy solution from PETSc Vec to my Vector
  solution.resize(ndofs);
  u.get_data(solution.memptr());
}


void Projection::transfer(Mesh & target_mesh, arma::vec & target_solution)
{
  cout << "\nL2 Projection to transfer the solution"
       << "from a Source mesh to a Target mesh" << endl;

  GmshIO gmshreader(target_mesh);
  gmshreader.read("/data/meshes/test/unit5.msh");

  petsc::Matrix P;

  // define the mesh for the FE space
  H1FESpace fespaceF;
  fespaceF.set_mesh(&target_mesh);
  fespaceF.config();

  FiniteElement & feC = fespace.createFE(0);
  FiniteElement & feF = fespaceF.createFE(0);

  int n = feF.get_ndof();
  int ndofsC = fespaceF.get_num_dofs();
  int ndofsF = fespaceF.get_num_dofs();

  arma::mat elmat(n,n);
  std::vector<int> dnumsF, dnumsC;

  P.create(ndofsF,ndofsC,30);

  for(int iF=0; iF<msh.get_n_elements(); iF++)
  {
    // GAMBIARRAAAAA-
    int iC = 0;
    /// -------------

    calc_elmat_mixed_mass(iF, iC, feF, feC, elmat);

    fespaceF.get_element_dofs(iF,dnumsF);
    fespace.get_element_dofs(iC,dnumsC);

    for(int k=0; k<n; k++)
      for(int l=0;l<n;l++)
      {
        int I = dnumsF[k];
        int J = dnumsF[l];
        P.add(I, J, elmat(k,l));
      }

  }

  P.assemble();
}