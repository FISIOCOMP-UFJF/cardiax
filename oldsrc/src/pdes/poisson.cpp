#include "poisson.hpp"
#include "../util/timer.hpp"

namespace {
  const double penalty = 1.e+10;
}

void Poisson::assemble_system()
{
  FiniteElement & fe = fespace.createFE(0);
  int n = fe.get_ndof();
  int ndofs = fespace.get_num_dofs();
  arma::mat elmat(n,n);
  arma::vec elvec(n);
  vector<int> dnums;

  ZeroFunction<double> zerofunc;

  // Please select the RHS function here
  RHSExample3D<double> rhsfunc;
  //RHSExample0<double> rhsfunc;

  // PETSC matrix and vector creation
  K.create(ndofs,ndofs,30);
  u.create(ndofs);
  f.create(ndofs);

  msg("Assembling PETSc Matrix");

  // Assemble matrix
  for(int i=0; i < mesh.get_n_elements(); i++)
  {

    fespace.get_element_dofs(i,dnums);
    calc_elmat_poisson(i, fe, elmat);
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

  // Apply Dirichlet boundary conditions u = gD on \gamma_D
  if(mesh.get_n_boundary_elements() > 0)
  {
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
        {
          int I = bdnums[k];
          int J = bdnums[l];
          K.add(I,J,belmat(k,l));
        }

      for(int k=0;k<m;k++)
        f.add(bdnums[k],belvec(k));
    }

    delete &bfe;
  }
  
  // fix boundary nodes
  if (fixed_nodes_map.size())
    FETools::apply_boundary_values(fixed_nodes_map, K, f);

  K.assemble();
  
  delete &fe;
}

void Poisson::config(std::string optfile)
{
  if( file_extension(optfile) == "xml")
  {
    read_xml_section(optfile, "poisson", "neumann", neumann_map);
    read_xml_section(optfile, "poisson", "dirichlet", dirichlet_map);
    read_xml_section(optfile, "poisson", "fix_node", fixed_nodes_map);    
  }
  else if( file_extension(optfile) == "par" )
  {
    InputFile ifile(optfile);
    ifile.read_section("neumann", neumann_map);
    ifile.read_section("dirichlet", dirichlet_map);
    ifile.read_section("fix_node", fixed_nodes_map);
    ifile.close();
  }
}

void Poisson::calc_elmat_poisson(const int iel, const FiniteElement & fe,
                                 arma::mat & elmat)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  double detJxW;
  double cond = 1.0;

  if (has_cond) cond = mesh.get_data(iel);

  Quadrature * qd = Quadrature::create(2*fe.order()-2,fe.type());
  arma::vec shape;
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);

  Mapping em = fe.get_mapping(iel);
  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    elmat += detJxW * cond * (gradn * gradn.t());
  }

  delete qd;
}

void Poisson::calc_elmat_mass(const int iel, const FiniteElement & fe,
			      arma::mat & elmat)
{
  int ndim = fe.get_ndim();
  double detJxW;
  arma::vec shape;
  arma::mat dshape;
  arma::mat jacinv(ndim,ndim);

  Mapping em = fe.get_mapping(iel);
  Quadrature * qd = Quadrature::create(2*fe.order(), fe.type());

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q), shape);
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    elmat += detJxW * (shape * shape.t());
  }
  
  delete qd;
}

void Poisson::calc_elvec_source (const int iel, const FiniteElement & fe, 
                                 const ScalarFunction<double> & rhs, 
                                 arma::vec & elvec)
{
  double detJxW, f;
  arma::vec shape;
  arma::mat dshape;

  Quadrature * qd = Quadrature::create(2*fe.order(), fe.type());
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

void Poisson::calc_robin_elmat (const int iel, const FiniteElement & fe,
                                arma::mat & elmat)
{
  double detJxW;
  arma::vec shape;
  arma::mat dshape;

  SurfaceMapping sm = fe.get_boundary_mapping(iel);
  Quadrature * qd = Quadrature::create(2*fe.order(),fe.type());

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q), shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    sm.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * sm.get_det_jacobian();

    double coeff = coeff_robin_mat(sm.get_index());
    double pdx = coeff * detJxW;
    elmat += pdx * (shape * shape.t());
  }
  delete qd;
}

void Poisson::calc_robin_elvec (const int iel, const FiniteElement & fe,
                                const ScalarFunction<double> & gD, 
                                const ScalarFunction<double> & gN, 
                                arma::vec & elvec)
{
  double detJxW;
  arma::vec shape;
  arma::mat dshape;

  SurfaceMapping sm = fe.get_boundary_mapping(iel);
  Quadrature * qd = Quadrature::create(fe.order(), fe.type());

  elvec.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q),shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    sm.calc_jacobian(dshape);

    arma::vec3 pt = qd->get_point(q);
    arma::vec3 xyz = sm.map_point(pt[0]);

    //ugD = gD(xyz);
    //ugN = gN(xyz);

    detJxW = qd->get_weight(q) * sm.get_det_jacobian();
    double coeff = coeff_robin_vec(sm.get_index());
    double pdx = coeff * detJxW;

    elvec += pdx * shape;
  }

  delete qd;
}

void Poisson::calc_mass_matrix(arma::sp_mat & mass)
{
  int n = mesh.get_nen();
  int ndofs = fespace.get_num_dofs();
  arma::mat elmat_m(n,n);
  std::vector<int> dnums;

  FiniteElement & fe = fespace.createFE(0);

  mass.set_size(ndofs,ndofs);
  
  for(int i=0; i < mesh.get_n_elements(); i++){

    calc_elmat_mass(i, fe, elmat_m);
    fespace.get_element_dofs(i,dnums);

    for(int k=0;k<n;k++)
      for(int l=0;l<n;l++)
        mass(dnums[k],dnums[l]) += elmat_m(k,l);
  }

  //M.compress();
  delete &fe;
}

double Poisson::calc_l2_error(ScalarFunction<double> & exact)
{
  FiniteElement & fe = fespace.createFE(0);
  Quadrature * qd = Quadrature::create(2*fe.order(), fe.type());
  
  int nen   = mesh.get_nen();
  int ndofs = fespace.get_num_dofs();
  double detJxW, error=0;  
  std::vector<int> dnums;
  arma::vec shape;
  arma::mat dshape;
  
  for(int iel=0; iel < mesh.get_n_elements(); iel++)
  {
    fespace.get_element_dofs(iel,dnums);

    Mapping em = fe.get_mapping(iel);
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe.calc_shape(qd->get_point(q), shape);
      fe.calc_deriv_shape(qd->get_point(q),dshape);
      em.calc_jacobian(dshape);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();

      // map gauss point to physical point and evaluate solution
      arma::vec3 pt = em.map_point(qd->get_point(q), shape);
      double ue = exact(pt);

      // interpolate FEM solution to gauss point
      double uh = 0.0;           
      for(int i=0; i<nen; i++)
        uh += solution(dnums[i]) * shape(i);      

      // sum error
      error += (ue-uh)*(ue-uh)*detJxW;
    }  
  }

  delete qd;
  delete &fe;

  return sqrt(error);  
}

double Poisson::coeff_robin_mat(int index)
{
  // Dirichlet case
  if (dirichlet_map.find(index) != dirichlet_map.end())
    return penalty;

  // Neumann case
  if (neumann_map.find(index) != neumann_map.end())
    return 0.0;

  cout << " Boundary (Dirichlet) index " << index << " not set." << endl;
  exit(0);
}

double Poisson::coeff_robin_vec(int index)
{
  std::map<int,double>::iterator it;

  // Dirichlet case
  it = dirichlet_map.find(index);
  if (it != dirichlet_map.end())
  {
    double val = (double) it->second;
    return penalty * val;
  }

  // Neumann case
  it = neumann_map.find(index);
  if (neumann_map.find(index) != neumann_map.end())
  {
    double val = (double) it->second;
    return val;
  }

  cout << " Boundary (Neumann) index " << index << " not set." << endl;
  exit(0);
}

void Poisson::init()
{
  if(!file_exists(filename))
    error("mesh file does not exist");

  std::string fext = file_extension(filename);
  if(fext == "msh")
  {
    cout << "Reading GMSH mesh" << endl;
    GmshIO gmshreader(mesh);
    gmshreader.read(filename);
  }
  else
  {
    cout << "Reading XML mesh" << endl;
    mesh.read_xml(filename);
  }

  // define the mesh for the FE space
  fespace.set_mesh(&mesh);
  fespace.config();
      
  if (mesh.get_data_size())
    has_cond = true;

  cout << mesh;
}

void Poisson::write_data(const string & filename)
{
  std::string vtuname = filename + ".vtu";

  // write output data
  writer.open(filename, 1, 1);
  writer.write_vm_step(0, solution.memptr());
  writer.close();
}

void Poisson::run(const string & name)
{
  filename = name;
  init();
  assemble_system();
  solve();
}

void Poisson::solve()
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
  ir = ls.solve(K,u,f,1.e-10);
  timer.stop();

  cout << " Number of iterations: " << ir.first << endl;
  cout << " Residual norm: " << ir.second << endl;
  cout << " Total time: " << timer.get_total_time() << endl;

  // Copy solution from PETSc Vec to my Vector
  // OLD CODE OK
  solution.resize(ndofs);
  u.get_data(solution.memptr());
  
  //
  // Post-processing step for 1D Cubic Hermite finite elements
  //
  CubicHermiteFESpace * fec = dynamic_cast<CubicHermiteFESpace*>(&fespace);
  if (fec != NULL)
  {
    cout << "FESpace detected: Cubic Hermite\n";
    cout << "Post-processing function...\n";
    arma::vec func = solution(arma::span(0,(ndofs/2)-1));
    cout << func << endl;

    cout << "Post-processing derivative...\n";
    arma::vec deriv = solution(arma::span(ndofs/2, ndofs-1));
    FETools::cubic_hermite_deriv(fespace, mesh, deriv);
    cout << deriv << endl;
  }
  
}
