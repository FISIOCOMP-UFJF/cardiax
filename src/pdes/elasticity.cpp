#include "elasticity.hpp"
#include "util/pugixml.hpp"

#define DOBRAR

Elasticity::Elasticity() : writer(&msh)
{
  // setup default parameter values
  parameters.add("young", 0.5);
  parameters.add("poisson", 0.14);
  parameters.add("num_components", 0);
  parameters.add("penalty", 1.e+10);
}

Elasticity::~Elasticity()
{
  // do nothing
}

void Elasticity::run(const string & name)
{
  int nc = parameters["num_components"];
  if(nc == 0) print_error("elasticity.cpp", "run()", "num_components==0");

  filename = name;

  init();
  assemble_system();
  solve();
  cout<<"Solved\n";
}

void Elasticity::init()
{
  std::string extension = file_extension(filename);

  if(extension == "msh")
  {
    cout << "Reading GMSH mesh" << endl;
    GmshIO gmshreader(msh);
    gmshreader.read(filename);
  }
  else if(extension == "xml")
  {
    cout << "Reading XML mesh" << endl;
    msh.read_xml(filename);
  }

  fespace.set_mesh(&msh);
  cout << msh;
}

void Elasticity::assemble_system()
{
#ifdef DOBRAR
  petsc::Matrix K_l;
  petsc::Vector u_l,f_l;

  int n = msh.get_n_dim() * msh.get_nen();
  int ndofs = fespace.get_ndofs_u();

  arma::mat elmat(n,n);
  arma::vec elvec(n);
  vector<int> dnums;
  Timer asm_time;
  asm_time.start();

  // PETSC matrix and vector creation
  K.create(2*ndofs,2*ndofs,240);
  u.create(2*ndofs);
  f.create(2*ndofs);

  K_l.create(ndofs,ndofs,120);
  u_l.create(ndofs);
  f_l.create(ndofs);

  cout << "Number of DoFs: " << 2*ndofs << endl;
  cout << "Assembling PETSc Matrix" << endl;

  //
  // Assemble bilinear form - loop over elements
  //
  MixedFiniteElement * fe = fespace.createFE();
  for(int i=0; i<msh.get_n_elements(); i++)
  {
    fespace.get_element_dofs_u(i,dnums);
    calc_elmat_stiffness(i, fe, elmat);
    calc_elvec_source(i, fe, elvec);

    for(int k=0;k<n;k++)
      for(int l=0;l<n;l++)
        K_l.add(dnums[k],dnums[l],elmat(k,l));

    for(int k=0;k<n;k++)
      f_l.add(dnums[k],elvec(k));
  }

  //
  // Assemble boundary bilinearforms to impose boundary conditions
  // Loop over surface elements
  //
  cout << "Assembling boundary matrices" << endl;
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe != NULL)
  {
    int m = bfe->get_ndofs_u();
    arma::mat belmat(m,m);
    arma::vec belvec(m);
    vector<int> bdnums;

    for(int i=0; i<msh.get_n_boundary_elements(); i++)
    {
      calc_robin_elmat(i, bfe, belmat);
      calc_robin_elvec(i, bfe, belvec);

      fespace.get_boundary_element_dofs_u(i,bdnums);

      for(int k=0;k<m;k++)
        for(int l=0;l<m;l++)
          K_l.add(bdnums[k],bdnums[l],belmat(k,l));

      for(int k=0;k<m;k++)
        f_l.add(bdnums[k],belvec(k));
    }
  }

  // Fix boundary nodes
  if (fixed_nodes_map.size())
  {
    cout << "Fixing boundary nodes \n";

    int nod, idx, dir;
    double val;
    NodalData inf;

    cout << " Number of prescribed displacements ";
    cout << fixed_nodes_map.size() << endl;

    std::map<int,double> boundary_values;
    std::map<int,NodalData>::iterator it;
    for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
    {
      nod = it->first;
      inf = it->second;
      dir = inf.first;
      val = inf.second;
      idx = (dir * msh.get_n_points()) + nod;
      boundary_values.insert( std::pair<int,double>(idx,val) );
    }

    FETools::apply_boundary_values(boundary_values, K_l, f_l);
  }

  //int ln = K_l.size();
  //int nnz  = K_l.get_nnz();

  // define matrix
  //int *ia = new int[ln+1];
  //int *ja = new int[nnz];
  //double *v = new double[nnz];

  //K_l.get_CSR(&ln, ia, ja, v);


  int ln, lnnz;
  PetscBool done;
  const int *pia, *pja;
  double * vals, *vf;
  vf = f_l.get_array();

  PetscErrorCode ierr = MatGetRowIJ(K_l.mat(), 0, PETSC_FALSE, PETSC_FALSE, &ln, &pia, &pja, &done);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJGetArray(K_l.mat(), &vals);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);



  for(int i=0;i<ln;i++)
  {
    for(int j=pia[i]; j<pia[i+1]; j++)
    {
      //std::cout<<i<<" "<<pja[j] <<" "<<vals[j] <<std::endl;
      K.add(i, pja[j], vals[j]);

      K.add(i+ln, pja[j]+ln, vals[j]);
    }

    f.add(i, vf[i]);
    f.add(i+ln, vf[i]);

  }

  ierr = MatRestoreRowIJ(K_l.mat(), 0, PETSC_FALSE, PETSC_FALSE, &ln, &pia, &pja, &done);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = MatSeqAIJRestoreArray(K_l.mat(), &vals);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if(!done) throw std::runtime_error("error in MatGetRowIJ");



  K.assemble();
  f.assemble();

  delete fe;
  delete bfe;

  asm_time.stop();
  cout << "Total matrix assembly time: " << asm_time.get_total_time() << endl;

#endif

#ifndef DOBRAR

  int n = msh.get_n_dim() * msh.get_nen();
  int ndofs = fespace.get_ndofs_u();

  arma::mat elmat(n,n);
  arma::vec elvec(n);
  vector<int> dnums;
  Timer asm_time;
  asm_time.start();

  // PETSC matrix and vector creation
  K.create(ndofs,ndofs,120);
  u.create(ndofs);
  f.create(ndofs);

  cout << "Number of DoFs: " << ndofs << endl;
  cout << "Assembling PETSc Matrix" << endl;

  //
  // Assemble bilinear form - loop over elements
  //
  MixedFiniteElement * fe = fespace.createFE();
  for(int i=0; i<msh.get_n_elements(); i++)
  {
    fespace.get_element_dofs_u(i,dnums);
    calc_elmat_stiffness(i, fe, elmat);
    calc_elvec_source(i, fe, elvec);

    for(int k=0;k<n;k++)
      for(int l=0;l<n;l++)
        K.add(dnums[k],dnums[l],elmat(k,l));

    for(int k=0;k<n;k++)
      f.add(dnums[k],elvec(k));
  }

  //
  // Assemble boundary bilinearforms to impose boundary conditions
  // Loop over surface elements
  //
  cout << "Assembling boundary matrices" << endl;
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe != NULL)
  {
    int m = bfe->get_ndofs_u();
    arma::mat belmat(m,m);
    arma::vec belvec(m);
    vector<int> bdnums;

    for(int i=0; i<msh.get_n_boundary_elements(); i++)
    {
      calc_robin_elmat(i, bfe, belmat);
      calc_robin_elvec(i, bfe, belvec);

      fespace.get_boundary_element_dofs_u(i,bdnums);

      for(int k=0;k<m;k++)
        for(int l=0;l<m;l++)
          K.add(bdnums[k],bdnums[l],belmat(k,l));

      for(int k=0;k<m;k++)
        f.add(bdnums[k],belvec(k));
    }
  }

  // Fix boundary nodes
  if (fixed_nodes_map.size())
  {
    cout << "Fixing boundary nodes \n";

    int nod, idx, dir;
    double val;
    NodalData inf;

    cout << " Number of prescribed displacements ";
    cout << fixed_nodes_map.size() << endl;

    std::map<int,double> boundary_values;
    std::map<int,NodalData>::iterator it;
    for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
    {
      nod = it->first;
      inf = it->second;
      dir = inf.first;
      val = inf.second;
      idx = (dir * msh.get_n_points()) + nod;
      boundary_values.insert( std::pair<int,double>(idx,val) );
    }

    FETools::apply_boundary_values(boundary_values, K, f);
  }


  K.assemble();
  f.assemble();

  delete fe;
  delete bfe;

  asm_time.stop();
  cout << "Total matrix assembly time: " << asm_time.get_total_time() << endl;

#endif


}

inline void Elasticity::calc_B_matrix (const arma::mat & gradn, arma::mat & B)
{
  B.zeros();

  int Nr = gradn.n_rows;   // ndof
  int Nd = msh.get_n_dim(); // ndim
  int Nb = Nr/Nd;          // nubf

  if (Nd == 2)
  {
    // e,x
    for(int i=0; i<Nb; i++)
      B(0,i) = gradn(i,0);

    // e,y
    for(int i=Nb; i<2*Nb; i++)
      B(1,i) = gradn(i,1);

    // ex,y
    for(int i=0; i<Nb; i++)
      B(2,i) = gradn(i,1);
    for(int i=Nb; i<2*Nb; i++)
      B(2,i) = gradn(i,0);
  }
  else if (Nd == 3)
  {
    // e,x
    for(int i=0; i<Nb; i++)
      B(0,i) = gradn(i,0);

    // e,y
    for(int i=Nb; i<2*Nb; i++)
      B(1,i) = gradn(i,1);

    // e,z
    for(int i=2*Nb; i<3*Nb; i++)
      B(2,i) = gradn(i,2);

    // ex,y
    for(int i=0; i<Nb; i++)
      B(3,i) = gradn(i,1);
    for(int i=Nb; i<2*Nb; i++)
      B(3,i) = gradn(i,0);

    // ey,z
    for(int i=Nb; i<2*Nb; i++)
      B(4,i) = gradn(i,2);
    for(int i=2*Nb; i<3*Nb; i++)
      B(4,i) = gradn(i,1);

    // ex,z
    for(int i=0; i<Nb; i++)
      B(5,i) = gradn(i,2);
    for(int i=2*Nb; i<3*Nb; i++)
      B(5,i) = gradn(i,0);
  }
}

inline void Elasticity::calc_D_matrix (const arma::vec3 & X, int index, arma::mat & D)
{
  double e = parameters["young"];
  double v = parameters["poisson"];

  D.zeros();

  if (elastype == PLANE_STRAIN)
  {
    const double c = e*(1.0-v)/((1.0+v)*(1.0-2.0*v));
    D(0,0) = c;
    D(0,1) = c*(v/(1.0-v));
    D(1,0) = c*(v/(1.0-v));
    D(1,1) = c;
    D(2,2) = c*((1.0-2.0*v)/(2.0*(1.0-v)));
  }
  else if (elastype == PLANE_STRESS)
  {
    const double c = e/(1.0-v*v);
    D(0,0) = c;
    D(0,1) = c*v;
    D(1,0) = c*v;
    D(1,1) = c;
    D(2,2) = c*(0.5*(1.0-v));
  }
  else if (elastype == THREE_DIM)
  {
    const double c = e*(1.0-v)/((1.0+v)*(1.0-2.0*v));
    D(0,0) = c;
    D(0,1) = c*(v/(1.0-v));
    D(0,2) = c*(v/(1.0-v));
    D(1,1) = c;
    D(1,2) = c*(v/(1.0-v));
    D(2,2) = c;
    D(3,3) = (1.0-2.0*v)/(2.0*(1.0-v));
    D(4,4) = (1.0-2.0*v)/(2.0*(1.0-v));
    D(5,5) = (1.0-2.0*v)/(2.0*(1.0-v));

    // symmetric
    D(1,0) = D(0,1);
    D(2,0) = D(0,2);
    D(2,1) = D(1,2);
  }
}

inline void Elasticity::calc_H_matrix (const arma::vec & shape, arma::mat & H)
{
  H.zeros();

  int Nr = shape.n_rows;    // ndof
  int Nd = msh.get_n_dim(); // ndim
  int Nb = Nr/Nd;           // nubf

  if (Nd == 1) 
  {
    H(0,0) = shape(0);
    H(0,1) = shape(1);
    H(1,2) = shape(2);
    H(1,3) = shape(3);
  }
  else if (Nd == 2 && Nb == 2) // line
  {
    H(0,0) = shape(0);
    H(0,1) = shape(1);
    H(1,2) = shape(2);
    H(1,3) = shape(3);
  }
  else if (Nd == 3 && Nb == 4) // quad
  {
    for(int i=0; i<Nb; i++)
      H(0,i) = shape(i%Nb);

    for(int i=Nb; i<2*Nb; i++)
      H(1,i) = shape(i%Nb);

    for(int i=2*Nb; i<3*Nb; i++)
      H(2,i) = shape(i%Nb);
  }
}

inline void Elasticity::calc_C_matrix (const arma::vec3 & a, arma::mat & C)
{
  //C.zeros();
  C(0,0) = a(0);
  C(1,1) = a(1);
  if(msh.get_n_dim() == 3)
    C(2,2) = a(2);
}

void Elasticity::solve()
{
  std::pair<PetscInt,PetscReal> ir;

  petsc::LinearSolver ls;
  ls.init();

  cout << "Solving linear system using PETSc KSP" << endl;
#ifndef DOBRAR
  ir = ls.solve(K,u,f,1.e-8);
#endif
#ifdef DOBRAR
  ir = ls.solveFieldSplit(K,u,f,1.e-8);
#endif
  cout << " Number of iterations: " << ir.first << endl;
  cout << " Residual norm: " << ir.second << endl;


  // Copy solution from PETSc Vec to my Vector
  solution.resize(fespace.get_ndofs());

  u.get_data(solution.memptr());

}

void Elasticity::calc_elmat_stiffness (const int eindex,
                                       const MixedFiniteElement * fe,
                                       arma::mat & elmat)
{
  int index;
  int numc = parameters["num_components"];
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  double detJxW;
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat B(numc,ndof);
  arma::mat D(numc,numc);

  // Create quadrature rule (need to check if order is ok)
  Quadrature * qd;
  qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type()); 
  
  Mapping em = fe->get_mapping(eindex);

  elmat.zeros();
  index = em.get_index();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe->calc_deriv_shape_u(qd->get_point(q),dshape);
    em.calc_jacobian(dshape, nubf);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();

    gradn = (dshape * jacinv);

    // compute B and D matrices
    arma::vec3 pt = em.map_point(qd->get_point(q));
    calc_B_matrix (gradn, B);
    calc_D_matrix (pt, index, D);

    // 2 mu eps(u) : eps(v) + lambda div(u) div(v)
    elmat += (B.t() * D * B) * detJxW;
  }

  delete qd;
}

void Elasticity::calc_elvec_source (const int eindex,
                                    const MixedFiniteElement * fe,
                                    arma::vec & elvec)
{
  int index;
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  double detJxW;
  arma::vec shape;
  arma::mat dshape;

  // Create quadrature rule (need to check if order is ok)
  Quadrature * qd = Quadrature::create(fe->get_order_u(), fe->get_type());
  Mapping em = fe->get_mapping(eindex);

  elvec.zeros();
  index = em.get_index();
  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe->calc_shape_u(qd->get_point(q),shape);
    fe->calc_deriv_shape_u(qd->get_point(q),dshape);
    em.calc_jacobian(dshape, nubf);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();

    arma::vec3 pt = em.map_point(qd->get_point(q));
    arma::vec3 f(arma::fill::zeros); // body forces - load volume

    for(int i=0; i<ndim; i++)
      for(int j=0; j<nubf; j++)
        elvec(j+i*nubf) += f(i) * shape(j+i*nubf) * detJxW;
  }

  delete qd;
}

void Elasticity::calc_robin_elmat (const int eindex,
                                   const MixedFiniteElement * fe,
                                   arma::mat & elmat)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  double detJxW;
  arma::mat H(ndim,ndof);
  arma::mat C(ndim,ndim);
  arma::vec shape;
  arma::mat dshape;

  // Create quadrature rule (need to check if order is ok)
  Quadrature * qd;
  qd = Quadrature::create(2*fe->get_order_u(), fe->get_type());

  SurfaceMapping sm = fe->get_boundary_mapping(eindex);

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe->calc_shape_u(qd->get_point(q),shape);
    fe->calc_deriv_shape_u(qd->get_point(q),dshape);

    sm.calc_jacobian(dshape, nubf);

    detJxW = qd->get_weight(q) * sm.get_det_jacobian();

    //arma::vec3 pt = sm.map_point(qd->get_point(q)[0]);
    arma::vec3 a = coeff_robin_mat(sm.get_index());
    
    calc_C_matrix (a, C);
    calc_H_matrix (shape, H);

    elmat += H.t() * C * H * detJxW;
  }

  delete qd;
}

void Elasticity::calc_robin_elvec (const int eindex,
                                   const MixedFiniteElement * fe,
                                   arma::vec & elvec)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;

  double detJxW;
  arma::vec shape;
  arma::mat dshape;

  // Create quadrature rule (need to check if order is ok)
  Quadrature * qd;
  qd = Quadrature::create(fe->get_order_u(), fe->get_type());

  SurfaceMapping sm = fe->get_boundary_mapping(eindex);
  elvec.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe->calc_shape_u(qd->get_point(q),shape);
    fe->calc_deriv_shape_u(qd->get_point(q),dshape);
    sm.calc_jacobian(dshape, nubf);

    detJxW = qd->get_weight(q) * sm.get_det_jacobian();

    arma::vec3 pt = sm.map_point(qd->get_point(q)[0]);
    arma::vec3 f = coeff_robin_vec(sm.get_index());

    /*
    elvec(0) += f(0) * shape(0) * detJxW;
    elvec(1) += f(0) * shape(1) * detJxW;
    elvec(2) += f(1) * shape(2) * detJxW;
    elvec(3) += f(1) * shape(3) * detJxW;
    */

    for(int i=0; i<ndim; i++)
      for(int j=0; j<nubf; j++)
        elvec(j+i*nubf) += f(i) * shape(j+i*nubf) * detJxW;
  }

  delete qd;
}

arma::vec3 Elasticity::coeff_robin_mat(int index)
{
  arma::vec3 alpha(arma::fill::zeros);

  // Dirichlet case
  std::multimap<int,int>::iterator it;
  for (it=dirichlet_map.begin(); it != dirichlet_map.end(); ++it)
  {
    if (it->first == index)
    {
      double penalty = parameters["penalty"];
      alpha(0) = penalty;
      alpha(1) = penalty;
      alpha(2) = penalty;
    }
  }
  return alpha;
}

arma::vec3 Elasticity::coeff_robin_vec(int index)
{
  arma::vec3 t(arma::fill::zeros);
  std::map<int,arma::vec3>::iterator it;

  // Neumann case
  it = neumann_map.find(index);
  if (neumann_map.find(index) != neumann_map.end())
  {
    t = it->second;
    return t;
  }

  return t;
}

void Elasticity::config(std::string & parfile)
{
  cout << "Reading parameters file" << endl;

  std::string etype;
  std::string extension = file_extension(parfile);

  if(extension == "xml")
  {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(parfile.c_str());

    if (!result)
    {
      std::cout << "XML parsed with errors, attr value: ["
      << doc.child("mesh").attribute("attr").value() << "]\n";
      std::cout << "Error description: " << result.description() << "\n";
      std::cout << "Error offset: " << result.offset;
      std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
    }

    // problem type
    pugi::xml_attribute type = doc.child("elasticity").attribute("type");
    if (std::strncmp(type.as_string(), "PLANE_STRAIN", 12) == 0)
      elastype = PLANE_STRAIN;
    else if (std::strncmp(type.as_string(), "PLANE_STRESS", 12) == 0)
      elastype = PLANE_STRESS;
    else if (std::strncmp(type.as_string(), "THREE_DIM", 9) == 0)
      elastype = THREE_DIM;

    set_type(elastype);

    //
    // reading boundary conditions
    //

    // neumann
    pugi::xml_node nodes = doc.child("elasticity").child("neumann");
    for(pugi::xml_node node = nodes.child("node") ; node ;
        node = node.next_sibling("node"))
    {
      //int idx = node.attribute("id").as_int();
      int marker = node.attribute("marker").as_int();
      arma::vec3 t;
      t(0) = node.attribute("t0").as_double();
      t(1) = node.attribute("t1").as_double();
      t(2) = node.attribute("t2").as_double();
      neumann_map.insert( std::pair<int,arma::vec3>(marker,t) );
    }

    // dirichlet
    pugi::xml_node dnodes = doc.child("elasticity").child("dirichlet");
    for(pugi::xml_node node = dnodes.child("node") ; node ;
        node = node.next_sibling("node"))
    {
      //int idx = node.attribute("id").as_int();
      int marker = node.attribute("marker").as_int();
      int dir = node.attribute("direction").as_int();
      dirichlet_map.insert( std::pair<int,int>(marker,dir) );
    }

    // dirichlet - fixed nodes
    pugi::xml_node fnodes = doc.child("elasticity").child("prescribed_displacement");
    for(pugi::xml_node node = fnodes.child("node") ; node ;
        node = node.next_sibling("node"))
    {
      int idx = node.attribute("id").as_int();
      int dir = node.attribute("direction").as_int();
      double val = node.attribute("value").as_double();
      NodalData nd(dir,val);
      fixed_nodes_map.insert( std::pair<int,NodalData>(idx,nd) );
    }

    // TODO: implementar Young , PoissonRation etc...
    //material coef -> via parameters
    // dirichlet - fixed nodes
    pugi::xml_node params = doc.child("elasticity").child("parameters");
    if(params)
    {
      double young = params.child("young").text().as_double();
      double poisson = params.child("poisson").text().as_double();
      parameters["young"] = young;
      parameters["poisson"] = poisson;
      cout << " Young modulus: " << young << endl;
      cout << " Poisson ratio: " << poisson << endl;
    }
  }

  if(extension == "par")
  {
    InputFile ifile(parfile);
    ifile.read("problem_type", etype);

    if (std::strncmp(etype.c_str(), "PLANE_STRAIN", 12) == 0)
      elastype = PLANE_STRAIN;
    else if (std::strncmp(etype.c_str(), "PLANE_STRESS", 12) == 0)
      elastype = PLANE_STRESS;
    else if (std::strncmp(etype.c_str(), "THREE_DIM", 9) == 0)
      elastype = THREE_DIM;

    set_type(elastype);
    ifile.read_section("dirichlet", dirichlet_map);
    ifile.read_section("neumann", neumann_map);
    ifile.read_section("prescribed_displacement", fixed_nodes_map);
    ifile.close();
  }
}

void Elasticity::calc_element_stress(Mapping &em, const MixedFiniteElement *fe,
                                     const arma::vec &u, const arma::vec3 &X,
                                     arma::vec &sigma)
{
  int index;
  int numc = parameters["num_components"];
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;

  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat B(numc,ndof), D(numc,numc);

  sigma.resize(numc);

  index = em.get_index();

  fe->calc_deriv_shape_u(X,dshape);
  em.calc_jacobian(dshape, nubf);
  jacinv = em.get_inv_jacobian();

  gradn = dshape * jacinv;
  
  calc_B_matrix (gradn, B);
  calc_D_matrix (X, index, D);

  sigma = D * B * u;
}

void Elasticity::get_elvec(const std::vector<int> & vdnums,
                           const arma::vec & u,
                           arma::vec & elu)
{
  int nvd = vdnums.size();
  elu.resize(nvd);
  for(int i=0; i<nvd; i++)
    elu(i) = u(vdnums[i]);
}

void Elasticity::set_type(ElasticityType type)
{
  elastype = type;
  switch (elastype)
  {
    case THREE_DIM:
      parameters["num_components"] = 6;
      break;
    case PLANE_STRAIN:
      parameters["num_components"] = 3;
      break;
    case PLANE_STRESS:
      parameters["num_components"] = 3;
      break;
    default:
      std::cerr << " Error in elasticity.hpp: unknown ElasticityType.\n";
      exit(1);
  }
  cout << "Elasticity problem: " << elastype << endl;
}

void Elasticity::von_mises(const arma::vec & u, arma::vec & sigma)
{
  int nels = msh.get_n_elements();
  vector<int> dnums,sdnums;

  // Scalar field associated with Von Mises Stress in each element
  sigma.resize(nels);
  sigma.zeros();

  for(int i=0;i<nels; i++)
  {
    MixedFiniteElement * fe = fespace.createFE();
    Mapping em = msh.get_element_mapping(i);
    arma::vec elvecu, elsigma;

    fespace.get_element_dofs_u(i,dnums);

    get_elvec(dnums, solution, elvecu);
    calc_element_stress(em, fe, elvecu, msh.get_element_centroid(i), elsigma);

    double v = 0.3;
    double sigma11 = elsigma(0);
    double sigma22 = elsigma(1);
    double sigma12 = elsigma(2);
    double sigma33 = v*(sigma11 + sigma22);

    double s1 = sigma11*sigma11;
    double s2 = sigma22*sigma22;
    double s3 = sigma33*sigma33;
    double s4 = -(sigma11*sigma22);
    double s5 = -(sigma22*sigma33);
    double s6 = -(sigma33*sigma11);
    double s7 = (3.0*sigma12*sigma12);

    sigma(i) = sqrt(s1 + s2 + s3 + s4 + s5 + s6 + s7);

    delete fe;
  }
}

void Elasticity::write_data(const string & filename)
{
  arma::vec vmstress(msh.get_n_elements());
  std::string vtuname = filename + ".vtu";
  std::string vename = filename + "_ve";

  // copy solution to displacements array
  // (u1, u2, ... un, v1, v2, ...vn) -> (u1, v1, u2, v2, ...)
  int np = msh.get_n_points();
  int nd = msh.get_n_dim();
  arma::vec displacements(3*np, arma::fill::zeros);

  for(int d=0; d<nd; d++)
    for(int i=0; i<np; i++)
      displacements(d + i * 3) = solution(i + d * np);

  // write output data
  writer.open(filename, 1, 1);
  writer.write_displ_step(0, displacements.memptr());

  //von_mises(solution,vmstress);
  //vtkout.write_def_mesh(displacement, vmstress, vtuname);
}
