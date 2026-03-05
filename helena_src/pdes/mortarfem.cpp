#include "mortarfem.hpp"
#include "../util/timer.hpp"

namespace {
  const double penalty = 1.e+10;
}

void Mortarfem::assemble_system(int time_step)
{
  //Create typical 2D FE for the problem triang, quad and so on
  FiniteElement & fe = fespace.createFE(0);
  int n = fe.get_ndof(); //int n = msh.get_nen();
  int ndofs = fespace.get_num_dofs();
  //We need a method that return the number of interface elements in the finest side
  int num_interface_eleme = count_interface_elements();
  //here has to increase conforming the number of
  int mult_DGOfs = num_interface_eleme;
  int nTotal_dofs = ndofs + mult_DGOfs;
  puts("");puts("");

  printf("Systems Dim Total, Primary %d %d\n", nTotal_dofs, ndofs);
  arma::mat elmatd(n,n); //Defining mat local of diff
  arma::mat w(mult_DGOfs, mult_DGOfs); //Matrix to make conditioning
  arma::mat lambda_mat(mult_DGOfs, mult_DGOfs);


  arma::vec lambda_vec(mult_DGOfs);
  lambda_vec.zeros();

  // arma::vec Vm_vec(mult_DGOfs);
  // Vm_vec.zeros();
  if (time_step == 0 ){
    B.resize(ndofs,mult_DGOfs);
    B.zeros();
    Bt.resize(mult_DGOfs, ndofs);
    B.zeros();
    Vm_vec.resize(mult_DGOfs);
    Vm_vec.zeros();

    L.create(mult_DGOfs, mult_DGOfs, 30);
    f_lambda.create(mult_DGOfs);
    sol_lambda.create(mult_DGOfs);

    global_A.resize(ndofs, ndofs);
    global_A.zeros();
    primary_vec.resize(ndofs);
    primary_vec.zeros();

    K.create(nTotal_dofs, nTotal_dofs,30);
    u.create(nTotal_dofs);
    f.create(nTotal_dofs);
  }

  arma::vec elvec(n);
  vector<int> dnums;

  ZeroFunction<double> zerofunc;
  RHSExample1<double> rhsfunc;     // 2D

  // PETSC matrix and vector creation
  // This is going to be called every time loop important to avoid acomulation
  // f.create(nTotal_dofs);
  // lambda.create(mult_DGOfs);

  cout << "Assembling Block PETSc Matrix " << endl;
  if (time_step == 0 ){
    fill_primary_matrix_vector();

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
          global_A(I, J) += belmat(k,l);
        }
        for(int k=0;k<m;k++){
          f.add(bdnums[k],belvec(k));
          primary_vec(bdnums[k]) += belvec(k);
        }
      }
      delete &bfe;
    }
  }

  /////For matching grids!
  fill_interfaces_contributions_matching(time_step, ndofs, mult_DGOfs);

  ///for nonMatching meshes
  // fill_interfaces_contributions_non_matching(time_step, ndofs, mult_DGOfs);

  //This has to be called every time iteration cause load vector
  //of mult has to be recomputed every time iteration
  // Apply Dirichlet boundary conditions u = gD on \gamma_D

  if (time_step == 0){
    A_inv.resize(ndofs, ndofs);
    A_inv.zeros();
    w.eye();
    double b_norm = arma::norm(B);
    double omega = arma::norm(global_A);

    printf("Matrix Norm %f\n", omega);
    // w *= omega;
    w *= 1.0;
    arma::mat conditioning(ndofs, ndofs);
    conditioning = B*w*Bt;
    arma::vec conditioning_vec(ndofs);
    conditioning_vec = B*w*Vm_vec;
    A_inv = arma::inv(global_A + conditioning);
    lambda_mat = Bt*A_inv*B;
    lambda_vec = Bt*A_inv*(primary_vec + conditioning_vec) - Vm_vec;
    // lambda_vec = Bt*A_inv*primary_vec - Vm_vec;
    puts("Matrix condition number \n");
    printf("%f\n", arma::cond(lambda_mat));
    puts("");
    puts("");
    puts("");
    //From this matrix fill the L matrix to solve with PETSc lib
    puts("Solution of the multiplier from LU decompose");
    for (int i = 0; i < mult_DGOfs; i++) {
      for (int j = 0; j < mult_DGOfs; j++) {
        if (lambda_mat(i,j) != 0.0) {
          L.add(i,j,lambda_mat(i,j));
        }
      }
      f_lambda.add(i, lambda_vec(i));
    }
  }
  else{
    w.eye();
    w *= 1.0;
    arma::vec conditioning_vec(ndofs);
    conditioning_vec = B*w*Vm_vec;
    lambda_vec = Bt*A_inv*(primary_vec + conditioning_vec) - Vm_vec;
    for (int i = 0; i < mult_DGOfs; i++) {
      f_lambda.set(i, lambda_vec(i));
    }
  }
  // fix boundary nodes
  if (fixed_nodes_map.size()){
    printf("%s\n", "Imposing by fixed nodes" );
    FETools::apply_boundary_values(fixed_nodes_map, K, f);
  }

  K.assemble();
  L.assemble(); //may be needed to do only once
  // Schur = K.schur_complement(ndofs, ndofs, nTotal_dofs, nTotal_dofs);
  delete &fe;
}

void Mortarfem::computing_potential_diff( std::map<int, vector<int> > map_multiplier, int t_step){
  double Vm0, Vmt, Vmt_sol, Im, Cm, Rm;
  Cm = 1.0;
  Rm = 1.0;
  Vmt = 0.0;
  Vmt_sol = 0.0;
  int element_i = 0;
  int count = 0;

  double edo_dt;
  double g_dt; //To solve in micro seg scale is 1e-3 therefore the dt is 0.01!
  //it gives 10 time steps per deltaT of diff problem
  //////////////////////PARAMS OF HH
  double epsilon = 0.01;
  double omega = 0.5;
  double alfa = 0.1;
  double b = 0.5;
  double t = 0.02e-3; //Time at the evaluation for the exact sol
  //Every two keys is the value of interface voltage to be averaged by element

  /////////////////Constants for the Hudking-Huxley cell model
  double gk, gna, gl, nt, mt, ht, n0, m0, h0, vk, vna, vl;
  double alfa_m, beta_m, alfa_n, beta_n, alfa_h, beta_h;
  double k1v, k2v, k1n, k2n, k1m, k2m, k1h, k2h;
  vector<double> state_var(3);
  gk = 36.0;
  gna = 120.0;
  gl = 0.3;
  vk = -12.0;
  vna = 115.0;
  vl = 10.6;
  for (std::map< int , vector<int> >::iterator it=map_multiplier.begin(); it!=map_multiplier.end(); ++it){
    Im = lambda_sol(it->first); //Current is in scale of microAmp
    Vm0 = primal_solution(it->second[0]) - primal_solution(it->second[1]);

    if (t_step < 101){
      edo_dt = 0.005e-3;
      g_dt  = 0.05e-3; //To solve in micro seg scale is 1e-3 therefore the dt is 0.01!

      Vm0 *= 1000.0; //this is needed cause the system is already in units of mV
      if (t_step == 0){ //Here is when this solver Enters!, may be change to Zero
        m0 = 0.01;
        n0 = 0.3;
        h0 = 0.6;
      }
      else{
        m0 = interface_state_var[it->first][0];
        n0 = interface_state_var[it->first][1];
        h0 = interface_state_var[it->first][2];
        // printf("%f %f %f ", m0, n0, h0);
      }
      for (double t = 0.0; t <= g_dt; t += edo_dt) {
        alfa_m = 0.1*(25.0 - Vm0)/(exp((25.0-Vm0)/10.0) - 1.0);
        beta_m = 4.0*exp(-Vm0/18.0);
        alfa_n = 0.01*(10.0 - Vm0)/(exp((10.0-Vm0)/10.0) - 1.0);
        beta_n = 0.125*exp(-Vm0/80.0);
        alfa_h = 0.07*exp(-Vm0/20.0);
        beta_h = 1.0/(exp((30.0-Vm0)/10.0) + 1.0);

        // printf("%f %f %f %f %f %f ", alfa_m, beta_m, alfa_n, beta_n, alfa_h, beta_h);
        k1v = (-gk*n0*n0*n0*n0*(Vm0 - vk)
               -gna*m0*m0*m0*h0*(Vm0 - vna)
               -gl*(Vm0 - vl) + Im) * edo_dt;
        k2v = (-gk*n0*n0*n0*n0*(Vm0 + k1v*0.5 - vk )
               -gna*m0*m0*m0*h0*(Vm0 + k1v*0.5 - vna)
               -gl*(Vm0 + k1v*0.5 - vl) + Im ) * edo_dt;

        Vmt = k2v + Vm0;

        k1m = (alfa_m*(1.0 - m0) - beta_m*m0) * edo_dt;
        k2m = (alfa_m*(1.0 - m0 + 0.5*k1m) - beta_m*(m0 + 0.5*k1m)) * edo_dt;
        mt = k2m + m0;

        k1n = (alfa_n*(1.0 - n0) - beta_n*n0) * edo_dt;
        k2n = (alfa_n*(1.0 - n0 + 0.5*k1n) - beta_n*(n0 + 0.5*k1n)) * edo_dt;
        nt =  k2n + n0;

        k1h = (alfa_h*(1.0 - h0) - beta_h*h0) * edo_dt;
        k2h = (alfa_h*(1.0 - h0 + 0.5*k1h) - beta_h*(h0+ 0.5*k1h)) * edo_dt;
        ht =  k2h + h0;
        Vm0 = Vmt;
        m0 = mt;
        n0 = nt;
        h0 = ht;
      }
      state_var[0] = mt;
      state_var[1] = nt;
      state_var[2] = ht;
      interface_state_var[it->first] = state_var;
      interface_potential[it->first] = Vmt/1000.0; //needed to be scaled for the resutl!!
      // interface_potential[it->first] = Vmt/1000.0;
      printf(" %f ", Vmt/1000.0);

      // Vmt = (Vm0-Im*Rm)*exp(-1.0/(Cm*Rm) * t) + Im*Rm;
      // interface_potential[it->first] = Vmt;
      // printf(" %f ", Vmt);
    }
    else{
      edo_dt = 0.005;
      g_dt  = 0.05; //To solve in micro seg scale is 1e-3 therefore the dt is 0.01!
      // Vm0 = primal_solution(it->second[0]) - primal_solution(it->second[1]);
      ////Solve the action potential with inner cell voltaje!
      // Vm0 = primal_solution(it->second[0]);
      // Vm0 = primal_solution(48); //The node just in the middle of the cell
      Vm0 *= 1000.0; //this is needed cause the system is already in units of mV

      if (t_step == 1){ //Here is when this solver Enters!, may be change to Zero
        m0 = 0.01;
        n0 = 0.1;
        h0 = 0.6;
      }
      else{
        m0 = interface_state_var[it->first][0];
        n0 = interface_state_var[it->first][1];
        h0 = interface_state_var[it->first][2];
        // printf("%f %f %f ", m0, n0, h0);
      }
      for (double t = 0.0; t <= g_dt; t += edo_dt) {
        alfa_m = 0.1*(25.0 - Vm0)/(exp((25.0-Vm0)/10.0) - 1.0);
        beta_m = 4.0*exp(-Vm0/18.0);
        alfa_n = 0.01*(10.0 - Vm0)/(exp((10.0-Vm0)/10.0) - 1.0);
        beta_n = 0.125*exp(-Vm0/80.0);
        alfa_h = 0.07*exp(-Vm0/20.0);
        beta_h = 1.0/(exp((30.0-Vm0)/10.0) + 1.0);

        // printf("%f %f %f %f %f %f ", alfa_m, beta_m, alfa_n, beta_n, alfa_h, beta_h);
        k1v = (-gk*n0*n0*n0*n0*(Vm0 - vk)
               -gna*m0*m0*m0*h0*(Vm0 - vna)
               -gl*(Vm0 - vl) + Im) * edo_dt;
        k2v = (-gk*n0*n0*n0*n0*(Vm0 + k1v*0.5 - vk )
               -gna*m0*m0*m0*h0*(Vm0 + k1v*0.5 - vna)
               -gl*(Vm0 + k1v*0.5 - vl) + Im ) * edo_dt;

        Vmt = k2v + Vm0;

        k1m = (alfa_m*(1.0 - m0) - beta_m*m0) * edo_dt;
        k2m = (alfa_m*(1.0 - m0 + 0.5*k1m) - beta_m*(m0 + 0.5*k1m)) * edo_dt;
        mt = k2m + m0;

        k1n = (alfa_n*(1.0 - n0) - beta_n*n0) * edo_dt;
        k2n = (alfa_n*(1.0 - n0 + 0.5*k1n) - beta_n*(n0 + 0.5*k1n)) * edo_dt;
        nt =  k2n + n0;

        k1h = (alfa_h*(1.0 - h0) - beta_h*h0) * edo_dt;
        k2h = (alfa_h*(1.0 - h0 + 0.5*k1h) - beta_h*(h0+ 0.5*k1h)) * edo_dt;
        ht =  k2h + h0;
        Vm0 = Vmt;
        m0 = mt;
        n0 = nt;
        h0 = ht;
      }
      state_var[0] = mt;
      state_var[1] = nt;
      state_var[2] = ht;
      interface_state_var[it->first] = state_var;
      interface_potential[it->first] = Vmt/1000.0 - primal_solution(it->second[1]); //needed to be scaled for the resutl!!
      // interface_potential[it->first] = Vmt/1000.0;
      printf(" %f ", Vmt/1000.0);
    }
  }
}

void Mortarfem::config(std::string optfile)
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

void Mortarfem::calc_elmat_poisson(const int iel, const FiniteElement & fe,
                                 arma::mat & elmat)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  double detJxW;
  double alfa;
  double cond = 1.0;

  if (has_cond) cond = mesh.get_data(iel);

  Quadrature * qd = Quadrature::create(2*fe.order()-2,fe.type());
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat gradnv(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  //This method make mapping ref in mesh to local element
  Mapping em = fe.get_mapping(iel);
  elmat.zeros();
  if (mesh.get_element_index(iel) == 100){ //Inner Cell elemenets
    alfa = 5.0;
  }
  else if(mesh.get_element_index(iel) == 99){ //Outter Cell elements
    alfa = 20.0;
  }
  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);
    elmat += alfa * detJxW * cond * (gradn * gradn.t()); //2x3 . 3x2 operation
  }
  delete qd;
}

void Mortarfem::calc_elmat_mass(const int iel, const FiniteElement & fe,
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


void Mortarfem::calc_elvec_source (const int iel, const FiniteElement & fe,
																 const ScalarFunction<double> & rhs,
																 arma::vec & elvec)
{
  double detJxW, f;
  arma::vec shape;
  arma::mat dshape;
  Quadrature * qd = Quadrature::create(2, fe.type());
  Mapping em = fe.get_mapping(iel);
  std::vector<arma::vec3> geo_points;
  geo_points = em.get_points();
  elvec.zeros();
  //Here for each element is evaluated spatially the function to
  //add contribution of the global vector
  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q),shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    em.calc_jacobian(dshape);
    // f = rhs( em.map_point(qd->get_point(q), shape) );
    f = 0.0;
    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    elvec += (f * detJxW) * shape;
  }

  delete qd;
}

void Mortarfem::calc_interface_elmat (const int iel, const FiniteElement & fe,
																arma::mat & elmat)
{
  double detJxW;
  arma::vec shape;
  arma::mat dshape;
  int ndim = fe.get_ndim(); //DIM of the element
  int ndof = fe.ndofs();    //Ndof of the element
  arma::mat gradn(ndof,ndim);
  arma::mat gradnv(ndof,ndim);
  arma::vec3 jacinv;
  SurfaceMapping sm = fe.get_interface_mapping(iel);
  Quadrature * qd = Quadrature::create(2*fe.order(),fe.type());
  elmat.zeros();
  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q), shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    sm.calc_jacobian(dshape);
    detJxW = qd->get_weight(q) * sm.get_det_jacobian();
    double pdx = detJxW;
    elmat += pdx * (shape * shape.t());
  }
  delete qd;
}

void Mortarfem::calc_interface_elvec (const int iel, const FiniteElement & fe,
                                arma::vec & elvec, double Vm)
{
  double detJxW;
  arma::vec shape;
  arma::mat dshape;
  SurfaceMapping im = fe.get_interface_mapping(iel);
  Quadrature * qd = Quadrature::create(fe.order(), fe.type());
  elvec.zeros();
  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q),shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    im.calc_jacobian(dshape);
    arma::vec3 pt = qd->get_point(q);
    arma::vec3 xyz = im.map_point(pt[0]);
    detJxW = qd->get_weight(q) * im.get_det_jacobian();
    double coeff = Vm;//coeff_robin_mat(sm.get_index());
    double pdx = coeff * detJxW;
    elvec += pdx * shape;
  }
  delete qd;
}


void Mortarfem::calc_robin_elmat (const int iel, const FiniteElement & fe,
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

void Mortarfem::calc_robin_elvec (const int iel, const FiniteElement & fe,
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

void Mortarfem::calc_mass_matrix(arma::sp_mat & mass)
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
  delete &fe;
}

double Mortarfem::coeff_robin_mat(int index)
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

double Mortarfem::coeff_robin_vec(int index)
{
  std::map<int,double>::iterator it;

  //This function are hard coded to impose Homogeneous Neumann boundary conditions
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
  //Impose Neumann homogeneous boundary conditions
  // return 20.0;
  cout << " Boundary (Neumann) index " << index << " not set." << endl;
  exit(0);
}

void Mortarfem::init()
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
  fespace.set_mesh(&mesh);   //Already defined contain all information of the elements
  fespace.config();          //Initialize the number of degree of freedoms ndofs

  if (mesh.get_data_size())
    has_cond = true;

  cout << mesh;
}

//Function to writte exact solution in file for visualization
void Mortarfem::write_data(const string & filename, arma::vec & u_exact)
{
                 //nsteps, step
  writer.open("Sol_exact", 1, 1); //Maybe it is already set to write several times
  writer.setup();
  writer.write_vm_step(0, u_exact.memptr());
  writer.close();
}

void Mortarfem::write_data(const string & filename)
{
  // write output data
                 //nsteps, step
  writer.open(filename, 1, 1); //Maybe it is already set to write several times
  writer.setup();
  writer.write_vm_step(0, primal_solution.memptr());
  writer.close();
}

void Mortarfem::run(const string & name)
{
  filename = name;
  init();
  //Each counter is a time step of 0.02 microseconds
  for (int i = 0; i <= 0; i++) {
    assemble_system(i);
    solve(i);
    std::ostringstream sstream;
    sstream << "cell_pot_" << i;
    std::string result = sstream.str();
    //Function to write txt result!
    write_interface_potential(result);
    computing_potential_diff(map_multiplier, i); //Solving edo system
  }
  std::string name_gradient_file = "gradient";
  write_gradient(name_gradient_file);
}

void Mortarfem::solve(int t_step)
{
  Timer timer;
  int ndofs = fespace.get_num_dofs();
  int num_interface_eleme = count_interface_elements();
  // int mult_DGOfs = num_interface_eleme * 2;
  int mult_DGOfs = num_interface_eleme;
  int nTotal_dofs = ndofs + mult_DGOfs;

  std::pair<PetscInt,PetscReal> ir;
  std::pair<PetscInt,PetscReal> ir1;

  msg("Solving linear system");
  timer.start();
  static petsc::LinearSolver ls;
  static petsc::LinearSolver ls1; //only create once the solver!

  if (t_step == 0){
    ls1.init();
    ls.init();
    primal_solution.resize(ndofs);
    primary_sol.resize(nTotal_dofs);
    lambda_sol.resize(mult_DGOfs);
  }

  ls.use_umfpack();
  // ls.use_fieldsplit_schur(ndofs, nTotal_dofs);
  ir = ls.solve(K,u,f,1.e-20);
  // ir = ls.solve(K,u,f, 1.e-16);
  ir1 = ls1.solve(L,sol_lambda,f_lambda, 1.e-16);
  timer.stop();

  cout << " Number of iterations: " << ir1.first << endl;
  cout << " Residual norm: " << ir1.second << endl;
  cout << " Total time: " << timer.get_total_time() << endl;

  cout << " Number of iterations second solver: " << ir.first << endl;
  cout << " Residual norm: " << ir.second << endl;
  cout << " Total time: " << timer.get_total_time() << endl;

  // Copy solution from PETSc Vec to my Vector
  // sol_lambda.get_data(lambda_sol.memptr());
  // primal_solution = A_inv*(primary_vec - B*lambda_sol);

  arma::vec global_sol(nTotal_dofs);
  u.get_data(primary_sol.memptr());

  for (int i = 0; i < ndofs; i++) {
    primal_solution(i) = primary_sol(i);
  }
  puts("Printing the multiplier solution Schur condensate System");
  for (int i = 0; i < mult_DGOfs; i++) {
    printf(" %0.10f ", lambda_sol(i));
  }
  puts("");
  puts("");
  puts("");


  puts("Printing the multiplier solution from global system");
  for (int i = 0; i < mult_DGOfs; i++) {
    printf(" %0.10f ", primary_sol(i + ndofs));
    lambda_sol(i) = primary_sol(i + ndofs);
    f.set(i + ndofs, 0.0);
  }
  puts("");
  puts("");
  puts("");

}

///////////////////FUNCTION TO WRITTE SOLUTION IN TIME OF POTENTIAL
void Mortarfem::write_interface_potential(string file_name)
{
  ofstream myfile_in;
  ofstream myfile_out;
  ofstream myfile_pot;

  std::string name_1 = "/home/cesar/Desktop/finit_elem/trab_tesis/HH_cell_potential/in" + file_name;
  std::string name_2 = "/home/cesar/Desktop/finit_elem/trab_tesis/HH_cell_potential/out" + file_name;
  std::string name_3 = "/home/cesar/Desktop/finit_elem/trab_tesis/HH_cell_potential/pot" + file_name;

  myfile_in.open(name_1.c_str());
  myfile_out.open(name_2.c_str());
  myfile_pot.open(name_3.c_str());

  //Here the angle in Degree is defined for post-processing porposes
  double x, y, r, cos_theta, angle;
  arma::vec3 points_coordinates_in;

  for (std::map< int , vector<int> >::iterator it=map_multiplier.begin(); it!=map_multiplier.end(); ++it){
    //here by default taking the inner point coordinates
    points_coordinates_in = mesh.get_point(it->second[0]);
    x = points_coordinates_in[0]-0.005;
    y = points_coordinates_in[1]-0.005;
    r = sqrt(x*x + y*y);
    cos_theta = x/r;
    if ( y > 0 ){
      angle = acos(x/r) * 180.0/M_PI;
    }
    else if (x < 0){
      angle = acos(abs(x)/r) * 180.0/M_PI + 180.0;
    }
    else{
      angle = 360.0 - acos(x/r) * 180.0/M_PI;
    }
    myfile_in << angle;
    myfile_in << ",";
    myfile_in << primal_solution[it->second[0]];
    myfile_in << "\n";

    myfile_out << angle;
    myfile_out << ",";
    myfile_out << primal_solution[it->second[1]];
    myfile_out << "\n";

    myfile_pot << angle;
    myfile_pot << ",";
    myfile_pot << primal_solution[it->second[0]] - primal_solution[it->second[1]]; //computing potential diff
    myfile_pot << "\n";
  }
  myfile_in.close();
  myfile_out.close();
  myfile_pot.close();
}

int Mortarfem::count_interface_elements(){
  int count97 = 0;
  int count98 = 0;
  for(int i=0; i < mesh.get_n_interface_elements(); i++)
  {
    //Save in a map structure the index
    if (mesh.get_interface_element_index(i) == 97){
      count97 += 1;
    }
    if(mesh.get_interface_element_index(i) == 98){
      count98 += 1;
    }
  }
  // count98 >= count97 ? return count98 : return count97;
  if (count97 <= count98){
    return count97;
  }
  else{
    return count98;
  }
}
///Fill Block mat A and
void Mortarfem::fill_primary_matrix_vector()
{
  FiniteElement & fe = fespace.createFE(0);
  int n = fe.get_ndof();
  arma::mat elmatd(n,n);
  arma::vec elvec(n);
  vector<int> dnums;
  // ZeroFunction<double> zerofunc;
  RHSExample1<double> rhsfunc;
  for(int i=0; i < mesh.get_n_elements(); i++)
  {
    fespace.get_element_dofs(i,dnums); //Here define which nodes are of each element
    calc_elmat_poisson(i, fe, elmatd); // local diff matrix
    calc_elvec_source(i, fe, rhsfunc, elvec);

    for(int k=0; k<n; k++) //dimension size x
    for(int l=0;l<n;l++) //dimension size y
    {
      int I = dnums[k];
      int J = dnums[l];
      K.add(I, J, elmatd(k,l));
      global_A(I,J) += elmatd(k,l);
    }
    for(int k=0; k<n; k++){
      f.add(dnums[k], elvec(k));
      primary_vec(dnums[k]) += elvec(k);
    }
  }
  delete &fe;

}

////This method make multiplier continous
////Nevertheless matrix initialization has to be
///Acording to chossen mult space
///Matching and non matching (for refined inner cell)
/////////BEWARE FOR CELL CLUSTER THE SYSTEM BLOWSUP
 void Mortarfem::fill_interfaces_contributions_non_matching(int time_step, int ndofs, int num_interface_eleme)
 {
   ///This is needed to refill the load vector in multiplier side
   Vm_vec.zeros();
   std::map<int, int> in_interface;
   std::map<int, int> out_interface;
   FiniteElement & ife = fespace.create_interface_FE(0);
   int m = ife.get_ndof();
   arma::mat belmat(m,m);
   arma::vec belvec(m);
   vector<int> idnums_i;
   vector<int> idnums_o;
   vector<int> hold(2); //This guy is to relate Multiplier with
   double Vm;
   for(int i=0; i < mesh.get_n_interface_elements(); i++)
   {
     if (mesh.get_interface_element_index(i) == 97){
       in_interface[i] = i;
     }
     if(mesh.get_interface_element_index(i) == 98){
       out_interface[i] = i;
     }
   }

   int i = 0;
   for (std::map< int, int >::iterator it=in_interface.begin(); it!=in_interface.end(); ++it)
   {
     fespace.get_interface_element_dofs(it->first, idnums_i);
     calc_interface_elmat(it->first, ife, belmat);
     if (time_step == 0) { //Initial condition
       // t = 1.0;
       // Vm = E*15e-4*(1.0-epsilon)*(cos_theta1+cos_theta2)/2.0*(1.0-exp(-t/tau));
       Vm = 0.0; //To test Convergence!
       for(int k=0;k<m;k++) //dimension size y
       {
         for(int l=0;l<m;l++) //dimension size x
         {
           // int I_i = i*2 + k + ndofs;    //idnums[k];
           fespace.get_interface_element_dofs(2*i + l + num_interface_eleme, idnums_o);
           int I_i = i + k + ndofs;    //idnums[k];
           int J_i = idnums_i[l];        //-> makes reference to Primary DGOFs
           int J_o, I_o;
           J_o = idnums_o[l];
           if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
             I_i = ndofs;
           }
           K.add(I_i, J_i, belmat(k,l));
           K.add(I_i, J_o, -belmat(k,l));
           I_i = J_i;         //beware it was k originally
           J_i = i + k + ndofs;     //idnums[l];//-> makes reference to Primary DGOFs
           I_o = J_o;
           if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
             J_i = ndofs;
             J_o = ndofs;
           }
           K.add(I_i, J_i, belmat(k,l) );
           K.add(I_o, J_i, -belmat(k,l) );
           ///0, 1
           //It is exclusive for non-matching grids taken 2 by 2
           I_i = i + k;
           J_i = idnums_i[l];
           J_o = idnums_o[l];
           if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
             I_i = 0;
           }
          //  Bt(I_i, J_i) += belmat(k,l);
          //  Bt(I_i, J_o) += -belmat(k,l);
           I_i = J_i;
           J_i = i + k;
           I_o = J_o;
           if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
             J_i = 0;
             J_o = 0;
           }
           B(I_i, J_i) += belmat(k,l);
           B(I_o, J_i) += -belmat(k,l);

           hold[0] = idnums_i[l];
           hold[1] = idnums_o[l];
           if (i != num_interface_eleme && l != 1){ //Means is the last element
             map_multiplier[i + l] = hold;
           }
         }
       }
       //Define matrix transpose
       Bt = B.t();
     }
     else{ //Solution of the EDO
       if (i != num_interface_eleme-1){
         Vm = (interface_potential[i] + interface_potential[i + 1])/2.0;
       }
       else{
         Vm = (interface_potential[i] + interface_potential[0])/2.0;
       }
     }
     calc_interface_elvec(it->first, ife, belvec, Vm);
     for(int k=0; k<m; k++){
       if (i == num_interface_eleme - 1 && k == 1){ //Means is the last element
         f.add(ndofs, belvec(k));
         Vm_vec(0) += belvec(k);
       }
       else{
         f.add(i + k + ndofs, belvec(k));
         Vm_vec(i + k) += belvec(k);
       }
     }
     i += 1;
   }
 }

//////This method is to fill interface contributions
//////for matching grids!
void Mortarfem::fill_interfaces_contributions_matching(int time_step, int ndofs, int num_interface_eleme)
{

  Vm_vec.zeros();

  FiniteElement & ife = fespace.create_interface_FE(0);
  int m = ife.get_ndof();
  arma::mat belmat(m,m);
  arma::vec belvec(m);
  vector<int> idnums_i;
  vector<int> idnums_o;
  vector<int> hold(2);
  double Vm;

  std::map<int, int> in_interface;
  std::map<int, int> out_interface;
  for(int i=0; i < mesh.get_n_interface_elements(); i++)
  {
    if (mesh.get_interface_element_index(i) == 97){
      in_interface[i] = i;
    }
    if(mesh.get_interface_element_index(i) == 98){
      out_interface[i] = i;
    }
  }
  int i = 0;
  for (std::map< int, int >::iterator it=out_interface.begin(); it!=out_interface.end(); ++it)
  {
    fespace.get_interface_element_dofs(it->first, idnums_o);
    fespace.get_interface_element_dofs(i, idnums_i);
    calc_interface_elmat(it->first, ife, belmat);
    if (time_step == 0) { //Initial condition
      // t = 1.0;
      // Vm = E*15e-4*(1.0-epsilon)*(cos_theta1+cos_theta2)/2.0*(1.0-exp(-t/tau));
      Vm = 0.0; //To test Convergence!
      for(int k=0;k<m;k++) //dimension size y
      {
        hold[0] = idnums_i[k];
        hold[1] = idnums_o[k];
        if (i != num_interface_eleme && k != 1){ //Means is the last element
          map_multiplier[i + k] = hold;
        }
        for(int l=0;l<m;l++) //dimension size x
        {
          int I_i = i + k + ndofs;    //idnums[k];
          int J_i = idnums_i[l];        //-> makes reference to Primary DGOFs
          int J_o, I_o;
          J_o = idnums_o[l];
          if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
            I_i = ndofs;
          }
          K.add(I_i, J_i, belmat(k,l));
          K.add(I_i, J_o, -belmat(k,l));
          I_i = J_i;         //beware it was k originally
          J_i = i + k + ndofs;     //idnums[l];//-> makes reference to Primary DGOFs
          I_o = J_o;
          if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
            J_i = ndofs;
            J_o = ndofs;
          }
          K.add(I_i, J_i, belmat(k,l) );
          K.add(I_o, J_i, -belmat(k,l) );

          I_o = i + k;
          J_o = idnums_o[l];
          J_i = idnums_i[l];
          if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
            I_o = 0;
          }
          // Bt(I_o, J_o) += -belmat(k,l);
          // Bt(I_o, J_i) += belmat(k,l);
          I_o = J_o;
          J_o = i + k;
          I_i = J_i;
          if (i == num_interface_eleme-1 && k == 1){ //Means is the last element
            J_o = 0;
          }
          B(I_o, J_o) += -belmat(k,l);
          B(I_i, J_o) += belmat(k,l);

          /////////////Solving with a penalisation in diagonal of C
          double omega = 1.e-16;
          I_i = i + ndofs;
          K.add(I_i, I_i, omega);
        }
      }
      Bt = B.t();
    }
    else{ //Solution of the EDO
      if (i != num_interface_eleme-1){
        Vm = (interface_potential[i] + interface_potential[i + 1])/2.0;
      }
      else{
        Vm = (interface_potential[i] + interface_potential[0])/2.0;
      }
    }
    calc_interface_elvec(it->first, ife, belvec, Vm);
    for(int k=0; k<m; k++){
      if (i == num_interface_eleme - 1 && k == 1){ //Means is the last element
        f.add(ndofs, belvec(k));
        Vm_vec(0) += belvec(k);
      }
      else{
        f.add(i + k + ndofs, belvec(k));
        Vm_vec(i + k) += belvec(k);
      }
    }
    i += 1;
  }
}

////Function to compute the gradient of the solution for each element
////For visualization in paraview
void Mortarfem::write_gradient(const string & filename)
{
  FiniteElement & fe = fespace.createFE(0);
  int n = fe.get_ndof();                //Number of nodes per element
  int ndim = fe.get_ndim();             //Spatial dimension
  int ndofs = fespace.get_num_dofs();    //Number total of dgof os the mesh
  // int nelement = mesh.get_n_elements();
  double detJxW;

  vector<int> dnums;

  arma::vec pos_vect(n, arma::fill::ones);
  arma::vec gradient(ndofs, arma::fill::zeros);

  arma::mat dshape;
  arma::vec grad(ndim);  //2 mat  linear triang
  arma::vec diff(ndim);  //2 mat  linear triang
  arma::vec sol_u(n);          //3X1 vect linear triang
  arma::mat jacinv(ndim,ndim);

  puts("");puts("");puts("");puts("");
  for(int iel=0; iel < mesh.get_n_elements(); iel++)
  {
    Mapping em = fe.get_mapping(iel);
    fespace.get_element_dofs(iel,dnums); //Here define which nodes are of each element

    for (int i = 0; i < n; i++) {
      sol_u(i) = primal_solution(dnums[i]);
      printf("%f ", sol_u(i));
    }
    puts("");
    ///Derivative of the function may have Spatial Dependance
    ///For linear triangular elements the gradient of the
    ///interpol functions is constant!
    fe.calc_deriv_shape(pos_vect, dshape);
    em.calc_jacobian(dshape);
    jacinv = em.get_inv_jacobian();

    // diff(0) = -sol_u(0) + sol_u(1);
    // diff(1) = -sol_u(0) + sol_u(2);
    grad = jacinv * (dshape.t() * sol_u);
    // grad = jacinv * diff;
    // printf("%f %f", grad(0), grad(1));
    // printf("%f %f %f", dshape(0,0), dshape(0,1), dshape(0,2));
    // printf("%f %f %f", dshape(1,0), dshape(1,1), dshape(1,2));

    // for(int d=0; d<ndim; d++)      //2
    for (int i = 0; i < n; i++)  //3
      gradient(dnums[i]) = sqrt(grad(0)*grad(0) + grad(1)*grad(1));

  }

  writer.open(filename, 1, 1);
  writer.setup();
  writer.write_vm_step(0, gradient.memptr());
  writer.close();
}


////////This is in case it is needed to compute
////////the transmembrane potential with spatial dependance
//     double x1, y1, x2, y2, r1, cos_theta1, cos_theta2, t;
//     double Sig_E, Sig_I, Cm, Rm, tau, epsilon, E;
//     arma::vec3 points_coordinates1;
//     arma::vec3 points_coordinates2;
//     E = 5.0;
//     Sig_E = 20.0;
//     Sig_I = 5.0;
//     Rm = 1.0;
//     Cm = 1.0;
//     tau = 1.0/(1.0/(Rm*Cm*1000.0) + 2.0*Sig_E*Sig_I/(Cm*1.5*(Sig_E+Sig_I)) );
//     epsilon = tau/(Cm*Rm*1000.0);

//     //   fespace.get_interface_element_dofs(it->first, idnums_o);
//     //   fespace.get_interface_element_dofs(i, idnums_i);
//     //   points_coordinates1 = mesh.get_point(idnums_i[0]);
//     //   points_coordinates2 = mesh.get_point(idnums_i[1]);
//     //   x1 = points_coordinates1[0]-0.005;
//     //   y1 = points_coordinates1[1]-0.005;
//     //   x2 = points_coordinates2[0]-0.005;
//     //   y2 = points_coordinates2[1]-0.005;
//     //   r1 = sqrt(x1*x1 + y1*y1);
//     //   cos_theta1 = x1/r1;
//     //   cos_theta2 = x2/r1;
//     //   calc_interface_elmat(it->first, ife, belmat);
//     //   if (time_step == 0) { //Initial condition
//     //     // t = 1.0;
//     //     // Vm = E*15e-4*(1.0-epsilon)*(cos_theta1+cos_theta2)/2.0*(1.0-exp(-t/tau));
//     //     Vm = 0.0; //To test Convergence!