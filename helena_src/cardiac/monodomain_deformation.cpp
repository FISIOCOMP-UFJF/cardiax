#include "monodomain_deformation.hpp"

MonodomainDeformation::MonodomainDeformation() :
        Monodomain(), condtype(S_ISOTROPIC)
{
  cout << "Monodomain Deformation" << endl;
  mesh = new Mesh();
  writer = new WriterHDF5(mesh);

  // *************************
  // Hardcoded
  // *************************
  //condtype = M_ISOTROPIC;
  //condtype = S_ISOTROPIC;
  //condtype = M_TRANSVERSE;
  //condtype = S_TRANSVERSE;
  //condtype = M_ORTHOTROPIC;
  //condtype = S_ORTHOTROPIC;
  // *************************
}

void MonodomainDeformation::advance(const std::vector<arma::mat33*> & vec_F)
{
  //static int step=1;

  if(vec_F.empty()) error("vecF is empty");

  if( !tip.finished() )
  {
    tip.increase_time();

    timer.enter("ODEs");
    solve_odes();
    timer.leave();

    timer.enter("Parabolic");
    solve_parabolic(vec_F);
    timer.leave();

    //write_data_text(vm, step);
    //write_data(vm, displ, "vm_vtu/vm", &step);
  }
}

void MonodomainDeformation::assemble_matrices(const std::vector<arma::mat33*> & vec_F)
{
  bool test = false;
  const int n = mesh->get_nen();
  const int nz = parameters["maxnz"];
  const double um2_to_cm2  = 1.0e-8;
  const double theta = parameters["theta_method"];
  const double surf_to_vol = parameters["surface_to_volume"];
  const double kappa = (surf_to_vol/timestep) * um2_to_cm2;
  const double dtkappa = 1.0/kappa;
  arma::mat elmat_k(n,n), elmat_m(n,n);
  std::vector<int> dnums;

  if(Ai.is_null())
  {
    std::cout << "\nCreating PETSc-FEM matrices\n";
    Ai.create(ndofs, ndofs, nz);
    Mi.create(ndofs, ndofs, nz);
    v0.create(ndofs);
    v1.create(ndofs);
    f.create(ndofs);
    test = true;
  }
  else
  {
    Ai = 0.0;
    Mi = 0.0;
  }
  
  FiniteElement & fe = fespace.createFE(0);

  for(int e=0; e<mesh->get_n_elements(); e++)
  { 
    calc_elmat_stiff (e, fe, elmat_k, vec_F);
    calc_elmat_mass  (e, fe, elmat_m);
    fespace.get_element_dofs(e,dnums);

    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
      {
        double aval = elmat_m(i,j) + theta * dtkappa * elmat_k(i,j);
        double mval = elmat_m(i,j) - theta * dtkappa * elmat_k(i,j);
        int ii = dnums[i];
        int jj = dnums[j];
        Ai.add(ii, jj, aval);
        Mi.add(ii, jj, mval);
      }
  }

  delete &fe;

  Ai.assemble();
  Mi.assemble();

  if(test) solver.init();
}

void MonodomainDeformation::calc_elmat_stiff(const int elid, 
                                             const FiniteElement & fe,
                                             arma::mat & elmat,
                                             const ArrayMat33 & vecF)
{
    int ndim = fe.get_ndim();
    int ndof = fe.ndofs();
    double detJxW;
    arma::mat dshape;
    arma::mat gradn(ndof,ndim);
    arma::mat jacinv(ndim,ndim);
    arma::mat sigma(ndim,ndim);
    arma::mat33 *F;

    // fibers direction
    const arma::vec3 & f0 = get_fiber(elid);
    const arma::vec3 & s0 = get_trans(elid);
    const arma::vec3 & n0 = get_normal(elid);

    Quadrature * qd = Quadrature::create(fe.order(),fe.type());
    Mapping em = fe.get_mapping(elid);

    int nint = qd->get_num_ipoints();
    elmat.zeros();

    calc_cond_tensor(elid, ndim, sigma);

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
        fe.calc_deriv_shape(qd->get_point(q), dshape);
        em.calc_jacobian(dshape);
        detJxW = qd->get_weight(q) * em.get_det_jacobian();
        jacinv = em.get_inv_jacobian();
        gradn = dshape * jacinv; // shape function deriv

        // Get pointer to deformation gradient tensor
        F = vecF[elid*nint + q];

        // 'New' conductivity ----------------------
        update_conductivity(*F, f0, s0, n0, sigma);
        // -----------------------------------------

        // Compute element matrix
        elmat += detJxW * (gradn * sigma) * gradn.t();
    }

    delete qd;
}
/*
void MonodomainDeformation::init()
{
    tip = TimeParameters(timestep, totaltime, printrate);

    mesh->read(mesh_filename);
    stimuli.read(stimuli_filename);

    fespace.set_mesh(mesh);
    ndofs = mesh->get_n_points();

    // setup data writer to write at every 1 ms
    // potential field and displacements
    std::string output = mesh_filename + "_output";
    int nsteps = tip.get_nsteps()/(1.0/timestep);
    writer->open(output, nsteps+1, timestep);

    // setup model and cells
    cellmodel = CellModel::create(cell_name);
    cellmodel->setup(odesolver, timestep, totaltime, 1.0);
    cells = new Cells(ndofs,cellmodel);

    vm.resize(ndofs);

    timer.enter("Assemble");
    Monodomain::assemble_matrices();
    timer.leave();
}
*/

void MonodomainDeformation::init()
{
    tip = TimeParameters(timestep, totaltime, printrate);

    mesh->read_xml(mesh_filename);
    stimuli.read_xml(stimuli_filename);

    fespace.set_mesh(mesh);
    ndofs = mesh->get_n_points();

    // setup data writer to write at every 1 ms
    // potential field and displacements
    std::size_t pos  = mesh_filename.find(".xml");
    std::string output = mesh_filename.substr(0,pos) + "_output";
    int nsteps = tip.get_size(); //tip.get_nsteps();// /(1.0/timestep); //troquei para conseguir rodar casos maiores (juvs)
    writer->open(output, nsteps+1, timestep);

    // setup model and cells
    cellmodel = CellModel::create(cell_name);
    cellmodel->setup(odesolver, timestep, totaltime, 1.0);
    cells = new Cells(ndofs,cellmodel);

    vm.resize(ndofs);

    timer.enter("Assemble");
    Monodomain::assemble_matrices();
    timer.leave();
}


void MonodomainDeformation::set_conductivity(int cond)
{
    CondTensorType tcond = static_cast<CondTensorType>(cond);
    condtype = tcond;
  
    cout << "Conductivity type: ";
    switch(condtype)
    {
        case S_ISOTROPIC:   cout << "S_ISOTROPIC"   << endl; break;
        case S_TRANSVERSE:  cout << "S_TRANSVERSE"  << endl; break;
        case S_ORTHOTROPIC: cout << "S_ORTHOTROPIC" << endl; break;
        case M_ISOTROPIC:   cout << "M_ISOTROPIC"   << endl; break;
        case M_TRANSVERSE:  cout << "M_TRANSVERSE"  << endl; break;
        case M_ORTHOTROPIC: cout << "M_ORTHOTROPIC" << endl; break;
    }
}

void MonodomainDeformation::solve_parabolic(const ArrayMat33 & vec_F)
{
  const double pcgtol = parameters["pcgtol"];

  // if needed, re-assemble the stiffness matrix
  if(re_assembly_mats)
  {
    assemble_matrices(vec_F);
    re_assembly_mats = false;
  }

  //
  // b = M * vm0 (sparse matrix vector multiplication)
  // A * v1 = b (solve)
  //
  std::pair<PetscInt,PetscReal> ir;
  Mi.mult(v0, f);
  ir = solver.solve(Ai, v1, f, pcgtol);

  if (tip.time2print())
    cout << "Time: " << tip.time() << " num its " << ir.first << " rnorm " << ir.second << endl;

  // copy solution from PETSc Vec to my Vector
  v1.get_data(vm.memptr());
  cells->set_var(0, vm);
}

void MonodomainDeformation::update_conductivity(const arma::mat33 & F, 
                                                const arma::vec3 & f0,
                                                const arma::vec3 & s0,
                                                const arma::vec3 & n0,
                                                arma::mat & sigma)
{
  const double sigma_l = parameters["sigma_l"];
  const double sigma_t = parameters["sigma_t"];
  const double sigma_n = parameters["sigma_n"];

  double J = arma::det(F);

  if (condtype == S_ISOTROPIC)
  {
      arma::mat33 Finv = arma::inv(F);
      arma::mat33 Cinv = Finv * Finv.t();

      //sigma = J * sigma_l * Cinv;
      sigma = J * sigma_l * F * F.t();
  }
  else if (condtype == S_TRANSVERSE)
  {
      double diso = sigma_t;
      double dani = (sigma_l-sigma_t);

      // CHECK THIS-----------------------------------------
      //double stretchf = arma::norm(F*f0, 2);
      //double dani = (sigma_l-sigma_t)/(stretchf*stretchf);
      // ---------------------------------------------------

      arma::mat33 Finv = arma::inv(F);
      arma::mat33 Cinv = Finv * Finv.t();
      arma::mat33 Mfib = f0 * f0.t();

      sigma = J * (diso*Cinv + dani*Mfib);
  }
  else if (condtype == S_ORTHOTROPIC)
  {
      arma::mat33 R, U, Ds, Dm, Fi;
      arma::vec3 f, s, n;
      polar_decomposition(F, R, U);
      f  = R * f0;
      s  = R * s0;
      n  = R * n0;
      Fi = arma::inv(F);
      Ds = sigma_l * f * f.t() + sigma_t * s * s.t() + sigma_n * n * n.t();
      Dm = Fi * Ds * Fi.t();
      sigma = J * Dm;
  }
  else if (condtype == M_ISOTROPIC)
  {
      arma::mat33 I = arma::eye(3,3);
      sigma = J * sigma_l * I;
  }
  else if (condtype == M_TRANSVERSE)
  {
      double diso = sigma_t;
      double dani = (sigma_l-sigma_t);
      arma::mat33 I    = arma::eye(3,3);
      arma::mat33 Finv = arma::inv(F);
      arma::mat33 Mfib = f0 * f0.t();
      sigma = J * (diso*I + dani*Mfib);
  }
  else if (condtype == M_ORTHOTROPIC)
  {
      //double stretchf = arma::norm(F*f0, 2);
      //double stretchs = arma::norm(F*s0, 2);
      //double stretchn = arma::norm(F*n0, 2);
      //double da = sigma_l/(stretchf*stretchf);
      //double db = sigma_t/(stretchs*stretchs);
      //double dc = sigma_n/(stretchn*stretchn);
      //double da = sigma_l;
      //double db = sigma_t;
      //double dc = sigma_n;
      arma::mat33 Mfib = f0 * f0.t();
      arma::mat33 Mshe = s0 * s0.t();
      arma::mat33 Mnor = n0 * n0.t();
      //sigma = J * (da*Mfib + db*Mshe + dc*Mnor);

      sigma = J * (sigma_l*Mfib + sigma_t*Mshe + sigma_n*Mnor);
  }
}
