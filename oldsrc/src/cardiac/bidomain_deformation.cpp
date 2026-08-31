#include "bidomain_deformation.hpp"

BidomainDeformation::BidomainDeformation() : Bidomain(), condtype(S_ISOTROPIC)
{
  cout << "Bidomain Deformation" << endl;
  mesh = new Mesh();
  writer = new WriterHDF5(mesh);
}

void BidomainDeformation::set_conductivity(int cond)
{
  CondTensorType tcond = static_cast<CondTensorType>(cond);
  condtype = tcond;
  
  cout << "Conductivity type (BidomainDeformation): ";
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

void BidomainDeformation::advance(const ArrayMat33 & vec_F)
{
  static int vmstep=0;
  static int vm_apd_step=0;
  static int vestep=0;

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
    
    timer.enter("Elliptic");
    solve_elliptic(vec_F);
    timer.leave();

    write_data_text(tmp_vm, &vm_apd_step);
    write_data(tmp_vm, "vm_vtu/vm", &vmstep);
    write_data(tmp_ve, "ve_vtu/ve", &vestep);
    
    //write_vtk_vm(tmp_vm, "output/vm_vtu/vm_tissue",&vmstep);
    //write_vtk_ve(tmp_ve, "output/ve_vtu/ve", &vestep);
   
    if (tip.time2print()) cout << endl;
  }
}

void BidomainDeformation::assemble_matrices(const ArrayMat33& vec_F)
{
  const int n = mesh->get_nen();
  const int nz = parameters["maxnz"];
  const double um2_to_cm2  = 1.0e-8;
  const double theta = parameters["theta_method"];
  const double surf_to_vol = parameters["surface_to_volume"];
  const double kappa = (surf_to_vol/timestep) * um2_to_cm2;
  const double dtkappa = 1.0/kappa;
  
  if(Ai.is_null())
  {
    cout << endl << "Assembling FEM matrices" << endl;
  
    // matrices
    Aie.create(ndof,ndof,nz);
    Ai.create(ndof_tissue,ndof_tissue,nz);
    Mi.create(ndof_tissue,ndof_tissue,nz);
    Ki.create(ndof_tissue,ndof_tissue,nz);

    // vectors
    vm0.create(ndof_tissue);
    vm1.create(ndof_tissue);
    b0.create(ndof_tissue);
    b1.create(ndof_tissue);
    aux_vec.create(ndof_tissue);
    ve.create(ndof);
    be.create(ndof);    
  }
  else
  {
    // reset matrices
    Aie = 0.0;
    Ai  = 0.0;
    Mi  = 0.0;
    Ki  = 0.0;
  }

  arma::mat elmat_ki(n,n), elmat_kie(n,n), elmat_m(n,n), elmat_b(n,n);
  std::vector<int> dnums;

  FiniteElement & fe = fespace.createFE(0);

  for(int e=0; e<mesh->get_n_elements(); e++)
  {
    calc_elmat_stiff_ie (e, fe, elmat_ki, elmat_kie, vec_F);
    calc_elmat_mass (e, fe, elmat_m);
    get_dof_nums (e,dnums);

    // tissue
    if (mesh->get_element_index(e) != 0)
    {
      // assembly of elliptic problem matrix
      for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
          Aie.add(dnums[i], dnums[j], elmat_kie(i,j));

      // assembly of parabolic problem matrices
      for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
        {
          double aval = elmat_m(i,j) + theta * dtkappa * elmat_ki(i,j);
          double mval = elmat_m(i,j) - theta * dtkappa * elmat_ki(i,j);
          double kval = timestep * elmat_ki(i,j);
          Ai.add(dnums[i], dnums[j], aval);
          Mi.add(dnums[i], dnums[j], mval);
          Ki.add(dnums[i], dnums[j], kval);
        }
    }
    // bath
    else
    {
      calc_elmat_bath(e, fe, elmat_b);

      // assembly of elliptic problem matrix
      for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
          Aie.add(dnums[i], dnums[j], elmat_kie(i,j));
    }
  }

  // FIX FIX FIX
  const double p = 1.e+10;
  if(ndof_bath != 0) 
    Aie.add(ndof_bath-1,ndof_bath-1,p);
  else
    Aie.add(0,0,p);

  Mi.assemble();
  Ai.assemble();
  Ki.assemble();
  Aie.assemble();

  cout << "End of FEM matrices assembly" << endl;

  par_solver.init();
  ell_solver.init();

  delete &fe;
}

void BidomainDeformation::calc_elmat_stiff_ie(const int iel, 
                                              const FiniteElement & fe,
                                              arma::mat & elmat_i,
                                              arma::mat & elmat_ie,
                                              const ArrayMat33 & vecF)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  int nint;
  double detJxW;

  Quadrature * qd = Quadrature::create(fe.order(),fe.type());
  nint = qd->get_num_ipoints();
  
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat sigma_i(ndim,ndim);
  arma::mat sigma_e(ndim,ndim);
  arma::mat sigma_aux(ndim,ndim);
  arma::mat33 *F;
  
  // fibers direction
  const arma::vec3 & f0 = get_fiber (iel);
  const arma::vec3 & s0 = get_trans (iel);
  const arma::vec3 & n0 = get_normal(iel);

  Mapping em = fe.get_mapping(iel);

  elmat_i.zeros();
  elmat_ie.zeros();   

  calc_cond_tensors(iel,ndim,sigma_i,sigma_e);

  for(int q=0; q<nint; q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();

    gradn = dshape * jacinv;

    // Get pointer to deformation gradient tensor
    F = vecF[iel*nint + q];
    
    // 'New' conductivity ----------------------
    sigma_aux = sigma_i + sigma_e;
    update_conductivity(*F, f0, s0, n0, sigma_i, sigma_aux);
    // -----------------------------------------
    
    elmat_i  += detJxW * (gradn * sigma_i  ) * gradn.t();
    elmat_ie += detJxW * (gradn * sigma_aux) * gradn.t();
  }

  delete qd;
}

void BidomainDeformation::solve_parabolic(const ArrayMat33 & vec_F)
{
  const double pcgtol = parameters["pcgtol"];

  // if needed, re-assemble the stiffness matrix
  if(re_assembly_mats)
  {
    cout << "Reassembling stiffness matrix" << endl;
    assemble_matrices(vec_F);
    re_assembly_mats = false;
  }  
    
  // Solve
  // b0 = M * vm0 (sparse matrix vector multiplication)
  // b = Ki * ve0
  // A * vm1 = b (solve)

  std::pair<PetscInt,PetscReal> ir;

  // NO BATH
  //Mi.mult(vm0, b0);
  //Ki.mult(ve, b1);

  // BATH
  Mi.mult(vm0, b0);
  aux_vec.copy_values(ndof_tissue,ve);
  Ki.mult(aux_vec, b1);

  b1.add_vec(b0);

  ir = par_solver.solve(Ai, vm1, b1,  pcgtol);

  if (tip.time2print())
  {
    cout << "Time " << tip.time();
    cout << " Parabolic its " << ir.first;
    cout << " rnorm " << scientific << ir.second;
  }

  // Copy solution from PETSc Vec to my Vector
  vm1.get_data(tmp_vm.memptr());

  // Copy data from tmp_vm to state variables array
  cells->set_var(0,tmp_vm);
}

void BidomainDeformation::solve_elliptic(const ArrayMat33 & vec_F)
{
  const double pcgtol = parameters["pcgtol"];

  // if needed, re-assemble the stiffness matrix
  if(re_assembly_mats)
  {
    assemble_matrices(vec_F);
    re_assembly_mats = false;
  }  
    
  // Solve
  // div((sigma_i + sigma_e) grad(Vm)) = - div(sigma_i grad(vm))
  std::pair<PetscInt,PetscReal> ir;

  // NO BATH
  //Ki.mult(vm1, be);
  //be.scale(-1.0);

  // BATH
  Ki.mult(vm1, aux_vec); // temp
  be.set(0.0); // form be = [be_tmp 0]
  be.copy_values(ndof_tissue, aux_vec);
  be.scale(-1.0);

  // FIX FIX FIX 
  if(ndof_bath != 0) 
    be.set(ndof_bath-1, 0.0);
  else
    be.set(0, 0.0);
    
  
  ir = ell_solver.solve(Aie, ve, be, pcgtol);

  if (tip.time2print())
  {
    cout << "  Elliptic its " << ir.first;
    cout << " rnorm " << ir.second;
  }

  // Copy solution from PETSc Vec to my Vector
  ve.get_data(tmp_ve.memptr());
}

void BidomainDeformation::update_conductivity(const arma::mat33 & F,
                                              const arma::vec3 & f0,
                                              const arma::vec3 & s0,
                                              const arma::vec3 & n0,
                                              arma::mat & sigma_i,
                                              arma::mat & sigma_ie)
{
  double J = arma::det(F);

  if (condtype == S_ISOTROPIC)
  {    
    arma::mat33 Finv = arma::inv(F);
    arma::mat33 Cinv = Finv * Finv.t();

    sigma_i  = J * sigma_il * F * F.t();
    sigma_ie = J * (sigma_il + sigma_el) * F * F.t();
  }
  else if (condtype == S_TRANSVERSE)
  {        
    double diso     = sigma_it;
    double diso_ie  = sigma_it + sigma_et;
    double dani     = (sigma_il-sigma_it);
    double dani_ie  = ((sigma_il+sigma_el)-(sigma_it+sigma_et));
    
    // CHECK THIS-----------------------------------------
    //double stretchf = arma::norm(F*f0, 2);
    //double dani = (sigma_l-sigma_t)/(stretchf*stretchf);
    // ---------------------------------------------------        
    
    arma::mat33 Finv = arma::inv(F);
    arma::mat33 Cinv = Finv * Finv.t();
    arma::mat33 Mfib = f0 * f0.t();
    
    sigma_i  = J * (diso    * Cinv + dani    * Mfib);
    sigma_ie = J * (diso_ie * Cinv + dani_ie * Mfib);
  }
  else if (condtype == S_ORTHOTROPIC)
  {
    double norm_Ff0 = arma::norm(F * f0);
    double norm_Fs0 = arma::norm(F * s0);
    double norm_Fn0 = arma::norm(F * n0);

    arma::mat33 Ds_i  = (sigma_il/(norm_Ff0*norm_Ff0)) * f0 * f0.t()
                      + (sigma_it/(norm_Fs0*norm_Fs0)) * s0 * s0.t()
                      + (sigma_in/(norm_Fn0*norm_Fn0)) * n0 * n0.t();

    arma::mat33 Ds_ie = ((sigma_il+sigma_el)/(norm_Ff0*norm_Ff0)) * f0 * f0.t()
                      + ((sigma_it+sigma_et)/(norm_Fs0*norm_Fs0)) * s0 * s0.t()
                      + ((sigma_in+sigma_en)/(norm_Fn0*norm_Fn0)) * n0 * n0.t();

    sigma_i  = Ds_i;
    sigma_ie = Ds_ie;

    /*
    arma::mat33 R, U, Ds_i, Ds_ie, Dm_i, Dm_ie, Fi;
    arma::vec3 f, s, n;
    polar_decomposition(F, R, U);
    f  = R * f0;
    s  = R * s0;
    n  = R * n0;
    Fi = arma::inv(F);

    Ds_i  = sigma_il * f * f.t()
          + sigma_it * s * s.t()
          + sigma_in * n * n.t();
    Ds_ie = (sigma_il+sigma_el) * f * f.t()
          + (sigma_it+sigma_et) * s * s.t()
          + (sigma_in+sigma_en) * n * n.t();

    Dm_i  = Fi * Ds_i  * Fi.t();
    Dm_ie = Fi * Ds_ie * Fi.t();

    sigma_i  = J * Dm_i;
    sigma_ie = J * Dm_ie;
    */
  }
  
  /*
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
    double stretchf = arma::norm(F*f0, 2);
    double stretchs = arma::norm(F*s0, 2);
    double stretchn = arma::norm(F*n0, 2);
    double da = sigma_l/(stretchf*stretchf);
    double db = sigma_t/(stretchs*stretchs);
    double dc = sigma_n/(stretchn*stretchn);
    arma::mat33 Mfib = f0 * f0.t();
    arma::mat33 Mshe = s0 * s0.t();
    arma::mat33 Mnor = n0 * n0.t();
    //sigma = J * (da*Mfib + db*Mshe + dc*Mnor);
    
    sigma = J * (sigma_l*Mfib + sigma_t*Mshe + sigma_n*Mnor);
  }*/
}
