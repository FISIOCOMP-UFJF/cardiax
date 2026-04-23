#include "monodomain.hpp"

/*!
 * Conductivities should be given in mS/um
 * 
 * My old setting: sigma_l(0.00012),   sigma_t(0.00006226), sigma_n(0.00002556)
 * Rodrigo's code: sigma_l(0.0001543), sigma_t(0.0000324),  sigma_n(0.0000120)
 * Benchmark     : sigma_l(0.0001334), sigma_t(0.0000176),  sigma_n(0.0000176)
 * Rabbit        : sigma_l(0.000204),  sigma_t(0.0000102),  sigma_n(0.000037)
 * Sundnes       : sigma_l(0.000300), sigma_t(0.000100), sigma_n(0.000031525)
*/

//#define MONO_ISOTROPIC

Monodomain::Monodomain()
  : CardiacProblem(),
    //sigma_l(0.0001334), sigma_t(0.0000176),  sigma_n(0.0000176),
    stim_apply_nodes(false)
{
  mesh = new Mesh();
  writer = new WriterHDF5(mesh);

  parameters.rename("monodomain_parameters");
  parameters.add("sigma_l", 0.0001334);
  parameters.add("sigma_t", 0.0000176);
  parameters.add("sigma_n", 0.0000176);
}

Monodomain::~Monodomain()
{
  delete cells;
  delete cellmodel;
}

void Monodomain::advance ()
{
  static int step = 0;
  //static int step_apd = 0;

  if( !tip.finished() )
  {
    tip.increase_time();

    timer.enter("ODEs");
    solve_odes();
    timer.leave();

    timer.enter("Parabolic");
    solve_parabolic();
    timer.leave();

    //write_data_text(vm, &step_apd);
    write_data(vm, "vm", &step);
  }
}

void Monodomain::assemble_matrices()
{
  bool test = false;
  const int n = mesh->get_nen();
  const int nz = parameters["maxnz"];
  const double um2_to_cm2  = 1.0e-8;
  const double theta = parameters["theta_method"];
  const double surf_to_vol = parameters["surface_to_volume"];
  const double kappa = (surf_to_vol/timestep) * um2_to_cm2;
  const double dtkappa = 1.0/kappa;

  if(Ai.is_null())
  {
    cout << "Creating FEM matrices" << endl;
    Ai.create(ndofs, ndofs, nz);
    Mi.create(ndofs, ndofs, nz);
    v0.create(ndofs);
    v1.create(ndofs);

    //vm = arma::vec(ndofs);
    //vm.fill(0);
    //v1.create(ndofs);
    //v1.place_array(vm.memptr());

    f.create(ndofs);
    test = true;
  }
  else
  {
    Ai = 0.0;
    Mi = 0.0;
  }

  arma::mat elmatA(n,n), elmatM(n,n);
  arma::mat elmat_k(n,n), elmat_m(n,n);
  std::vector<int> dnums;

  FiniteElement & fe = fespace.createFE(0);

  for(int e=0; e<mesh->get_n_elements(); e++)
  {
    calc_elmat_stiff (e, fe, elmat_k);
    calc_elmat_mass  (e, fe, elmat_m);
    fespace.get_element_dofs(e,dnums);
    
    // Slower assembly, easier to understand
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
      {
        double aval = elmat_m(i,j) + theta * dtkappa * elmat_k(i,j);
        double mval = elmat_m(i,j) - theta * dtkappa * elmat_k(i,j);
        int ii = dnums[i];
        int jj = dnums[j];
        Ai.add(ii, jj, aval);
        Mi.add(ii, jj, mval);
      }
    
    //
    // Fast assembly inserting several values at once
    //
    //int * pidx = &dnums[0];
    //elmatA = elmat_m + (theta*dtkappa)*elmat_k;
    //elmatM = elmat_m - (theta*dtkappa)*elmat_k;
    //elmatA = elmatA.t();
    //elmatM = elmatM.t();
    //Ai.add(n, n, pidx, pidx, elmatA.memptr());
    //Mi.add(n, n, pidx, pidx, elmatM.memptr());
  }

  delete & fe;

  Ai.assemble();
  Mi.assemble();

  if(test)
    solver.init();
  
  string ksptype = "cg";
  string pctype = "jacobi";
  solver.set_solver_type(ksptype);
  solver.set_preconditioner(pctype);
}

void Monodomain::calc_cond_tensor(const int index, const int ndim,
                                  arma::mat & sigma)
{
  const double sigma_l = parameters["sigma_l"];
  const double sigma_t = parameters["sigma_t"];
  const double sigma_n = parameters["sigma_n"];
  const arma::vec3 f = get_fiber(index);
  const arma::vec3 s = get_trans(index);
  const arma::vec3 n = get_normal(index);
  arma::mat tmp;
  arma::mat33 I = arma::eye(3,3);
  sigma = I;

  if(mesh->get_prop_type() == ISOTROPIC)
  {
    int nr = sigma.n_rows;
    for(int i=0; i<nr; i++)
      sigma(i,i) = sigma_l;
  }
  else if (mesh->get_prop_type() == TRANSVERSELY_ISOTROPIC)
  {
    tmp = sigma_t*I + (sigma_l-sigma_t) * (f * f.t());
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        sigma(i,j) = tmp(i,j);
  }
  else if (mesh->get_prop_type() == ORTHOTROPIC)
  {
    for(int k=0; k<ndim; k++)
      for(int i=0; i<ndim; i++)
        sigma(k,i) = sigma_l*f[k]*f[i] + sigma_t*s[k]*s[i] + sigma_n*n[k]*n[i];
  }
}

void Monodomain::calc_elmat_stiff(const int eindex,
				  const FiniteElement & fe,
                                  arma::mat & elmat)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  double detJxW;
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat sigma(ndim,ndim);

  Quadrature * qd = Quadrature::create(fe.order(),fe.type());
  Mapping em = fe.get_mapping(eindex);

  elmat.zeros();
  calc_cond_tensor(eindex,ndim,sigma);

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();

    gradn = dshape * jacinv;
    elmat += (((gradn * sigma) * gradn.t()) * detJxW);
  }

  delete qd;
}

void Monodomain::calc_elmat_mass(const int eindex,
				 const FiniteElement & fe, 
                                 arma::mat & elmat)
{
  double detJxW;
  arma::vec shape;
  arma::mat dshape;

  Quadrature * qd = Quadrature::create(2*fe.order(),fe.type());
  Mapping em = fe.get_mapping(eindex);

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q), shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    em.calc_jacobian(dshape);
    detJxW = qd->get_weight(q) * em.get_det_jacobian();

    elmat += detJxW * (shape * shape.t()) ;
  }

  delete qd;
}

void Monodomain::init()
{
  tip = TimeParameters(timestep, totaltime, printrate);

  std::string output = mesh_filename + "_output";
  std::string extension = file_extension(mesh_filename);

  if(extension == "xml")
  {
    mesh->read_xml(mesh_filename);
    stimuli.read_xml(stimuli_filename);
    std::size_t pos = mesh_filename.find(".xml");
    std::string mesh_name = mesh_filename.substr(0,pos);
    output = mesh_name + "_output";
  }
  else
  {
    mesh->read(mesh_filename);
    stimuli.read(stimuli_filename);
  }

  // stimulus initialization
  stim_values.resize(mesh->get_n_points());
  stim_values.fill(0.0);
  
  fespace.set_mesh(mesh);
  ndofs = mesh->get_n_points();

  // setup parameters
  if(extension == "xml")
    parameters.parse_xml(mesh_filename, "electrophysiology");
  cout << parameters.str();

  // setup data writer
  int nsteps = tip.get_nsteps()/(printrate/timestep);
  cout << "HDF5 Writer" << endl;
  cout << " HDF5 save rate = " << printrate << endl;
  cout << " Number of steps = " << nsteps << endl;
  cout << " Timestep = " << timestep << endl;
  writer->open(output, nsteps, timestep);

  // setup model and cells
  cellmodel = CellModel::create(cell_name);
  cellmodel->setup(odesolver, timestep, totaltime, 1.0);
  cells = new Cells(ndofs,cellmodel);

  vm.resize(ndofs);

  timer.enter("Assemble");
  assemble_matrices();
  timer.leave();
}

void Monodomain::initial_conditions()
{
  tip.reset();

  cells->init();
  cells->get_var(0,v1);
}

void Monodomain::set_stimulus_value(int index, double val)
{
  stim_apply_nodes = true;
  stim_values(index) = val;
}

void Monodomain::solve()
{
  cout << "\nSimulating" << endl;

  initial_conditions();
    
  // loop in time
  int step=0;
  //int step_apd=0;

  stimuli.check(tip.time(), *mesh, stim_nodes, &stim_val, &stim_apply);
  cells->advance(tip.time(), timestep, stim_val, stim_nodes);
  cells->get_var(0, v0);

  while( !tip.finished() )
  {
    tip.increase_time();
    tip.show_time();
    
    timer.enter("ODEs");
    solve_odes();
    timer.leave();

    timer.enter("Parabolic");
    solve_parabolic();
    timer.leave();
        
    //write_data_text(vm, &step_apd);
    write_data(vm, "vm_vtu/vm", &step);
  }  

  timer.summary();
}

void Monodomain::solve_odes()
{
  stimuli.check(tip.time(), *mesh, stim_nodes, &stim_val, &stim_apply);
  
  if (stim_apply)
  {
    cells->advance(tip.time(), timestep, stim_val, stim_nodes);
    stim_nodes.clear();
  }
  else if(stim_apply_nodes)
  {
    cout << "Aplicando estimulos " << tip.time() << endl;
    cells->advance(tip.time(), timestep, stim_values);
    stim_values.fill(0);
    stim_apply_nodes = false;
  }
  else
  {
    cells->advance(tip.time(), timestep);
  }

  cells->get_var(0, v0);  
  v0.assemble();
}

void Monodomain::solve_parabolic()
{
  const double pcgtol = parameters["pcgtol"];

  // b = M * v0 (sparse matrix vector multiplication)
  // A * v1 = b (solve)

  std::pair<PetscInt,PetscReal> ir;
  Mi.mult(v0, f);

#ifdef ELECTRIC_FIELD
  if(tip.time() > 10 && tip.time() < 15){
    ZeroFunction<double> zerofunc;

    if(mesh->get_n_boundary_elements() > 0)
    {
      FiniteElement & bfe = fespace.create_boundary_FE(0);
      int m = bfe.get_ndof();
      arma::vec belvec(m);
      vector<int> bdnums;

      // Assemble boundary bilinear form to impose boundary conditions
      for(int i=0; i < mesh->get_n_boundary_elements(); i++)
      {
        calc_robin_elvec (i, bfe, zerofunc, zerofunc, belvec);

        fespace.get_boundary_element_dofs(i,bdnums);
        for(int k=0;k<m;k++) {
          f.add(bdnums[k],belvec(k));
        }
      }

      delete &bfe;
    }
  }
#endif

  ir = solver.solve(Ai, v1, f, pcgtol);
  
  if (tip.time2print())
    cout << " num its " << ir.first << " rnorm " << ir.second << endl;
   
  // copy solution from PETSc Vec to my Vector
  v1.get_data(vm.memptr());
  cells->set_var(0, vm);
}

void Monodomain::update_coords(const arma::mat & um)
{
  cout << "Updating coordinates (Monodomain)" << endl;   
  for(uint i=0; i<mesh->get_n_points(); i++)
  {
    arma::vec3 upt = (um.row(i)).t();
    mesh->update_point(i, upt);
  }
}

