#include "bidomain_bath.hpp"

//#define USE_ISOTROPIC

BidomainBath::BidomainBath() 
  : sigma_il(0.000300), sigma_it(0.000100), sigma_in(0.000031525), 
    sigma_el(0.000200), sigma_et(0.000165), sigma_en(0.000135140)
{
  mesh = new Mesh();
  writer = new WriterHDF5(mesh);
}

BidomainBath::~BidomainBath()
{
  delete cells;
  delete cellmodel;
}

void BidomainBath::assemble_matrices()
{
  const int n  = mesh->get_nen();
  const int nz = parameters["maxnz"];
  const double um2_to_cm2  = 1.0e-8;
  const double theta = parameters["theta_method"];
  const double surf_to_vol = parameters["surface_to_volume"];
  const double kappa = (surf_to_vol/timestep) * um2_to_cm2;
  const double dtkappa = 1.0/kappa;
  
  cout << endl << " Assembling PETSc-FEM matrices" << endl;

  // Matrices
  Aie.create(ndof,ndof,nz);
  Ai.create(ndof_tissue,ndof_tissue,nz);
  Mi.create(ndof_tissue,ndof_tissue,nz);
  Ki.create(ndof_tissue,ndof_tissue,nz);

  // RHS vectors
  vm0.create(ndof_tissue);
  vm1.create(ndof_tissue);
  b0.create(ndof_tissue);
  b1.create(ndof_tissue);
  aux_vec.create(ndof_tissue);
  ve.create(ndof);
  be.create(ndof);

  arma::mat elmat_ki(n,n), elmat_kie(n,n), elmat_m(n,n), elmat_b(n,n);
  std::vector<int> dnums;

  FiniteElement & fe = fespace.createFE(0);

  for(int e=0; e<mesh->get_n_elements(); e++)
  {
    calc_elmat_stiff_ie(e, fe, elmat_ki, elmat_kie);
    calc_elmat_mass(e, fe, elmat_m);

    get_dof_nums(e,dnums);

    // tissue
    if (mesh->get_element_index(e) != 0)
    {
      // assembly of elliptic problem matrix
      for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
          Aie.add(dnums[i], dnums[j], elmat_kie(i,j));

      // assembly of parabolic problem matrices
      for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){
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
      this->calc_elmat_bath(e, fe, elmat_b);

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

  cout << " End of PETSc-FEM matrices assembly" << endl;

  par_solver.init();
  ell_solver.init();

  delete &fe;
}

void BidomainBath::calc_cond_tensors(const int index, const int ndim,
                                     arma::mat & sigma_i, arma::mat & sigma_e)
{
  const arma::mat I = arma::eye(3,3);
  arma::mat tmpi, tmpe;
  arma::vec3 l, t, n;

  Element e = mesh->get_element(index);
  l = e.get_fiber();
  t = e.get_trans();
  n = e.get_normal();

  // TODO: treat isotropic case ?  

  if (t[0]<0 && t[1]<0 && t[2]<0 && n[0]<0 && n[1]<0 && n[2]<0)
  {
    // tranversely isotropic case
    arma::vec3 al = l;

    tmpi = sigma_it * I + (sigma_il-sigma_it) * (al * al.t() );
    tmpe = sigma_et * I + (sigma_el-sigma_it) * (al * al.t() );

    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
      {
        sigma_i(i,j) = tmpi(i,j);
        sigma_e(i,j) = tmpe(i,j);
      }
  }
  else
  {
    // orthotropic case
    for(int k=0; k<ndim; k++)
      for(int i=0; i<ndim; i++)
      {
        sigma_i(k,i) = l[k]*sigma_il*l[i] + t[k]*sigma_it*t[i] + n[k]*sigma_in*n[i];
        sigma_e(k,i) = l[k]*sigma_el*l[i] + t[k]*sigma_et*t[i] + n[k]*sigma_en*n[i];
      }
  }

#ifdef USE_ISOTROPIC
  sigma_i.zeros();
  for(int i=0; i<ndim; i++)
    sigma_i(i,i) = 0.0003;

  sigma_e.zeros();
  for(int i=0; i<ndim; i++)
    sigma_e(i,i) = 0.0002;
#endif

}

void BidomainBath::calc_elmat_bath(const int iel, const FiniteElement & fe,
                                   arma::mat & elmat_b)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  double detJxW;

  Quadrature * qd = Quadrature::create(fe.order(),fe.type());
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  //Matrix sigma_i(ndim,ndim);
  //Matrix sigma_e(ndim,ndim);
  Mapping em = fe.get_mapping(iel);

  elmat_b.zeros();

  const double sigmab = 0.0010; // mS/um

  //calc_cond_tensors(iel,ndim,sigma_i,sigma_e);

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();

    gradn = dshape * jacinv;

    elmat_b += (detJxW * sigmab) * gradn * gradn.t();
  }

  delete qd;
}

void BidomainBath::calc_elmat_mass(const int iel, const FiniteElement & fe,
                                   arma::mat & elmat)
{
  double detJxW;
  Mapping em = fe.get_mapping(iel);
  arma::vec shape;
  arma::mat dshape;

  Quadrature * qd = Quadrature::create(2*fe.order(),fe.type());

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_shape(qd->get_point(q), shape);
    fe.calc_deriv_shape(qd->get_point(q), dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();

    elmat += detJxW * shape * shape.t();
  }

  delete qd;
}

void BidomainBath::calc_elmat_stiff_ie(const int iel, const FiniteElement & fe,
                                       arma::mat & elmat_i, arma::mat & elmat_ie)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  double detJxW;

  Quadrature * qd = Quadrature::create(fe.order(),fe.type());
  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat sigma_i(ndim,ndim);
  arma::mat sigma_e(ndim,ndim);

  Mapping em = fe.get_mapping(iel);

  elmat_i.zeros();
  elmat_ie.zeros();

  calc_cond_tensors(iel,ndim,sigma_i,sigma_e);

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();

    gradn = dshape * jacinv;

    elmat_i  += detJxW * (gradn * sigma_i) * gradn.t();
    elmat_ie += detJxW * (gradn * (sigma_i + sigma_e)) * gradn.t();
  }

  delete qd;
}

void BidomainBath::get_dof_nums(int i, std::vector<int> & dnums)
{
  fespace.get_element_dofs(i,dnums);
  for (uint i=0; i<dnums.size(); i++)
    dnums[i] = bidomap[(uint)dnums[i]];
}

void BidomainBath::init()
{
  cout << "BidomainBath solver" << endl;

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

  // renumber DoFs if necessary
  setup_nodes();

  fespace.set_mesh(mesh);

  ndof = mesh->get_n_points();

  // setup model and cells
  cellmodel = CellModel::create(cell_name);
  cellmodel->setup(odesolver, timestep, totaltime, 1.0);
  cells = new Cells(ndof_tissue, cellmodel);
  
  neq = cells->get_ode_size();

  tmp_vm.resize(ndof_tissue);
  tmp_ve.resize(ndof);

  // Allocate memory for state variables array
  //sv = new double[ndof_tissue*neq];

  assemble_matrices();
}

void BidomainBath::initial_conditions()
{
  tip.reset();

  // initial conditions
  cells->init();
  cells->get_var(0,vm1);

  //sv = cells->get_state_vars();
}

void BidomainBath::setup_types(std::string & f)
{        
  if( !file_exists(f) )
    cout << "BidomainBath cells: all cells are of the same type" << endl;
  else
  {
    cout << "BidomainBath cells: configuring cell types" << endl;
   
    int aux, idx, size;
    int * vtypes;
    ifstream in(f.c_str());
    in >> size;
    vtypes = new int[size];
    for(int i=0; i<size; i++)
    {
      in >> aux;
      idx = bidomap[i];
      vtypes[idx] = aux;
    }
    cells->set_cell_types(size, vtypes);       
    
    // write in file 
    ofstream out("bidomain_cell_types.typ");
    out << size << endl;
    for(int i=0; i<size; i++)
    {
      out << vtypes[i] << endl;
    }
    out.close();
    
    delete [] vtypes;
  }
}

void BidomainBath::setup_mapping()
{
   
}

void BidomainBath::setup_nodes()
{
  int index;
  uint k;
  std::vector<int> ptnum;

  // tissue
  k = 0;
  for(int i=0; i<mesh->get_n_elements(); i++)
  {
    index = mesh->get_element(i).get_index();
    ptnum = mesh->get_element(i).get_pt_nums();
    if (index != 0)
    {
      for(uint j=0; j<ptnum.size(); j++)
        tissueset.insert(ptnum[j]);
      k++;
    }
  }
  num_tissue_elems = k;

  // bath
  k = 0;
  for(int i=0; i<mesh->get_n_elements(); i++)
  {
    index = mesh->get_element(i).get_index();
    ptnum = mesh->get_element(i).get_pt_nums();
    if (index == 0)
    {
      for(uint j=0; j<ptnum.size(); j++)
        bathset.insert(ptnum[j]);
      k++;
    }
  }
  num_bath_elems = k;

  // remove tissue nodes present at bathset
  std::unordered_set<uint>::iterator b;
  for(b=bathset.begin(); b!=bathset.end(); ++b)
  {
    if (tissueset.find(*b) != tissueset.end())
      bathset.erase(bathset.find(*b));
  }

  // setup bidomain dofs
  ndof_bath   = bathset.size();
  ndof_tissue = tissueset.size();
  
  if (ndof_bath==0)
  {
    cout << " BidomainBath without bath/torso" << endl;
    
    for(uint k=0; k<ndof_tissue; k++)
      bidomap[k] = k;
  }
  else
  {
    // renumbering of dofs
    cout << " BidomainBath with bath/torso -> Renumber DoFs" << endl;
    k=0;
    for(b=tissueset.begin(); b!=tissueset.end(); ++b)
    {
      bidomap[*b] = k;
      k++;
    }

    for(b=bathset.begin(); b!=bathset.end(); ++b)
    {
      bidomap[*b] = k;
      k++;
    }
  }
  
  // reserve memory in std::vectors
  connec_vm.resize(num_tissue_elems);
  connec.resize(num_tissue_elems + num_bath_elems);

  // now fill in tissue-bath connectivities
  int l=0, lt=0;
  for (int i=0; i<mesh->get_n_elements(); i++)
  {
    Element e = mesh->get_element(i);
    std::vector<int> ptnums = e.get_pt_nums();
    std::vector<uint> ptnums_tmp(ptnums.size());
    std::vector<uint> ptnums_tmp_vm(ptnums.size());
    
    // tissue
    if (e.get_index() != 0)
    {
      connec[l].reserve(ptnums.size());
      connec_vm[lt].reserve(ptnums.size());
      for(uint j=0; j<ptnums.size(); j++)
      {
        uint glob = ptnums[j];
        ptnums_tmp[j] = glob;
        ptnums_tmp_vm[j] = bidomap[glob];
      }
      connec.push_back( ptnums_tmp );
      connec_vm.push_back( ptnums_tmp_vm );
      l++;
      lt++;
    }    
  }

  for (int i=0; i<mesh->get_n_elements(); i++)
  {
    Element e = mesh->get_element(i);
    std::vector<int> ptnums = e.get_pt_nums();
    std::vector<uint> ptnums_tmp(ptnums.size());

    // bath
    if (e.get_index() == 0)
    {
      connec[l].reserve(ptnums.size());
      for(uint j=0; j<ptnums.size(); j++)
      {
        uint glob = ptnums[j];
        ptnums_tmp[j] = glob;
      }
      connec.push_back( ptnums_tmp );
      l++;
    }
  }

  cout << endl;
  cout << " Number of vm (tissue) dofs      " << ndof_tissue << endl;
  cout << " Number of ve (tissue+bath) dofs " << ndof_bath << endl;
  cout << " BidomainBath dofs map size          " << bidomap.size() << endl;
  cout << " Bathset size                    " << bathset.size() << endl;
  cout << " Tissueset size                  " << tissueset.size() << endl;
  cout << " Number of bath elements         " << num_bath_elems << endl;
  cout << " Number of tissue elements       " << num_tissue_elems << endl;
  
  //for(uint k=0; k<bidomap.size(); k++)
  //    cout << "k= " << k << " bido= " << bidomap[k] << endl; 
}

void BidomainBath::solve ()
{
  cout << "\n Simulating\n\n" << flush;

  // setup initial conditions on cells
  initial_conditions();

  // loop in time
  int vm_step=0;
  int vm_apd_step=0;
  int ve_step=0;

  while( !tip.finished() )
  {
    tip.increase_time();
    
    timer.enter("ODEs");
    solve_odes();
    timer.leave();
    
    timer.enter("Parabolic");
    solve_parabolic();
    timer.leave();
    
    timer.enter("Elliptic");
    solve_elliptic();
    timer.leave();
    
    write_data_text(tmp_vm, &vm_apd_step);
    write_data(tmp_vm, "vm_vtu/vm", &vm_step);
    write_data(tmp_ve, "ve_vtu/ve", &ve_step);
    
    // old
    //write_vtk_vm(tmp_vm, "output/vm_vtu/vm_tissue",&vm_step);
    //write_vtk_ve(tmp_ve, "output/ve_vtu/ve", &ve_step);
   
    if (tip.time2print()) cout << endl;
  }
  
  timer.summary();  
}

void BidomainBath::solve_odes()
{
  stimuli.check_bidomain(tip.time(), *mesh, bidomap, si_nodes,
                         &stimval, &intra_apply);

  if (intra_apply)
    cells->advance(tip.time(), timestep, stimval, si_nodes);
  else
    cells->advance(tip.time(), timestep);

  cells->get_var(0,vm0);
  vm0.assemble();
}

void BidomainBath::solve_parabolic()
{
  const double pcgtol = parameters["pcgtol"];
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
    cout << " Time " << tip.time();
    cout << " Parabolic its " << ir.first;
    cout << " rnorm " << scientific << ir.second;
  }

  // Copy solution from PETSc Vec to my Vector
  vm1.get_data(tmp_vm.memptr());

  // Copy data from tmp_vm to state variables array
  cells->set_var(0,tmp_vm);
}

void BidomainBath::solve_elliptic()
{
  const double pcgtol = parameters["pcgtol"];
  // Solve
  // div((sigma_i + sigma_e) grad(Vm)) = - div(sigma_i grad(vm))
  std::pair<PetscInt,PetscReal> ir;

  // NO BATH
  //Ki.mult(vm1, be);
  //be.scale(-1.0);

  // BATH
  Ki.mult(vm1, aux_vec); // temp
  be.set(0.0);           // form be = [be_tmp 0]
  be.copy_values(ndof_tissue, aux_vec);
  be.scale(-1.0);

  // FIX FIX FIX 
  if(ndof_bath != 0) 
    be.set(ndof_bath-1, 0.0);
  else
    be.set(0, 0.0);
    
  ir = ell_solver.solve(Aie, ve, be, pcgtol);

  if (tip.time2print())
    cout << "  Elliptic its " << ir.first << " rnorm " << ir.second;

  // copy solution from PETSc Vec to my Vector
  ve.get_data(tmp_ve.memptr());
}

/*
void BidomainBath::write_vtk_vm(const arma::vec & u, const std::string & s, 
			    int * step)
{
  if (tip.time2print())
  {
    char buf[256];
    sprintf(buf,"%s_%03d.vtu",s.c_str(),*step);
    *step = *step + 1;
    std::string file(buf);

    ofstream vtu;
    vtu.open(file.c_str());
    vtkout.set_points(tissueset);
    vtkout.set_cells(connec_vm);
    vtkout.write_point_data(u, "vm");
    vtkout.write_vtu(file.c_str());
    vtu.close();
  }
}

void BidomainBath::write_vtk_ve(const arma::vec & u, const std::string & s, 
			    int * step)
{
  if (tip.time2print())
  {
    // use bidomap to create the correct Ve vector to output
    arma::vec ve;
    ve.resize(u.size());
    std::map<uint, uint>::iterator p;
    for(uint k=0; k<bidomap.size(); k++)
      ve(k) = u(bidomap[k]);

    // prepare to save
    char buf[256];
    sprintf(buf,"%s_%03d.vtu",s.c_str(),*step);
    *step = *step + 1;
    std::string file(buf);

    ofstream vtu;
    vtu.open(file.c_str());
    vtkout.set_points(mesh->get_points());
    vtkout.set_cells(connec);
    
    //vtkout.write_point_data(u, "vm"); // there was a bug here
    vtkout.write_point_data(ve, "ve");
    
    vtkout.write_vtu(file.c_str());
    vtu.close();
  }
}

*/