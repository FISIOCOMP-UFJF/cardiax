#include "monodomain_purkinje.hpp"
#include "mesh/purkinje_mesh.hpp"

MonodomainPurkinje::MonodomainPurkinje()
  : CardiacProblem(), sigma_l(2e-3), radius(10.0)
{
  // Purkinje fiber conduction velocity: 2.4 m/s
  // Tests for different conductivities
  // sigma_l=3e-4 mS/um and dx=100um ---> CV: 0.85 m/s
  // sigma_l=1e-1 mS/um and dx=100um ---> CV: 12.8 m/s
  // sigma_l=2e-3 mS/um and dx=100um ---> CV: 2.50 m/s

  // create the Purkinje network mesh
  mesh = new PurkinjeMesh();

  // create writers
  writer = new WriterHDF5(mesh);
}

MonodomainPurkinje::~MonodomainPurkinje()
{
  delete cells;
  delete cellmodel;
}

void MonodomainPurkinje::advance()
{
  static int step = 0;

  if( !tip.finished() )
  {
    tip.increase_time();

    timer.enter("ODEs");
    solve_odes();
    timer.leave();

    timer.enter("Parabolic");
    solve_parabolic();
    timer.leave();

    write_data(vm, "vm", &step);
  }
}

void MonodomainPurkinje::assemble_matrices()
{
  const int nz = parameters["maxnz"];
  const double um2_to_cm2  = 1.0e-8;
  const double theta = parameters["theta_method"];
  const double surf_to_vol = parameters["surface_to_volume"];
  const double kappa = (surf_to_vol/timestep) * um2_to_cm2;
  const double dtkappa = 1.0/kappa;

	cout << "Creating PETSc-FEM matrices" << endl;
	
	Ai.create(ndofs, ndofs, nz);
	Mi.create(ndofs, ndofs, nz);
	v0.create(ndofs);
	v1.create(ndofs);
	f.create(ndofs);

  // Alternative to use pointer to memory instead of copying all the time
  //vm = arma::vec(ndofs);
  //vm.fill(0);
  //v1.create(ndofs);
  //v1.place_array(vm.memptr());

  FiniteElement & fe = pk_fespace.createFE(0);
	const int n = fe.get_ndof();

	arma::mat elmatA(n,n), elmatM(n,n);
  arma::mat elmat_k(n,n), elmat_m(n,n);
	std::vector<int> dnums;

  // Computing bifurcations in order to impose
  // conservation of charge (Kirchhoff's law)
  int ndofs_u = pk_fespace.get_num_dofs()/2;

  PurkinjeMesh::BifucartionIterator it;
  PurkinjeMesh * pk_mesh = static_cast<PurkinjeMesh*>(mesh);
  pk_mesh->set_root(0);
  pk_mesh->setup();

  for(it = pk_mesh->begin(); it != pk_mesh->end(); ++it)
  {
    const int bif_node = (*it).first;
    bif_dofs.push_back(bif_node + ndofs_u);
  }

  // Assemble FEM matrices of the Purkinje network
  for(int e=0; e<pk_mesh->get_n_elements(); e++)
  {
    calc_elmat_stiff (e, fe, elmat_k);
    calc_elmat_mass  (e, fe, elmat_m);

    pk_fespace.get_element_dofs(e,dnums);
    
    // Slower assembly, easier to understand
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
      {
        double aval = elmat_m(i,j) + theta * dtkappa * elmat_k(i,j);
        double mval = elmat_m(i,j) - theta * dtkappa * elmat_k(i,j);
        int ii = dnums[i];
        int jj = dnums[j];

        // we only assemble the components which are not bifucartions
        if( ! (std::find(bif_dofs.begin(), bif_dofs.end(), ii) != bif_dofs.end()) )
        {
          Ai.add(ii, jj, aval);
          Mi.add(ii, jj, mval);
        }
      }
  }

  // Assembly the equations of each bifurcation
  for(it = pk_mesh->begin(); it != pk_mesh->end(); ++it)
  {
    std::pair< int, std::vector<int> > bif = *it;
    std::vector<int> neighbors = bif.second;
    const int bif_node = bif.first + ndofs_u;

    Ai.add(bif_node, bif_node, 1.0);
    for(std::vector<int>::iterator pn = neighbors.begin(); pn != neighbors.end(); ++pn)
    {
      const int j = *pn + ndofs_u;
      Ai.add(bif_node, j, -1.0);
    }
  }

  delete & fe;

  Ai.assemble();
  Mi.assemble();

  string ksptype = "gmres";
  string pctype = "asm";
  solver.init();
  solver.set_solver_type(ksptype);
  solver.set_preconditioner(pctype);
}

void MonodomainPurkinje::calc_elmat_stiff(const int eindex,
                                          const FiniteElement & fe,
																					arma::mat & elmat)
{
  int ndim = fe.get_ndim();
  int ndof = fe.ndofs();
  double detJxW;
  const double cond_sigma = sigma_l;

  arma::mat dshape;
  arma::mat gradn(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat ji(1,3);

  Quadrature * qd = Quadrature::create(fe.order(),fe.type());
  Mapping em = fe.get_mapping(eindex);

  elmat.zeros();

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    fe.calc_deriv_shape(qd->get_point(q),dshape);
    em.calc_segment_jacobian(dshape);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();

		ji(0,0) = jacinv(0,0);
		ji(0,1) = jacinv(1,0);
		ji(0,2) = jacinv(2,0);
		gradn = (dshape * ji);
    elmat += (((gradn * cond_sigma) * gradn.t()) * detJxW);
  }

  delete qd;
}

void MonodomainPurkinje::calc_elmat_mass(const int eindex,
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
    em.calc_segment_jacobian(dshape);		
    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    elmat += detJxW * (shape * shape.t()) ;
  }

  delete qd;
}

void MonodomainPurkinje::init()
{
  tip = TimeParameters(timestep, totaltime, printrate);

  mesh->read(mesh_filename);
  stimuli.read(stimuli_filename);

  pk_fespace.set_mesh(mesh);
	pk_fespace.config();

	ndofs = pk_fespace.get_num_dofs();

  // setup data writer
  std::string output = mesh_filename + "_output";

  int nsteps = tip.get_nsteps()/(printrate/timestep);
  cout << "HDF5 Writer" << endl;
  cout << " HDF5 save rate = " << printrate << endl;
  cout << " Number of steps = " << nsteps << endl;
  cout << " Timestep = " << timestep << endl;

  writer->open(output, nsteps, timestep, true);
  writer->add_ve();

  // setup model and cells
  cellmodel = CellModel::create(cell_name);
  cellmodel->setup(odesolver, timestep, totaltime, 1.0);

  // hermite FE
  int ncells = ndofs/2;
  cells = new Cells(ncells,cellmodel);

  vm.resize(ndofs);

  timer.enter("Assemble");
  assemble_matrices();
  timer.leave();
}

void MonodomainPurkinje::initial_conditions()
{
  tip.reset();

  cells->init();
  cells->get_var(0,v1);
}

void MonodomainPurkinje::solve()
{
  cout << "\nSimulating" << endl << endl << flush;

  initial_conditions();
    
  // loop in time
  int step=0;
  //int step_apd=0;
  int step_dv=0;

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
    write_data(vm, "vm_vtu/dv", &step_dv);
  }  

  timer.summary();
}

void MonodomainPurkinje::solve_odes()
{
  bool apply;
  double value;

  stimuli.check(tip.time(), *mesh, stim_nodes, &value, &apply);
  
  if (apply)
  {
    cells->advance(tip.time(), timestep, value, stim_nodes);
    stim_nodes.clear();
  }
  else
  {
    cells->advance(tip.time(), timestep);
  }

  cells->get_var(0, v0);  
  v0.assemble();
}

void MonodomainPurkinje::solve_parabolic()
{
  const double pcgtol = parameters["pcgtol"];

  // b = M * v0 (sparse matrix vector multiplication)
  // A * v1 = b (solve)
  std::pair<PetscInt,PetscReal> ir;
  Mi.mult(v0, f);

  ir = solver.solve(Ai, v1, f, pcgtol);

  if (tip.time2print())
    cout << " num its " << ir.first << " rnorm " << ir.second << endl;

  /*
  // Post-processing step for 1D Cubic Hermite finite elements
  CubicHermiteFESpace * fec = dynamic_cast<CubicHermiteFESpace*>(&pk_fespace);
  if (fec != NULL)
  {
    cout << "Post-processing derivative...\n";
    arma::vec deriv = solution(arma::span(ndofs/2, ndofs-1));
    FETools::cubic_hermite_deriv(pk_fespace,msh,deriv);
    cout << deriv << endl;
  }
  */

  // TODO: need to improve these copies!!
  // copy solution from PETSc Vec to my Vector
  v1.get_data(vm.memptr());
  cells->set_var(0, vm);
}


