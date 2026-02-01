#include "nonlinear_elasticity.hpp"
#include "util/pugixml.hpp"

//#define DEBUG_PRESSURE
#define INCs_TA 0

NonlinearElasticity::NonlinearElasticity()
  : timer(), lc(), output_step(false), material(0),
    log(), writer(&msh), first_step(true)
{
  // define default parameters
  parameters.add("tol_energy",1.0e-8);
  parameters.add("tol_residual",1.0e-6);
  parameters.add("tol_displacement",1.0e-16);

  //parameters.add("tol_energy",1.0e-6);
  //parameters.add("tol_residual",1.0e-5);
  //parameters.add("tol_displacement",1.0e-16);

  parameters.add("energy_norm0",0.0);
  parameters.add("residual_norm0",0.0);

  // create log file using the format:
  // its1 its2 ... itsN volume
  log.create("log_elasticity.dat");
}

NonlinearElasticity::~NonlinearElasticity()
{
  if (material != NULL)
    delete material;

  for(uint i=0; i<vecF.size(); i++)
    delete vecF[i];
}

void NonlinearElasticity::apply_boundary(petsc::Matrix & Ks)
{
  int nd = msh.get_n_dim();
  int np = msh.get_n_points();
  int nod, idx, dir;
  double val, dlamb = lc.load_step();
  NodalData inf;

  // prescribe dirichlet boundary conditions
  if (fixed_nodes_map.size())
  {
    std::map<int,double> boundary_values;
    std::map<int,NodalData>::iterator it;
    for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
    {
      nod = it->first;
      inf = it->second;
      dir = inf.first;
      val = dlamb * inf.second;
      //idx = (dir * np) + nod;
      idx = (nod*nd) + dir;
      assert(dir+1 <= nd);
      boundary_values.insert( std::pair<int,double>(idx,val) );
    }
    FETools::apply_boundary_values(boundary_values, Ks);
  }
}

void NonlinearElasticity::assemble_traction()
{
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe != NULL && neumann_map.size() > 0)
  {
    cout << "Computing traction boundary conditions" << endl;
    int nu  = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    arma::mat belmat(nu,nu);
    arma::vec belvec(nu);
    vector<int> bdof;
    for(int i=0; i<nb; i++)
    {
      calc_neumann_elvec(i, bfe, belvec);
      fespace.get_boundary_element_dofs_u(i,bdof);

      for(int k=0; k<nu; k++)
        fext(bdof[k]) += belvec(k);
    }
  }

  delete bfe;
}

void NonlinearElasticity::assemble_pressure()
{
  MixedFiniteElement * bfe = fespace.create_boundary_FE();

  if (bfe != NULL && pressure_map.size() > 0)
  {
    cout << " Assembling pressure component of stiffness matrix" << endl;
    int nu  = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    double xlamb = lc.load();
    arma::mat belmat(nu,nu);
    arma::vec belvec(nu);
    vector<int> bdof;

    for(int i=0; i<nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i,bdof);
      calc_pforce_kpress(i,bfe,bdof,belvec,belmat);

      // assembles the pressure component of the stiffness matrix
      for(int j=0; j<nu; j++)
        for(int k=0; k<nu; k++)
          K.add(bdof[j], bdof[k], belmat(j,k));

      // assembles the nodal forces due to normal pressure
      for(int k=0; k<nu; k++)
      {

        if(ldgof[bdof[k]])
        {
          r.add(bdof[k], xlamb * belvec(k));
          //fext(bdof[k]) += belvec(k);
          tload(bdof[k]) += belvec(k);
        }
        else
          react.add(bdof[k], -xlamb * belvec(k));
      }
    }
  }

  K.assemble();
  r.assemble();
  apply_boundary(K);


  delete bfe;
}

double NonlinearElasticity::total_volume_cavity()
{
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  double total_endo_volume=0;

  if (bfe != NULL && pressure_map.size() > 0)
  {
    int nb = msh.get_n_boundary_elements();
    vector<int> bdof;
    for(int i=0; i<nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i, bdof);
      total_endo_volume += calc_cavity_volume(i, bfe, bdof);
    }
  }
  delete bfe;

  return total_endo_volume*1e6;
}


void NonlinearElasticity::body_forces()
{
  if(bforce.size() != 0)
  {
    int ndim  = msh.get_n_dim();
    int ndof  = msh.get_nen() * ndim;
    int nelem = msh.get_n_elements();
    int nubf  = ndof/ndim;
    arma::vec elvec(ndof);
    std::vector<int> dnums;

    MixedFiniteElement * fe = fespace.createFE();
    Quadrature * qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //Quadrature * qd = Quadrature::create(2, fe->get_type());
    cout << " Computing body forces (or forcing term)" << endl;

    for(int i=0; i<nelem; i++)
    {
      double detJxW;
      arma::vec shape;
      arma::mat dshape;
      arma::vec3 pt, qpt;

      Mapping em = fe->get_mapping(i);

      elvec.zeros();
      fespace.get_element_dofs_u(i,dnums);

      for(int q=0; q<qd->get_num_ipoints(); q++)
      {
        qpt = qd->get_point(q);
        fe->calc_shape_u(qpt, shape);
        fe->calc_deriv_shape_u(qpt, dshape);
        em.calc_jacobian(dshape, nubf);
        detJxW = qd->get_weight(q) * em.get_det_jacobian();
        pt = em.map_point(qpt, shape);

        // compute forcing vector
        //arma::vec f = forcing_term(pt, *material);
        arma::vec3 f;
        for(int i=0; i<3; i++) f(i) = bforce[i];

        for(int i=0; i<ndim; i++)
          for(int j=0; j<nubf; j++)
            elvec(j+i*nubf) += f(i) * shape(j+i*nubf) * detJxW;
      }

      // assemble the body force into the external load vector
      //   if dof is active, adds the forces due
      //   to gravity/body force to the global vector fext
      for(int k=0; k<ndof; k++)
        if(ldgof[dnums[k]])
          fext(dnums[k]) += elvec(k);
    }
    delete qd;
    delete fe;
  }
}

double NonlinearElasticity::calc_energy()
{
  const int ndofs = fespace.get_ndofs();
  double energy_norm = 0.0;
  arma::vec resid(ndofs);
  arma::vec displ(ndofs);
  r.get_data(resid.memptr());
  u.get_data(displ.memptr());

  for(int i=0; i<ndofs; i++)
    if(ldgof[i])
      energy_norm += resid(i) * displ(i);

  return abs(energy_norm);
}

double NonlinearElasticity::calc_volume(bool update)
{
  int ndof, nnode;
  double totv = 0.0;
  MixedFiniteElement * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());
  //Quadrature * qd = Quadrature::create(2, fe->get_type());
  ndof  = fe->get_ndofs_u();
  nnode = fe->get_nnode();

  for(int i=0; i<msh.get_n_elements(); i++)
  {
    double detJxW, evol = 0.0;
    arma::mat dshape;
    arma::vec3 qpt;
    std::vector<arma::vec3> xe(ndof);
    get_element_x (i, xe);

    Mapping em = fe->get_mapping(i, xe);

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      qpt = qd->get_point(q);
      fe->calc_deriv_shape_u(qpt,dshape);
      em.calc_jacobian(dshape, nnode);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();
      evol = evol + detJxW;
    }

    totv = totv + evol;

    // initial element volume
    if(update) vol0[i] = evol;
  }

  delete qd;
  delete fe;

  return totv;
}

void NonlinearElasticity::calc_neumann_elvec(const int eindex, const MxFE * fe,
                                             arma::vec & elvec)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  int index;
  double detJxW;
  arma::vec shape;
  arma::mat dshape;
  arma::vec3 f;
  std::map<int,arma::vec3>::iterator it;

  SurfaceMapping sm = fe->get_boundary_mapping(eindex);
  index = sm.get_index();
  it = neumann_map.find(index);

  f.zeros();
  elvec.zeros();

  // get traction vector on surface with number=index
  if(neumann_map.find(index) != neumann_map.end())
  {
    Quadrature * qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //Quadrature * qd = Quadrature::create(2, fe->get_type());

    f = it->second;

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, nubf);

      detJxW = qd->get_weight(q) * sm.get_det_jacobian();

      arma::vec3 pt = sm.map_point(qd->get_point(q)[0]);

      for(int i=0; i<ndim; i++)
        for(int j=0; j<nubf; j++)
	        elvec(j+i*nubf) += f(i) * shape(j+i*nubf) * detJxW;
    }
    delete qd;
  }
}

void NonlinearElasticity::elem_pforce(const int elem_id, const MxFE * fe,
                                      const std::vector<int> & bdof,
                                      arma::vec & belvec)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW, press;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();
  it = pressure_map.find(index);

  if(pressure_map.find(index) != pressure_map.end())
  {
    belvec.zeros();
    press = it->second;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);
      // trick
      detJxW = qd->get_weight(q);

      // evaluates dx/dxi, dx/deta
      // for 2D set dx/deta = (0,0,-1)
      dxis.zeros();
      if(ndim==2) dxis(2,1) = -1.0;
      for(int id=0; id<ndim; id++)
        for(int jd=0; jd<ndim-1; jd++)
          for(int in=0; in<neln; in++)
          {
            //int ip = bdof[in];
	    int ip = bdof[in]/ndim;
            dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
          }

      //
      // computes pressure force (nodal forces - subroutine pforce)
      //
      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
      
      //std::cout << xnorm <<std::endl;

      for(int in=0; in<neln; in++)
      {
        double vl = press * shape(in) * detJxW;
        // assemble boundary force to elemental vector (in + id*neln)
        // belvec data layout
        // [ bloco_x, bloco_y, bloco_z]
        // [n1_x, n2_x, ..., nN_x, n1_y, ..., nN_y, n1_z, ..., nN_z]
        for(int id=0; id<ndim; id++)
          belvec(in + id*neln) += xnorm(id) * vl;
      }
    }
    // end of integration loop
    delete qd;
  }
  else
  {
    it = spring_map.find(index);
    if(it != spring_map.end())
    {
      belvec.zeros();
      press = it->second;

      // create quadrature rule (need to check if order is ok)
      Quadrature * qd;
      qd = Quadrature::create(fe->get_order_u(), fe->get_type());
      //qd = Quadrature::create(2, fe->get_type());

      // quadrature loop
      for(int q=0; q<qd->get_num_ipoints(); q++)
      {
        fe->calc_shape_u(qd->get_point(q),shape);
        fe->calc_deriv_shape_u(qd->get_point(q),dshape);
        sm.calc_jacobian(dshape, neln);
        // trick
        detJxW = qd->get_weight(q);

        // evaluates dx/dxi, dx/deta
        // for 2D set dx/deta = (0,0,-1)
        dxis.zeros();
        if(ndim==2) dxis(2,1) = -1.0;
        for(int id=0; id<ndim; id++)
          for(int jd=0; jd<ndim-1; jd++)
            for(int in=0; in<neln; in++)
            {
              //int ip = bdof[in];
              int ip = bdof[in]/ndim;
              dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
            }

        //
        // computes pressure force (nodal forces - subroutine pforce)
        //
        arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
        for(int in=0; in<neln; in++)
        {
          double vl = press * shape(in) * detJxW;
          double vl2 = 0;
          for(int jn=0; jn<neln; jn++)
            vl += shape(jn)*xnorm(jn);
          // assemble boundary force to elemental vector (in + id*neln)
          // belvec data layout
          // [ bloco_x, bloco_y, bloco_z]
          // [n1_x, n2_x, ..., nN_x, n1_y, ..., nN_y, n1_z, ..., nN_z]
          
          // Tentativa de fazer J*(uh*n)*n*vh
          for(int id=0; id<ndim; id++)
            belvec(in + id*neln) += xnorm(id) * vl*vl2;
        }
      }
      // end of integration loop
      delete qd;
    }
    else
    {
      error("Boundary index for pressure load not found");  
    }
  }
}

void NonlinearElasticity::elem_kpress(const int elem_id, const MxFE * fe,
                                      const std::vector<int> & bdof,
                                      arma::mat & belmat)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW, press;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();
  //std::cout << "\n" << index << " ";
  it = pressure_map.find(index);

  if(it != pressure_map.end())
  {
    belmat.zeros();
    press = it->second;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);
      // trick
      detJxW = qd->get_weight(q);

      // evaluates dx/dxi, dx/deta
      // for 2D set dx/deta = (0,0,-1)
      dxis.zeros();
      if(ndim==2) dxis(2,1) = -1.0;
      for(int id=0; id<ndim; id++)
        for(int jd=0; jd<ndim-1; jd++)
          for(int in=0; in<neln; in++)
          {
            //int ip = bdof[in];
	          int ip = bdof[in]/ndim;
            dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
          }

      //
      // computes pressure component for stiffness matrix (subroutine kpress)
      //
      double sum, val, val2, val3;
      double apress = press * lc.load_step();
      val = -(apress * detJxW)/2.0;

      for(int in=0; in<neln; in++)
        for(int jn=0; jn<neln; jn++)
        {
          val2 = val * ( shape(in)*dshape(jn,0)
                       - shape(jn)*dshape(in,0) );
          val3 = 0.0;
          if(ndim==3) val3 = val * ( shape(in)*dshape(jn,1)
                                   - shape(jn)*dshape(in,1) );
          for(int id=0; id<ndim; id++)
            for(int jd=0; jd<ndim; jd++)
            {
              sum = 0.0;
              for(int kd=0; kd<3; kd++)
                sum += levi(id,jd,kd) * (val2*dxis(kd,1) - val3*dxis(kd,0));

              // write in the elemental pressure stiffness matrix (in + id*neln)
              // belmat data layout
              // [ bloco_x, bloco_y, bloco_z] , [..]
              // [n1_x, n2_x, ..., nN_x, n1_y, ..., nN_y, n1_z, ..., nN_z], [..]
              belmat(in + id*neln , jn + jd*neln) += sum;
            }
        }
    }
    
    // end of integration loop
    delete qd;
  }
  else
  {
    it = spring_map.find(index);
    
    if(it != spring_map.end())
    {
      belmat.zeros();
      press = it->second;

      // create quadrature rule (need to check if order is ok)
      Quadrature * qd;
      qd = Quadrature::create(fe->get_order_u(), fe->get_type());
      //qd = Quadrature::create(2, fe->get_type());

      // quadrature loop
      for(int q=0; q<qd->get_num_ipoints(); q++)
      {
        fe->calc_shape_u(qd->get_point(q),shape);
        fe->calc_deriv_shape_u(qd->get_point(q),dshape);
        sm.calc_jacobian(dshape, neln);
        // trick
        detJxW = qd->get_weight(q);
        //
        // computes spring component for stiffness matrix (subroutine kpress)
        //
        // evaluates dx/dxi, dx/deta
        // for 2D set dx/deta = (0,0,-1)
        dxis.zeros();
        if(ndim==2) dxis(2,1) = -1.0;
        for(int id=0; id<ndim; id++)
          for(int jd=0; jd<ndim-1; jd++)
            for(int in=0; in<neln; in++)
            {
              //int ip = bdof[in];
              int ip = bdof[in]/ndim;
              dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
            }
          arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));

          //xnorm = xnorm/arma::norm(xnorm,2);

          arma::mat33 TensorProd = xnorm*xnorm.t();

          //std::cout << xnorm <<  "\n" << std::endl;
          //std::cout << TensorProd << std::endl;
        double sum;
        //double apress = press * lc.load_step();

        for(int in=0; in<neln; in++)
          for(int jn=0; jn<neln; jn++)
          {
            for(int id=0; id<ndim; id++)
              for(int jd=0; jd<ndim; jd++)
              {
                //std::cout << in << " " << jn << " " << id << " " << jd << " " << ndim << " " << neln <<std::endl;
                sum = -press*shape(jn)*TensorProd(id,jd)*shape(in)*detJxW;
                
                // write in the elemental pressure stiffness matrix (in + id*neln)
                // belmat data layout
                // [ bloco_x, bloco_y, bloco_z] , [..]
                // [n1_x, n2_x, ..., nN_x, n1_y, ..., nN_y, n1_z, ..., nN_z], [..]
                belmat(in + id*neln , jn + jd*neln) += sum;
              }
          }
      }
      //std::cout << "belmat: \n" << belmat <<std::endl;
      // end of integration loop
      delete qd;
    }
    else
    {
      error("Boundary index for pressure load not found");
    }
  }
}

void NonlinearElasticity::calc_pforce_kpress(const int elem_id, const MxFE * fe,
                                             const std::vector<int> & bdof,
                                             arma::vec & belvec,
                                             arma::mat & belmat)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW, press;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // third order alternating tensor
  arma::cube eps = arma::zeros(3,3,3);
  eps(0,1,2) = 1; eps(1,2,0) = 1; eps(2,0,1) = 1;
  eps(0,2,1) = -1; eps(1,0,2) = -1; eps(2,1,0) = -1;

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();
  it = pressure_map.find(index);

  if(pressure_map.find(index) != pressure_map.end())
  {
    belvec.zeros();
    belmat.zeros();
    press = it->second;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    //qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);

      // Ajuste !
      detJxW = qd->get_weight(q);

      //
      // evaluates dx/dxi, dx/deta
      //
      dxis.zeros();

      // for 2D set dx/deta = (0,0,-1)
      if(ndim==2) dxis(2,1) = -1.0;

      for(int id=0; id<ndim; id++)
	    for(int jd=0; jd<ndim-1; jd++)
	      for(int in=0; in<neln; in++)
	      {
	        //int ip = bdof[in];
		int ip = bdof[in]/ndim;
	        dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
	      }

      //
      // computes pressure force (nodal forces)
      // subroutine pforce
      //
      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));
      for(int in=0; in<neln; in++)
      {
	    double vl = press * shape(in) * detJxW;
	    for(int id=0; id<ndim; id++)
          belvec(in + id*neln) += xnorm(id) * vl;
      }

      //
      // computes the pressure component of the stiffness matrix
      // subroutine kpress
      //
      double sum, val, val2, val3;
      double apress = press * lc.load();

      val = -(apress * detJxW)/2.0;

      for(int in=0; in<neln; in++)
	    for(int jn=0; jn<neln; jn++)
	    {
	      val2 = val * (shape(in)*dshape(jn,0) - shape(jn)*dshape(in,0));
	      val3 = 0.0;

	      if(ndim==3)
	        val3 = val * (shape(in)*dshape(jn,1) - shape(jn)*dshape(in,1));

	      for(int id=0; id<ndim; id++)
	        for(int jd=0; jd<ndim; jd++)
	        {
	          sum = 0.0;
	          for(int kd=0; kd<3; kd++)
		        sum = sum + eps(id,jd,kd) * (val2 * dxis(kd,1) - val3 * dxis(kd,0));

              // write in the elemental pressure stiffness matrix
	          belmat(in + id*neln, jn + jd*neln) += sum;
	        }
	    }
    }
    // end of integration loop
    std::cout << "belmat: \n" << belmat <<std::endl;
    delete qd;
  }

}

double NonlinearElasticity::calc_cavity_volume(const int elem_id, const MxFE * fe,
					                                     const std::vector<int> & bdof)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int neln = ndof/ndim;
  int index;
  double detJxW;
  arma::vec shape;
  arma::mat dshape, dxis(3,2);

  // get current coordinates of the element
  std::vector<arma::vec3> xe(neln);
  get_boundary_element_x (elem_id, xe);

  // boundary element mapping
  SurfaceMapping sm = fe->get_boundary_mapping(elem_id, xe);

  std::map<int,double>::iterator it;
  index = sm.get_index();

  //Return 0 to markers different from LV (30)
  if(index!=30)
	  return 0;

  double endo_volume;


    endo_volume = 0;

    // create quadrature rule (need to check if order is ok)
    Quadrature * qd;
    //qd = Quadrature::create(fe->get_order_u(), fe->get_type());
    qd = Quadrature::create(2, fe->get_type());

    // quadrature loop
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      fe->calc_shape_u(qd->get_point(q),shape);
      fe->calc_deriv_shape_u(qd->get_point(q),dshape);
      sm.calc_jacobian(dshape, neln);
      detJxW = qd->get_weight(q);

      //
      // evaluates dx/dxi, dx/deta
      //
      dxis.zeros(); 

      // for 2D set dx/deta = (0,0,-1)
      if(ndim==2) dxis(2,1) = -1.0;

      for(int id=0; id<ndim; id++)
        for(int jd=0; jd<ndim-1; jd++)
          for(int in=0; in<neln; in++)
	        {
	          //int ip = bdof[in];
		        int ip = bdof[in]/ndim;
	          dxis(id,jd) = dxis(id,jd) + x[ip][id] * dshape(in,jd);
	        }

      arma::vec3 xnorm = arma::cross(dxis.col(0), dxis.col(1));

      for(int in=0; in<neln; in++)
      {
	      arma::vec3 vl = xe[in] * detJxW * shape(in);

        for(int id=0; id<ndim; id++)
            endo_volume += xnorm(id) * vl[id]/3.;
        //int id = 0;
        //endo_volume += xnorm(id) * vl[id];
      }

    }
    // end of integration loop
    delete qd;

  return endo_volume;
}

void NonlinearElasticity::config(const string & mshfile, const string & parfile)
{
  cout << "Setup of non-linear elasticity problem" << endl;
  int num_increments;
  std::string mtype, etype;
  std::string extension = file_extension(parfile);
  std::vector<double> matprop;

  if(extension == "par")
  {
    cout << "Reading parameters from .par file" << endl;

    InputFile ifile(parfile);
    ifile.read("num_increments", num_increments);
    ifile.read("problem_type", etype);
    ifile.read("material_type", mtype);
    ifile.read_array("material_coef", matprop);
    ifile.read_array("body_force", bforce);
    ifile.read_section("neumann", neumann_map);
    ifile.read_section("dirichlet", dirichlet_map);
    ifile.read_section("loads", nodal_loads_map);
    ifile.read_section("pressure", pressure_map);
    ifile.read_section("prescribed_displacement", fixed_nodes_map);
    ifile.close();

    // setup load control
    lc.set_nincs(num_increments);

    // setup elasticity problem type
    ElasticityType elastype;
    if (std::strncmp(etype.c_str(),"PLANE_STRAIN",12)==0)
      elastype = PLANE_STRAIN;
    else if(std::strncmp(etype.c_str(),"PLANE_STRESS",12)==0)
      elastype = PLANE_STRESS;
    else if(std::strncmp(etype.c_str(),"THREE_DIM",9)==0)
      elastype = THREE_DIM;
    else
      throw runtime_error("Unknown elasticity type.");

    // setup material type
    material = HyperelasticMaterial::create(mtype, elastype, matprop);
    cout << "Hyperelastic material: " << mtype << endl;

    // setup mesh filename
    assert(file_exists(mshfile.c_str()));
    std::string::size_type idx = mshfile.rfind(".msh");
    filename = mshfile;
    basename = filename.substr(0,idx);
  }

  if(extension == "xml")
  {
    cout << "Reading parameters from XML file" << endl;

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(parfile.c_str());

    if (!result) {
      std::cout << "XML parsed with errors\n";
      std::cout << "Error description: " << result.description() << "\n";
      std::cout << "Error offset: " << result.offset;
      std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
    }

    // setup elasticity problem type
    ElasticityType elastype;
    pugi::xml_attribute type = doc.child("elasticity").attribute("type");
    if (std::strncmp(type.as_string(), "PLANE_STRAIN", 12) == 0)
      elastype = PLANE_STRAIN;
    else if (std::strncmp(type.as_string(), "PLANE_STRESS", 12) == 0)
      elastype = PLANE_STRESS;
    else if (std::strncmp(type.as_string(), "THREE_DIM", 9) == 0)
      elastype = THREE_DIM;

    // number of increments

    pugi::xml_node params = doc.child("elasticity").child("parameters");
    pugi::xml_node reg = params.child("regions");

    if (params.attribute("num_materials"))
    {
      std::vector<int> map_matAHA;
      //std::cout << "Regions: " << std::endl;
      for (pugi::xml_node m = reg.child("marker"); m;
           m = m.next_sibling("marker"))
      {
        int id = m.attribute("id").as_int();
        int mt = m.attribute("material").as_int();
        map_matAHA.push_back(mt);
        //std::cout << "reg: " << id << " mat: " << mt << std::endl;
      }

        // load control
        int num_increments = params.child("ninc").text().as_int();
        lc.set_nincs(num_increments);

        int num_materials = params.attribute("num_materials").as_int();
        std::vector<std::vector<double>> vec_matprop;

        for (pugi::xml_node m = params.child("material"); m;
             m = m.next_sibling("material"))
        {
          int id = m.attribute("id").as_int();
          mtype = m.attribute("type").as_string();
          
          std::string strprop = m.child("coefficients").text().as_string();
          matprop.clear();
          parse_to_vector(strprop, matprop);
          std::cout << "Material " << id << std::endl;

          std::vector<double>::iterator it;
          for (it = matprop.begin(); it != matprop.end(); ++it)
            cout << *it << " ";
          cout << endl;

          vec_matprop.push_back(matprop);
        }

      // setup material
      std::cout << mtype << std::endl;
      material = HyperelasticMaterial::create(mtype, elastype, vec_matprop, num_materials, map_matAHA);
      cout << " Hyperelastic material: " << mtype << endl;
      cout << " Material properties: ";
      cout << " Multiple materials were defined. ";
    }
    else
    {
      if (params)
      {
        // load control
        int num_increments = params.child("ninc").text().as_int();
        lc.set_nincs(num_increments);

        // material type
        mtype = params.child("material").text().as_string();

        // material properties
        std::string strprop = params.child("coefficients").text().as_string();
        parse_to_vector(strprop, matprop);

        // setup material
        material = HyperelasticMaterial::create(mtype, elastype, matprop);
        cout << " Hyperelastic material: " << mtype << endl;
        cout << " Material properties: ";
        std::vector<double>::iterator it;
        for (it = matprop.begin(); it != matprop.end(); ++it)
          cout << *it << " ";
        cout << endl;
      }
    }

    // body force (ex: gravity)
    pugi::xml_node bf = doc.child("elasticity").child("body_forces");
    if (bf) {
      std::string strbf = bf.text().as_string();
      parse_to_vector(strbf, bforce);
    }

    //
    // reading boundary conditions
    //

    // neumann
    pugi::xml_node nodes = doc.child("elasticity").child("neumann");
    for (pugi::xml_node node = nodes.child("node"); node;
         node = node.next_sibling("node")) {
      //int idx = node.attribute("id").as_int();
      int marker = node.attribute("marker").as_int();
      arma::vec3 t;
      t(0) = node.attribute("t0").as_double();
      t(1) = node.attribute("t1").as_double();
      t(2) = node.attribute("t2").as_double();
      neumann_map.insert(std::pair<int, arma::vec3>(marker, t));
    }

    // dirichlet
    pugi::xml_node dnodes = doc.child("elasticity").child("dirichlet");
    for (pugi::xml_node node = dnodes.child("node"); node;
         node = node.next_sibling("node")) {
      int marker = node.attribute("marker").as_int();
      int dir = node.attribute("direction").as_int();
      dirichlet_map.insert(std::pair<int, int>(marker, dir));
    }

    // dirichlet - fixed nodes
    pugi::xml_node fnodes = doc.child("elasticity").child("prescribed_displacement");
    for (pugi::xml_node node = fnodes.child("node"); node;
         node = node.next_sibling("node"))
    {
      int idx = node.attribute("id").as_int();
      int dir = node.attribute("direction").as_int();
      double val = node.attribute("value").as_double();
      NodalData nd(dir, val);
      fixed_nodes_map.insert(std::pair<int, NodalData>(idx, nd));
    }

    // pressure (normal following load) boundary condition
    pugi::xml_node pnodes = doc.child("elasticity").child("pressure");
    for (pugi::xml_node node = pnodes.child("node"); node;
         node = node.next_sibling("node")) {
      int marker = node.attribute("marker").as_int();
      double val = node.attribute("value").as_double();
      pressure_map.insert(std::pair<int, double>(marker, val));
    }

    // pressure (normal following load) boundary condition
    pugi::xml_node snodes = doc.child("elasticity").child("spring");
    for (pugi::xml_node node = snodes.child("node"); node;
         node = node.next_sibling("node")) {
      int marker = node.attribute("marker").as_int();
      double val = node.attribute("value").as_double();
      spring_map.insert(std::pair<int, double>(marker, val));
    }

    // nodal loads
    pugi::xml_node lnodes = doc.child("elasticity").child("loads");
    for (pugi::xml_node node = lnodes.child("node"); node;
         node = node.next_sibling("node")) {
      int idx = node.attribute("id").as_int();
      int dir = node.attribute("direction").as_int();
      double val = node.attribute("value").as_double();
      NodalData nd(dir, val);
      nodal_loads_map.insert(std::pair<int, NodalData>(idx, nd));
    }

    // setup mesh filename
    assert(file_exists(mshfile.c_str()));
    std::string::size_type idx = mshfile.rfind(".xml");
    filename = mshfile;
    basename = filename.substr(0,idx);
  }

  // print some info
  cout << " Number of nodal loads: " << nodal_loads_map.size() << endl;
  cout << " Number of prescribed displ.: " << fixed_nodes_map.size() << endl;
  cout << " Number of traction (Neumann) loads: " << neumann_map.size() << endl;
  cout << " Number of dirichlet boundary conds:" << dirichlet_map.size() << endl;
  cout << " Number of normal pressure loads: " << pressure_map.size() << endl;
  cout << " Number of spring boundary conds: " << spring_map.size() << endl;
}

void NonlinearElasticity::elem_resid (const int iel, const MxFE * fe,
                                      const Quadrature * qd, arma::vec & Re)
{
  error("NonlinearElasticity is an abstract class.");
}

void NonlinearElasticity::elem_stiff (const int iel, const MxFE * fe,
                                      const Quadrature * qd, arma::mat & Ke)
{
  error("NonlinearElasticity is an abstract class.");
}

void NonlinearElasticity::evaluate_forces(petsc::Vector & R)
{
  int ndofs = fespace.get_ndofs();
  double dlamb = lc.load_step();

  //tload = fext;

  for(int i=0; i<ndofs; i++)
    if(ldgof[i])
      R.add(i, dlamb * fext[i]);
}

int NonlinearElasticity::get_num_nz_prescribed()
{
  int cont=0;
  double val;
  NodalData inf;
  std::map<int,NodalData>::const_iterator it;

  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    inf = it->second;
    val = inf.second;
    if (val != 0.0) cont++;
  }

  return cont;
}

int NonlinearElasticity::get_num_integration_points()
{
  int nint = 0;

  if (msh.get_nen() == 3 && msh.get_n_dim() == 2)
    nint = 3;
  else if (msh.get_nen() == 4 && msh.get_n_dim() == 2)
    nint = 4;
  else if (msh.get_nen() == 4 && msh.get_n_dim() == 3)
    nint = 1;
  else if (msh.get_nen() == 8 && msh.get_n_dim() == 3)
    nint = 8;

  return nint;
}

void NonlinearElasticity::get_element_x(const int e,
                                        std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_element_pt_nums(e, ptnums);

  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x[k];
  }
}

void NonlinearElasticity::get_element_x0(const int e,
                                         std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_element_pt_nums(e, ptnums);

  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x0[k];
  }
}

void NonlinearElasticity::get_boundary_element_x(const int e,
                                                 std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_boundary_element_pt_nums(e, ptnums);
  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x[k];
  }
}

void NonlinearElasticity::get_boundary_element_x0(const int e, std::vector<arma::vec3> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_boundary_element_pt_nums(e, ptnums);
  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x0[k];
  }
}

void NonlinearElasticity::get_displacements(arma::vec & u)
{
  for(uint i=0; i<msh.get_n_points(); i++)
  {
    double dx = x[i][0] - x0[i][0];
    double dy = x[i][1] - x0[i][1];
    double dz = x[i][2] - x0[i][2];
    u(0 + i*3) = dx;
    u(1 + i*3) = dy;
    u(2 + i*3) = dz;
  }
}

void NonlinearElasticity::get_displacements(arma::mat & umat)
{
  const int np = msh.get_n_points();
  const int nd = msh.get_n_dim();
  umat.reshape(np,nd);
  for(int d=0; d<nd; d++)
    for(int i=0; i<np; i++)
      umat(i,d) = x[i][d];
}

void NonlinearElasticity::init()
{
    init_mesh();
    init_matvecs();
}

void NonlinearElasticity::init_mesh()
{
  // initialize mesh
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

  cout << msh;
  fespace.set_mesh(&msh);

  // setup size of Voigt vec/mat
  nvoig = (msh.get_n_dim() == 2 ? 3 : 6);

  // initialize coordinates
  const std::vector<arma::vec3> & pts = msh.get_points();
  for(uint i=0; i<pts.size(); i++)
  {
    const arma::vec3 & p = pts[i];
    x.push_back(p);
    x0.push_back(p);
  }

  std::cout << "Output step: " << output_step << std::endl;
  // setup data writer
  if(output_step)
  {
    std::string output = filename.substr(0, filename.length() - 4) + "_output_nl";
    writer.open(output, lc.get_nincs() + 1 + INCs_TA, 1);
  }
}

void NonlinearElasticity::setup_data_writer(int size)
{
  std::cout << "Output step size: " << size << std::endl;
  std::string output = filename.substr(0, filename.length() - 4) + "_output_nl";
  writer.open(output, size + 1, 1);
}

void NonlinearElasticity::init_matvecs()
{
  int nod, dir, idx;
  int N = fespace.get_ndofs();
  int n_elem = msh.get_n_elements();

  NodalData inf;
  std::map<int,NodalData>::const_iterator it;

  num_dofs = N;

  fext.resize(N);  fext.zeros();
  fext0.resize(N); fext0.zeros();
  udisp.resize(N); udisp.zeros();
  U.resize(N); U.zeros();

  tload.resize(N); tload.zeros();

  // initialize ldgof array
  //  --> ldgof[i] = true, if dof i is active
  //  --> ldgof[i] = false, if dof i is prescribed
  for(int i=0; i<N; i++) ldgof.push_back(true);

  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    //idx = nod + (dir * msh.get_n_points());
    idx = (nod*msh.get_n_dim()) + dir;
    ldgof[idx] = false;
  }

  // initialize nodal loads vector
  for(it=nodal_loads_map.begin(); it!=nodal_loads_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    //idx = nod + (dir * msh.get_n_points());
    idx = (nod*msh.get_n_dim()) + dir;
    fext(idx) = inf.second;
  }

  // Allocate matrices and vectors
  cout << "N: " << N << endl;
  K.create(N,N,120); K = 0.0;
  u.create(N);       u = 0.0;
  r.create(N);       r = 0.0;
  react.create(N);   react = 0.0;

  // initialize stress array
  int nint = get_num_integration_points();
  stressdb.set_size(msh.get_n_elements(), nint, 9);
  straindb.set_size(msh.get_n_elements(), nint, 9);

  // initialize vector of F's
  for(int i=0; i<msh.get_n_elements(); i++)
    for(int j=0; j<nint; j++)
    {
      arma::mat33 * F = new arma::mat33();
      F->eye(3,3);
      vecF.push_back(F);
    }

  // resize vectors for 3 field formulation and ALG
  ep.resize(n_elem);
  eJ.resize(n_elem);
  eL.resize(n_elem);
  eL0.resize(n_elem);
  eps.resize(n_elem);
  eps0.resize(n_elem);
  vol0.resize(n_elem);

  if( material->is_incompressible() )
  {
    double penalt;
    penalt = (static_cast<IncompressibleMaterial*>(material))->get_kappa();

    eL0.fill(1.0);
    eL.zeros();
    eps.fill(penalt);
    eps0.fill(penalt);
  }
}

void NonlinearElasticity::init_resid_stiff()
{
  K = 0.0;
  react = 0.0;
  double xlamb = lc.load();
  for(int i=0; i<fespace.get_ndofs(); i++)
  {
    r.set(i, fext0(i) + xlamb * fext(i) );
  }
}

void NonlinearElasticity::line_search(double & eta0, double & eta,
                                      double & rtu0, double & rtu)
{
  double rtu1, alfa, q;

  rtu = calc_energy();
  if (rtu == 0.0) return;

  rtu1 = (rtu - rtu0 * (1 - eta))/(eta*eta);
  alfa = rtu0/rtu1;
  eta0 = eta;

  // Take different actions, depending on the value of alfa
  if (alfa < 0)
  {
    q = (alfa - sqrt(alfa * (alfa - 4)))/2.0;
    eta = alfa / q;
  }
  else if (alfa < 2)
  {
    eta = alfa/2.0;
  }
  else
  {
    eta = 1.0;
    rtu = 0.0;
  }
}

void NonlinearElasticity::output_vtk(const int cont, const int step)
{
  int np = msh.get_n_points();
  //int ne = msh.get_n_elements();
  std::stringstream ss;
  std::string name;
  std::string vtuname;

  //arma::vec t11(ne), t22(ne), t33(ne);
  //arma::vec t12(ne), t13(ne), t23(ne);
  //t11.zeros(); t22.zeros(); t33.zeros();
  //t12.zeros(); t13.zeros(); t23.zeros();

  ss << cont << "_" << step;
  name = this->basename + "_" + ss.str();
  vtuname = name + ".vtu";

  // Displacements
  //cout << "Aloca vetor" << endl;
  double *u_field = new double[np*3];
  //cout << "Preenche vetor" << endl;
  for(int i=0; i<np; i++)
  {
    double dx = x[i][0] - x0[i][0];
    double dy = x[i][1] - x0[i][1];
    double dz = x[i][2] - x0[i][2];
    u_field[i*3 + 0] = dx;
    u_field[i*3 + 1] = dy;
    u_field[i*3 + 2] = dz;
  }


  //if(output_step)
  //{
  //cout << "Salva hdf5" << endl;
    writer.write_displ_step(step, u_field);
  //}

  //cout << "Libera vetor" << endl;
  delete [] u_field;

  /*
  // Average stress output
  int nint = get_num_integration_points();
  for(int i=0; i<msh.get_n_elements(); i++)
  {
    for(int j=0; j<nint; j++)
    {
      t11(i) += stressdb(i,j,0);
      t12(i) += stressdb(i,j,1);
      t13(i) += stressdb(i,j,2);
      t22(i) += stressdb(i,j,4);
      t23(i) += stressdb(i,j,5);
      t33(i) += stressdb(i,j,8);
    }
    t11(i) /= nint;
    t12(i) /= nint;
    t13(i) /= nint;
    t22(i) /= nint;
    t23(i) /= nint;
    t33(i) /= nint;
  }

  vtkout.write_point_data(ux, "ux");
  vtkout.write_point_data(uy, "uy");
  vtkout.write_point_data(uz, "uz");
  vtkout.write_point_data(um, "umagnitude");
  vtkout.write_cell_data(ep,  "pressure");
  vtkout.write_cell_data(t11, "cauchy11");
  vtkout.write_cell_data(t22, "cauchy22");
  vtkout.write_cell_data(t33, "cauchy33");
  vtkout.write_cell_data(t12, "cauchy12");
  vtkout.write_cell_data(t13, "cauchy13");
  vtkout.write_cell_data(t23, "cauchy23");
  vtkout.write_def_mesh(x, vtuname);
  */
}

void NonlinearElasticity::output_vtk(const int step, const arma::vec & v, const arma::vec & displ)
{


  cout << "Escrevendo DADOS\n";
  writer.write_vm_step(step, v.memptr());
  writer.write_displ_step(step, displ.memptr());


}

void NonlinearElasticity::output_vtk(const int cont,
                                     const int step,
                                     const std::string & name,
                                     const arma::vec & v,
                                     const ArrayMat33 & matfib)
{
  int np = msh.get_n_points();

  //stringstream ss;
  //ss << "_" << step;
  //int idx = basename.rfind('/');
  //string outname = "output/" + basename.substr(idx + 1) + ss.str() + ".vtu";

  // First fiber directions
  arma::vec3 f, s, n;
  std::vector<arma::vec3> vec_f, vec_s, vec_n;

  for(uint i=0; i<matfib.size(); i++)
  {
    f = (matfib[i])->col(0);
    s = (matfib[i])->col(1);
    n = (matfib[i])->col(2);
    vec_f.push_back(f);
    vec_s.push_back(s);
    vec_n.push_back(n);
  }

  double *u_field = new double[np*3];
  for(int i=0; i<np; i++)
  {
    double dx = x[i][0] - x0[i][0];
    double dy = x[i][1] - x0[i][1];
    double dz = x[i][2] - x0[i][2];
    u_field[i*3 + 0] = dx;
    u_field[i*3 + 1] = dy;
    u_field[i*3 + 2] = dz;
  }

  if(output_step)
  {
    cout << "Escrevendo DADOS\n";
    writer.write_vm_step(step, v.memptr());
    writer.write_displ_step(step, u_field);
  }

  delete [] u_field;

  //vtkout.write_point_data(v, name);
  //vtkout.write_cell_vector(vec_f, "fiber");
  //vtkout.write_cell_vector(vec_s, "sheet");
  //vtkout.write_cell_vector(vec_n, "snorm");

  // then everything else
  //vtkout.write_def_mesh(x, outname);
}

void NonlinearElasticity::storeLVvolumes(string basename)
{
  ofstream file;
  string aux = basename + string("_volumes.dat");
  file.open(aux.c_str());

  file << total_volume_cavity() << "\n";
  file << calc_volume();

  file.close();
}

void NonlinearElasticity::storeStress(int step)
{
  //ofstream sfile, strainfile;
  //string aux = basename + string("_stress.dat");
  //sfile.open(aux.c_str());
  //aux = basename + string("_strain.dat");
  //strainfile.open(aux.c_str());
  arma::mat33 sig, E;
  int nint = get_num_integration_points();
  std::string vtuname;
  int ne = arma::size(stressdb, 0);

  arma::vec t11(ne), t22(ne), t33(ne);
  arma::vec t12(ne), t13(ne), t23(ne);
  arma::vec fiber_stress(ne), fiber_strain(ne), long_strain(ne), circ_strain(ne), rad_strain(ne);
  t11.zeros(); t22.zeros(); t33.zeros();
  t12.zeros(); t13.zeros(); t23.zeros();
  fiber_stress.zeros(); fiber_strain.zeros();
  long_strain.zeros(); circ_strain.zeros(); rad_strain.zeros();

  //cout << "Pontos de integração: " << nint << endl;

  vtuname = basename + "_stress.vtu";

  for (int id_el = 0; id_el < ne; id_el++)
  {
    for (int kk = 0; kk < 9; kk++) // stress component
    {
      double mediastress = 0;
      double mediaE = 0;
      for (int ii = 0; ii < nint; ii++) // integration point
      {
        mediastress += stressdb(id_el, ii, kk);
        mediaE += straindb(id_el, ii, kk);
      }
      mediastress = mediastress / (double) nint;
      //sfile.setf(ios::scientific);
      //sfile << mediastress << "\t";
      sig[kk] = mediastress;

      mediaE = mediaE / (double) nint;
      //strainfile.setf(ios::scientific);
      //strainfile << mediaE << "\t";
      E[kk] = mediaE;
    }

    t11(id_el) = sig(0,0);
    t12(id_el) = sig(0,1);
    t13(id_el) = sig(0,2);
    t22(id_el) = sig(1,1);
    t23(id_el) = sig(1,2);
    t33(id_el) = sig(2,2);

    arma::mat33 F, *auxF;
    F.zeros();
    for (int ii = 0; ii < nint; ii++) // integration point
      {
        auxF = vecF[id_el*nint + ii];
        F = F + *auxF;
      }
    F = F/nint;

    arma::vec3 fib0 = msh.get_element(id_el).get_fiber();
    arma::vec3 fib = msh.get_element(id_el).get_fiber();
    arma::vec3 f_l = msh.get_element(id_el).get_long();
    arma::vec3 f_c = msh.get_element(id_el).get_circ();
    arma::vec3 f_r = msh.get_element(id_el).get_rad();


    //std::cout << f_l << " " << f_c << " " << f_r << std::endl;

    fib = F*fib; //F*fib; //Alterei aqui para testar stress longitudinal
    fib = fib/arma::norm(fib,2);

    fiber_stress(id_el) = fib(0) * (fib(0)*sig(0,0) + fib(1)*sig(0,1) + fib(2)*sig(0,2)) +
                          fib(1) * (fib(0)*sig(0,1) + fib(1)*sig(1,1) + fib(2)*sig(1,2)) +
                          fib(2) * (fib(0)*sig(0,2) + fib(1)*sig(1,2) + fib(2)*sig(2,2));
    //sfile << fiber_stress(id_el) << "\t";
    //sfile << endl;

    fiber_strain(id_el) = fib0(0) * (fib0(0)*E(0,0) + fib0(1)*E(0,1) + fib0(2)*E(0,2)) +
                          fib0(1) * (fib0(0)*E(0,1) + fib0(1)*E(1,1) + fib0(2)*E(1,2)) +
                          fib0(2) * (fib0(0)*E(0,2) + fib0(1)*E(1,2) + fib0(2)*E(2,2));

    long_strain(id_el) =  f_l(0) * (f_l(0)*E(0,0) + f_l(1)*E(0,1) + f_l(2)*E(0,2)) +
                          f_l(1) * (f_l(0)*E(0,1) + f_l(1)*E(1,1) + f_l(2)*E(1,2)) +
                          f_l(2) * (f_l(0)*E(0,2) + f_l(1)*E(1,2) + f_l(2)*E(2,2));

    circ_strain(id_el) =  f_c(0) * (f_c(0)*E(0,0) + f_c(1)*E(0,1) + f_c(2)*E(0,2)) +
                          f_c(1) * (f_c(0)*E(0,1) + f_c(1)*E(1,1) + f_c(2)*E(1,2)) +
                          f_c(2) * (f_c(0)*E(0,2) + f_c(1)*E(1,2) + f_c(2)*E(2,2));

    rad_strain(id_el) =   f_r(0) * (f_r(0)*E(0,0) + f_r(1)*E(0,1) + f_r(2)*E(0,2)) +
                          f_r(1) * (f_r(0)*E(0,1) + f_r(1)*E(1,1) + f_r(2)*E(1,2)) +
                          f_r(2) * (f_r(0)*E(0,2) + f_r(1)*E(1,2) + f_r(2)*E(2,2));

    //std::cout << fiber_strain(id_el) << " " << long_strain(id_el) << " " << circ_strain(id_el) << " " << rad_strain(id_el) << std::endl;

    //strainfile << fiber_strain(id_el) << "\t";
    //strainfile << endl;
  }
  //sfile.close();
  //strainfile.close();

  writer.write_cell_field_step(step, fiber_stress.memptr(), string("stress"));
  writer.write_cell_field_step(step, fiber_strain.memptr(), string("strain"));
  writer.write_cell_field_step(step, long_strain.memptr(), string("long_strain"));
  writer.write_cell_field_step(step, circ_strain.memptr(), string("circ_strain"));
  writer.write_cell_field_step(step, rad_strain.memptr(), string("rad_strain"));

}


void NonlinearElasticity::prescribe_displacements()
{
  int nod, dir, ndim = msh.get_n_dim();
  double val;
  NodalData inf;
  std::map<int,double> boundary_values;
  std::map<int,NodalData>::iterator it;

  // loop over nodes with prescribed displacement
  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    val = inf.second;
    assert(dir+1 <= ndim);
    x[nod][dir] = x0[nod][dir] + lc.load() * val;
  }
}

void NonlinearElasticity::reset()
{
  // Gambiarra
  K.assemble();
  //--------------
  //for(int i=0; i<x.size(); i++) x0[i] = x[i];
  //r=0.0;
  //fext.zeros();
  //fext0.zeros();

  K = 0.0;
  u = 0.0;
  udisp.zeros();
  lc.reset();
}

void NonlinearElasticity::set_pressure_Ta(int mlv, double plv, int mrv, double prv, const arma::vec &ta, const arma::vec &dta)
{
  pressure_map[mlv] = plv;
  pressure_map[mrv] = prv;
  material->set_Ta(ta);
  material->set_dTa(dta);
}

void NonlinearElasticity::run()
{
  init();
  pre_solve();
  cout << "Initial cavity volume: " << total_volume_cavity() << "\n";
  solve();
  storeStress(lc.get_nincs());
  storeLVvolumes(this->basename);
  //cout << "Final cavity volume: " << total_volume_cavity() << "\n";
  timer.summary();
}

void NonlinearElasticity::update_geometry(double eta)
{
  const int np = msh.get_n_points();
  const int nd = msh.get_n_dim();
  arma::mat umat;

  // get data from PETSc vector and copy to udisp
  u.get_data( udisp.memptr() );
  umat = arma::reshape(udisp, np, nd);
  U = U + udisp;

  // update geometry
  for(int d=0; d<nd; d++)
    for(int i=0; i<np; i++)
    {
      //int idx = i + (d*np);
      int idx = (i*nd) + d;
      // update if dof is free (not fixed)
      if( ldgof[idx] )
        x[i][d] = x[i][d] + eta * umat(i,d);
    }
}

void NonlinearElasticity::update_vectors(const ArrayMat33 & matfib0,
                                         ArrayMat33 & matfib)
{
  int ndim, ndofs, nnode, nelem;
  MixedFiniteElement * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());

  //Quadrature * qd = Quadrature::create(2, fe->get_type());

  ndim  = msh.get_n_dim();
  nelem = msh.get_n_elements();
  ndofs = fe->get_ndofs_u();
  nnode = fe->get_nnode();

  for(int iel=0; iel<nelem; iel++)
  {
    arma::mat dshape;
    arma::mat gradn(ndofs,ndim);
    arma::mat jacinv(ndim,ndim);
    arma::vec3 centroid;
    std::vector<arma::vec3> xe(ndofs);
    std::vector<arma::vec3> x0(ndofs);
    get_element_x(iel, xe);
    get_element_x0(iel, x0);

    Mapping em = fe->get_mapping(iel, x0);

    // compute centroid
    //   works for quad
    //   works for hex
    centroid.zeros();

    // TODO: check for each reference element
    // TODO: need to map c to the reference element ?

    // compute shape and shape derivatives at centroid
    fe->calc_deriv_shape_u(centroid, dshape);
    em.calc_jacobian(dshape, nnode);
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    // compute F
    arma::mat33 F(arma::fill::zeros);

    F(2,2) = (ndim==2) ? 1.0 : 0.0;

    // F = dx/dX
    for(int id=0; id<ndim; id++)
      for(int jd=0; jd<ndim; jd++)
	      for(int k=0; k<nnode; k++)
	        F(id,jd) += xe[k][id] * gradn(k,jd);

    // update vector directions
    arma::vec3 f,s,n;

    f = F * (matfib0[iel])->col(0);
    s = F * (matfib0[iel])->col(1);
    n = F * (matfib0[iel])->col(2);

    // normalize vectors
    f = f / arma::norm(f,2);
    s = s / arma::norm(s,2);
    n = n / arma::norm(n,2);

    //arma::mat33 R,U;
    //polar_decomposition(F,R,U);
    //f = R * (matfib0[iel])->col(0);
    //s = R * (matfib0[iel])->col(1);
    //n = R * (matfib0[iel])->col(2);

    (matfib[iel])->col(0) = f;
    (matfib[iel])->col(1) = s;
    (matfib[iel])->col(2) = n;
  }

  delete qd;
  delete fe;
}

// ---------- NONLINEAR PROBLEM ------------------------------------------------

bool NonlinearElasticity::converged(petsc::Vector & du)
{
  double rnorm=0, enorm=0, dnorm=0;
  double etol   = parameters["tol_energy"];
  double rtol   = parameters["tol_residual"];
  double dtol   = parameters["tol_displacement"];
  //double enorm0 = parameters["energy_norm0"];
  //double rnorm0 = parameters["residual_norm0"];

  // TROCANDO fext por tload
  // tload - external load including pressure
  // fext - external load

  if(first_step)
  {
    double f0 = react.l2norm() + arma::norm(lc.load()*tload,2);
    u.copy_values(du.size(), du);
    parameters["energy_norm0"]   = fabs(calc_energy());
    parameters["residual_norm0"] = f0;

    rnorm = r.l2norm();
    rnorm = rnorm/f0;
    enorm = fabs(calc_energy());
    dnorm = du.l2norm();

    first_step = false;
  }
  else
  {
    rnorm = r.l2norm() / ( react.l2norm() + arma::norm(tload,2) ); // R/F (F=force+react)
    enorm = fabs(calc_energy());
    dnorm = du.l2norm();

    cout << scientific << setprecision(3);
    cout << " rnorm " << rnorm;
    cout << " enorm " << enorm;
    cout << " dunorm " << dnorm;
    cout << endl;
 }

  //cout << endl;
  //cout << "fext norm = " << arma::norm(fext,2) << endl;
  //cout << "tload norm = " << arma::norm(tload,2) << endl;

  // convergence criterion
  bool check;
  //rnorm = rnorm/(rnorm0 + 1e-10);
  //check  = fabs(enorm) < (etol * fabs(enorm0));
  //check |= fabs(rnorm) < (rtol * fabs(rnorm0));
  //check |= (dnorm < dtol);
  check  = enorm < etol;
  check |= rnorm < rtol;
  check |= dnorm < dtol;
  return check;
}

void NonlinearElasticity::evaluate(petsc::Vector & resid)
{
  timer.enter("Residual");

  int n_dofs = msh.get_n_dim() * msh.get_nen();
  int n_elem = msh.get_n_elements();

  // init resid
  double xlamb = lc.load();
  for(int i=0; i<fespace.get_ndofs(); i++)
    r.set(i, fext0(i) + xlamb * fext(i));

  // set tload (external loads incl. pressure) to current fext (external loads)
  tload = fext;

  //
  // assemble residual vector r
  //
  arma::vec Re(n_dofs);
  std::vector<int> dnums;
  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());
  //Quadrature * qd = Quadrature::create(2, fe->get_type());
  for(int i=0; i<n_elem; i++)
  {
    elem_resid (i, fe, qd, Re);
    fespace.get_element_dofs_u (i, dnums);

    for(int k=0; k<n_dofs; k++)
    {
      if (ldgof[dnums[k]])
      {
        r.add(dnums[k], -Re(k)); // add -R
      }
      else
        react.add(dnums[k], Re(k));

    }
  }
  delete qd;
  delete fe;

  //
  // pressure forces contribution
  //
  MixedFiniteElement * bfe = fespace.create_boundary_FE();
  if (bfe != NULL && pressure_map.size() > 0)
  {
    //cout << endl << "Assemble pressure forces" << endl;
    int nu  = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    double xlamb = lc.load();
    arma::vec belvec(nu);
    vector<int> bdof;

    for(int i=0; i<nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i,bdof);
      elem_pforce(i,bfe,bdof,belvec);

      // assembles the nodal forces due to normal pressure
      for(int k=0; k<nu; k++)
      {
        if(ldgof[bdof[k]])
        {
          r.add(bdof[k], xlamb * belvec(k));
          tload(bdof[k]) += belvec(k);
        }
        else
          react.add(bdof[k], -xlamb * belvec(k));
      }
    }
  }

  r.assemble();

  // copy from r to resid
  resid.copy_values(r.size(), r);
  resid.assemble();

  timer.leave();
}

void NonlinearElasticity::jacobian(petsc::Matrix & Kstiff)
{

  timer.enter("Stiffness");

  int n_dofs = msh.get_n_dim() * msh.get_nen();
  int n_elem = msh.get_n_elements();
  arma::mat Ke(n_dofs,n_dofs);
  std::vector<int> dnums;

  Kstiff = 0.0;

  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());
  //Quadrature * qd = Quadrature::create(2, fe->get_type());

  for(int i=0; i<n_elem; i++)
  {

    elem_stiff (i, fe, qd, Ke);
    fespace.get_element_dofs_u (i, dnums);
#ifndef USE_BFGS
    // Fast assembling
    int * pidx;
    pidx = &dnums[0];
    Ke = Ke.t();
    Kstiff.add(n_dofs, n_dofs, pidx, pidx, Ke.memptr());
#endif

#ifdef USE_BFGS
    for(int j=0; j<n_dofs; j++)
    {
      for(int k=0; k<n_dofs; k++)
      {
        int I = dnums[j];
        int J = dnums[k];
        if(J >= I) Kstiff.add(I, J, Ke(j,k));
      }
    }
#endif
  }

  //
  // *** PRESSURE COMPONENT OF THE STIFFNESS MATRIX ***
  //
  MxFE * bfe = fespace.create_boundary_FE();
  if (bfe != NULL && (pressure_map.size() > 0 || spring_map.size()>0))
  {
    //cout << " Assembling pressure component of stiffness matrix" << endl;
    int nu = bfe->get_ndofs_u();
    int nb = msh.get_n_boundary_elements();
    arma::mat belmat(nu, nu);
    vector<int> bdof;

    for (int i = 0; i < nb; i++)
    {
      fespace.get_boundary_element_dofs_u(i, bdof);
      elem_kpress(i, bfe, bdof, belmat);
      // assembles the pressure component of the stiffness matrix
      for (int j = 0; j < nu; j++)
        for (int k = 0; k < nu; k++)
          Kstiff.add(bdof[j], bdof[k], belmat(j, k));
    }
  }

  delete qd;
  delete fe;

  Kstiff.assemble();
  apply_boundary(Kstiff);

  /*
  int tam = 3*x.size();
  double *values = new double[tam], *p;
  int *idx = new int[tam], id=0;

  for(std::vector<arma::vec3>::iterator j=x.begin(); j<x.end(); j++)
  {

    idx[3*id + 0] = 3*id + 0;
    idx[3*id + 1] = 3*id + 1;
    idx[3*id + 2] = 3*id + 2;

    p = j->memptr();
    values[3*id + 0] = p[0];
    values[3*id + 1] = p[1];
    values[3*id + 2] = p[2];

    id++;
  }

  petsc::Vector *coord = new petsc::Vector();
  coord->create(tam);
  coord->set(tam, idx, values);
  Kstiff.setNullSpace(coord);
 */

  timer.leave();
}

void NonlinearElasticity::update(petsc::Vector & uu, double s)
{
  const int n_dim   = msh.get_n_dim();
  const int n_nodes = msh.get_n_points();
  arma::mat umat;

    arma::vec udispB;

    // get data from PETSc vector and copy to udisp
    uu.get_data( udisp.memptr() );
    udispB.resize(fespace.get_ndofs());
    udispB.zeros();
    for(int i=0; i<n_nodes; i++)
    {
        udispB(i) = udisp(n_dim*i);
        udispB(i + n_nodes) = udisp(n_dim*i + 1);
        udispB(i + 2*n_nodes) = udisp(n_dim*i + 2);
    }
    umat = arma::reshape(udispB, n_nodes, n_dim);

    // update total displacement vector
    U = U + udisp;

  // copy uu to internal u
  u.copy_values(uu.size(), uu);

  // update geometry
  for(int d=0; d<n_dim; d++)
    for(int i=0; i<n_nodes; i++)
    {
      //int idx = i + (d * n_nodes);
      int idx = (i*n_dim) + d;

      // if degree of freedom is free (not fixed/prescribed)
      // then update x = x + u
      if( ldgof[idx] )
        x[i][d] = x[i][d] + s * umat(i,d);
    }
}




/*
void NonlinearElasticity::update(petsc::Vector & uu, double s)
{
  const int n_dim   = msh.get_n_dim();
  const int n_nodes = msh.get_n_points();
  arma::mat umat;

  // get data from PETSc vector and copy to udisp
  uu.get_data( udisp.memptr() );
  umat = arma::reshape(udisp, n_nodes, n_dim);

  // update total displacement vector
  U = U + udisp;

  // copy uu to internal u
  u.copy_values(uu.size(), uu);

  // update geometry
  for(int d=0; d<n_dim; d++)
    for(int i=0; i<n_nodes; i++)
    {
      //int idx = i + (d * n_nodes);
      int idx = (i*n_dim) + d;

      // if degree of freedom is free (not fixed/prescribed)
      // then update x = x + u
      if( ldgof[idx] )
        x[i][d] = x[i][d] + s * umat(i,d);
    }
}
*/
