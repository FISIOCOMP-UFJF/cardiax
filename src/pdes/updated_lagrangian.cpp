#include "updated_lagrangian.hpp"

double UpdatedLagrangian::al_augment(double tol)
{
  int ne = msh.get_n_elements();
  double he;
  double barJ;
  double alnorm;

  IncompressibleMaterial *imat = (IncompressibleMaterial *)material;

  alnorm = 0.0;
  for (int ie = 0; ie < ne; ie++)
  {
    barJ = eJ(ie);
    he = imat->h(barJ);

    // Calculate |h(Jbar)| to check tol
    // alnorm += eL(ie) * eL(ie);
    alnorm += he * he;

    // Update Lambda if the elements violates |detF-1.0| >= tol
    // See Brinkhues 2013
    if (fabs(barJ - 1.0) > tol)
    {
      eL(ie) += eps(ie) * (barJ - 1.0);
    }

    // Old
    // eL(ie) += eps(ie) * he; // or eps(ie) * (barJ - 1.0);
  }

  alnorm = sqrt(alnorm);
  return alnorm;
}

void UpdatedLagrangian::al_update(int increment)
{
  if (material->is_incompressible())
  {
    int ne = msh.get_n_elements();

    if (increment == 1)
    {
      // Store Lambda^1
      for (int ie = 0; ie < ne; ie++)
        eL0(ie) = eL(ie);
    }
    else
    {
      // Update the value of the penalty parameter
      for (int ie = 0; ie < ne; ie++)
      {
        // kappa adaptativo
        // eps(ie) = fabs( eps0(ie) * ( eL(ie) / eL0(ie) ) );
      }
    }
  }
}

void UpdatedLagrangian::assemble_const()
{
  int n = msh.get_n_dim() * msh.get_nen();
  int ne = msh.get_n_elements();
  int *pidx;
  std::vector<int> dnums;
  arma::mat Ke(n, n);

  MxFE *fe = fespace.createFE();
  Quadrature *qd = Quadrature::create(0, fe->get_type());
  for (int i = 0; i < ne; i++)
  {
    calc_elmat_const(i, fe, qd, Ke);
    fespace.get_element_dofs_u(i, dnums);
    // Fast assembling
    pidx = &dnums[0];
    Ke = Ke.t();
    K.add(n, n, pidx, pidx, Ke.memptr());

  }

  K.assemble();
  apply_boundary(K);
  
  delete qd;
  delete fe;
}

void UpdatedLagrangian::assemble_stiff()
{
  int n = msh.get_n_dim() * msh.get_nen();
  int ne = msh.get_n_elements();
  int *pidx;
  std::vector<int> dnums;
  arma::mat Ke(n, n);
  arma::vec Re(n);

  MixedFiniteElement *fe = fespace.createFE();
  Quadrature *qd = Quadrature::create(0, fe->get_type());

  for (int i = 0; i < ne; i++)
  {
    calc_elmatvec(i, fe, qd, Ke, Re);
    fespace.get_element_dofs_u(i, dnums);

    // fast assembling
    pidx = &dnums[0];
    Ke = Ke.t();
    K.add(n, n, pidx, pidx, Ke.memptr());

    // assemble residual vector and reaction vector
    for (int k = 0; k < n; k++)
    {
      if (ldgof[dnums[k]])
        r.add(dnums[k], -Re(k)); // -R
      else
        react.add(dnums[k], Re(k));
    }
  }

  // Assemble matrix and vector and impose presc. displacements
  K.assemble();
  r.assemble();
  apply_boundary(K);

  delete qd;
  delete fe;
}

void UpdatedLagrangian::assemble_active(const arma::vec &is,
                                        std::vector<arma::mat33 *> &vstrs,
                                        std::vector<arma::mat33 *> &vfibcoords)
{
  MixedFiniteElement *fe = fespace.createFE();
  Quadrature *qd = Quadrature::create(0, fe->get_type());

  int n = msh.get_n_dim() * msh.get_nen();
  int ne = msh.get_n_elements();
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nint = qd->get_num_ipoints();
  int nnode = ndof / ndim;

  arma::vec tmp_fext(fespace.get_ndofs());
  arma::vec elvec(n);
  std::vector<int> dnums;

  tmp_fext.zeros();

  for (int iel = 0; iel < ne; iel++)
  {
    // Begin of element calculation
    double detJxW, J;
    arma::vec shape;
    arma::mat sigma_voigt(nvoig, 1), B(nvoig, ndof);
    arma::mat dshape, gradn(ndof, ndim), jacinv(ndim, ndim);
    arma::mat33 b, *F;
    arma::vec3 qpt;
    std::vector<arma::vec3> xe(ndof);
    std::vector<arma::vec3> x0(ndof);
    get_element_x(iel, xe);
    get_element_x0(iel, x0);

    // Tava assim!!!//////////////////////////////////
    Mapping em = fe->get_mapping(iel, x0);

    // TESTE...estava com o x0 antes
    // Mapping em = fe->get_mapping(iel, xe);
    //////////////////////////////////////////////////

    fespace.get_element_dofs_u(iel, dnums);

    elvec.zeros();

    for (int q = 0; q < nint; q++)
    {
      qpt = qd->get_point(q);
      fe->calc_shape_u(qpt, shape);
      fe->calc_deriv_shape_u(qpt, dshape);
      em.calc_jacobian(dshape, nnode);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();
      jacinv = em.get_inv_jacobian();
      gradn = (dshape * jacinv);

      F = vecF[iel * nint + q];
      lcg_tensor(gradn, x0, &J, F, b);
      calc_B_matrix(gradn, B);

      // interpolation of active stress from nodal to integration points
      arma::mat33 T;                       // active stress tensor
      arma::mat33 Tl(arma::fill::zeros);   // local active stress tensor
      arma::mat33 &M = *(vfibcoords[iel]); // change of basis matrix
      arma::mat33 &Ta = *(vstrs[iel]);     // anisotropic Cauchy active stress
      arma::mat33 Fi = arma::inv(*F);

      const arma::vec3 f = (*F) * M.col(0);

      Ta.zeros();

      for (int j = 0; j < nnode; j++)
      {
        int k = dnums[j]; // global node number

        // Tl(0,0) = is(k);     // local active tension in fiber direction
        // T = M * Tl * M.t();  // global active tension in fiber direction
        T = is(k) * (f * f.t());

        // interpolate
        for (int ii = 0; ii < 3; ii++)
          for (int jj = 0; jj < 3; jj++)
            Ta(ii, jj) += shape(j) * T(ii, jj);
      }

      // UL uses Sigma!!!
      // voigtvec(ndim, Ta, sigma_voigt);

      // Ta is Cauchy stress, convert to SPK
      arma::mat33 S = J * (Fi * Ta * Fi.t());
      // arma::mat33 S = J * ((*F) * Ta * F->t());
      voigtvec(ndim, S, sigma_voigt);

      // indicial notation
      // for(int i=0; i<nnode; i++)
      //  for(int id=0; id<ndim; id++)
      //    for(int jd=0; jd<ndim; jd++)
      //      elvec(i+id*nnode) += Ta(id,jd) * gradn(i,jd) * detJxW;

      // matrix form
      elvec += (B.t() * sigma_voigt) * detJxW;
    }

    // assemble into vector
    for (int k = 0; k < n; k++)
      if (ldgof[dnums[k]])
        tmp_fext(dnums[k]) -= elvec(k);
  }

  // copy
  fext += tmp_fext;

  delete qd;
  delete fe;
}

void UpdatedLagrangian::calc_elmat_const(const int iel, const MxFE *fe,
                                         const Quadrature *qd, arma::mat &Ke)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof / ndim;
  double detJxW, press;

  arma::mat dshape;
  arma::mat gradn(ndof, ndim);
  arma::mat jacinv(ndim, ndim);
  arma::mat B(nvoig, ndof), D(nvoig, nvoig);
  arma::vec3 qpt;
  arma::mat F = arma::eye<arma::mat>(3, 3);

  std::vector<arma::vec3> xe(nubf);

  get_element_x(iel, xe);

  Mapping em = fe->get_mapping(iel, xe);

  Tensor4 elastensor, elastensorM, elastensorA, elastensorM_global;
  elastensor.zero();
  elastensorM.zero();
  elastensorA.zero();

  Ke.zeros();
  press = 0;

  // Mean dilatation method
  if (material->is_incompressible())
    mean_dilatation(iel, fe, qd, press, Ke);


  // loop over integration points
  for (int q = 0; q < qd->get_num_ipoints(); q++)
  {
    qpt = qd->get_point(q);
    fe->calc_deriv_shape_u(qpt, dshape);
    em.calc_jacobian(dshape, nubf);
    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    calc_B_matrix(gradn, B);

    MaterialData *md = new MaterialData(msh.get_element(iel), F);
    if (material->is_incompressible())
    {

      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);
      // elasticity tensor - deviatoric part
      //

      if (im->get_name() == "Guccione")
      {
        im->deviatoric_elastensor(md, elastensorM);
        im->map_elas_to_global(md, elastensorM, elastensorM_global);

        if (md->get_marker() == 0)
          im->active_stress_elastensor(iel, lc.get_nincs(), md->fiber(), elastensorM_global); //Deu pau aqui

        im->push_forward(F, elastensorM_global, elastensor);
      }
      else
      {
        im->calc_fd_elastensor(iel, md, elastensorM);
        // im->active_stress_elastensor(lc.get_nincs(), md, elastensorM);
        im->push_forward(F, elastensorM, elastensor);
      }

      // converts to Voigt notation
      elastensor.get_matrix(D);
    }
    else
    {
      material->sp_elastensor(md, D);
    }

    // calculates the initial contribution to the stiffness matrix
    Ke += (B.t() * D * B) * detJxW;

    delete md;
  }
}

void UpdatedLagrangian::calc_elmatvec(const int iel, const MxFE *fe,
                                      const Quadrature *qd, arma::mat &Ke,
                                      arma::vec &Re)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof / ndim;
  int nint = qd->get_num_ipoints();
  double detJxW, detF, press;

  arma::vec shape;
  arma::mat dshape;
  arma::mat gradn(ndof, ndim);
  arma::mat jacinv(ndim, ndim);
  arma::mat sigma_voigt(nvoig, 1);
  arma::mat sigma2(ndim * ndim, ndim * ndim);
  arma::mat B(nvoig, ndof), D(nvoig, nvoig), H(ndim * ndim, ndof);
  arma::mat33 btens, *F, S, Sa, sigma;
  arma::vec3 pt, qpt;

  Tensor4 elastensorM, elastensorA, elastensor, elastensorM_global;

  elastensorA.zero();
  elastensorM.zero();

  std::vector<arma::vec3> elem_xe(nubf);
  std::vector<arma::vec3> elem_x0(nubf);
  get_element_x(iel, elem_xe);
  get_element_x0(iel, elem_x0);

  Mapping em = fe->get_mapping(iel, elem_xe);

  Ke.zeros();
  Re.zeros();
  press = 0;

  // mean dilatation method
  if (material->is_incompressible())
    mean_dilatation(iel, fe, qd, press, Ke);

  // loop over integration points
  for (int q = 0; q < qd->get_num_ipoints(); q++)
  {
    qpt = qd->get_point(q);
    fe->calc_shape_u(qpt, shape);
    fe->calc_deriv_shape_u(qpt, dshape);
    em.calc_jacobian(dshape, nubf);
    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    F = vecF[iel * nint + q];

    lcg_tensor(gradn, elem_x0, &detF, F, btens);
    calc_B_matrix(gradn, B);
    calc_H_matrix(gradn, H);

    // compute stress and elasticity tensor
    MaterialData *md = new MaterialData(msh.get_element(iel), *F);

    if (material->is_incompressible())
    {
      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);

      // im->deviatoric_stress(md, sigma);
      // im->add_pressure(press, sigma);
      // im->deviatoric_elastensor(md, elastensor);
      // im->sp_volumetric_elastensor(press, elastensor);
      // elastensor.get_matrix(D);

      // stress
      // im->calc_fd_stress(md, S);
      if (im->get_name() == "Guccione")
      {
        im->piola2_stress(md, S);
        if (md->get_marker() == 0)
          im->set_active_stress(iel, S, lc.load());
        im->map_to_global(md, S);
        im->push_forward(*F, S, sigma);
      }
      else
      {
        im->calc_fd_stress(iel, md, S);
        // im->set_active_stress(md, S, lc.load());
        im->push_forward(*F, S, sigma);
        // im->set_active_stress(md, sigma, lc.load());
      }

      im->add_pressure(press, sigma);

      // elasticity tensor - deviatoric part
      // im->calc_fd_elastensor(md, elastensorM);

      // TODO: improve this, needs to be general
      if (im->get_name() == "Guccione")
      {
        im->deviatoric_elastensor(md, elastensorM);
        im->map_elas_to_global(md, elastensorM, elastensorM_global);
        if (md->get_marker() == 0)
          im->active_stress_elastensor(iel, lc.get_nincs(), md->fiber(), elastensorM_global);
        im->push_forward(*F, elastensorM_global, elastensor);
      }
      else
      {
        im->calc_fd_elastensor(iel, md, elastensorM);
        // im->active_stress_elastensor(lc.get_nincs(), md, elastensorM);
        im->push_forward(*F, elastensorM, elastensor);
      }

      // elasticity tensor - volumetric part
      im->sp_volumetric_elastensor(press, elastensor);
      // converts to Voigt notation
      elastensor.get_matrix(D);
    }
    else
    {
      material->cauchy_stress(md, sigma);
      material->sp_elastensor(md, D);
    }

    // create and convert sigma stress tensor to
    // appropriate formats for matricial manipulation
    sigma2 = arma::kron(arma::eye(ndim, ndim), sigma.submat(0, 0, ndim - 1, ndim - 1));
    voigtvec(ndim, sigma, sigma_voigt);

    // internal nodal forces
    Re += (B.t() * sigma_voigt) * detJxW;

    // material stiffness matrix - constitutive contribution (elmat_const)
    Ke += (B.t() * D * B) * detJxW;

    // geometrical stiffness matrix - stress contribution (elmat_sigma)
    Ke += (H.t() * sigma2 * H) * detJxW;

    // store stresses
    stressdb.tube(iel, q) = arma::vectorise(sigma, 1);
    straindb.tube(iel, q) = arma::vectorise(md->lagrangian_strain(), 1);
    // S = arma::det(*F) * arma::inv(*F) * sigma * arma::inv(*F).t();
    // stressdb.tube(iel, q) = arma::vectorise(S,1);

    delete md;
  }
}

inline void UpdatedLagrangian::calc_B_matrix(const arma::mat &gradn,
                                             arma::mat &B) const
{
  B.zeros();

  int Nr = gradn.n_rows;    // ndof
  int Nd = msh.get_n_dim(); // ndim
  int Nb = Nr / Nd;         // nubf

  if (Nd == 2)
  {
    // e,x
    for (int i = 0; i < Nb; i++)
      B(0, i) = gradn(i, 0);
    // e,y
    for (int i = Nb; i < 2 * Nb; i++)
      B(1, i) = gradn(i, 1);
    // ex,y
    for (int i = 0; i < Nb; i++)
      B(2, i) = gradn(i, 1);
    for (int i = Nb; i < 2 * Nb; i++)
      B(2, i) = gradn(i, 0);
  }
  else if (Nd == 3)
  {
    // e,x
    for (int i = 0; i < Nb; i++)
      B(0, i) = gradn(i, 0);
    // e,y
    for (int i = Nb; i < 2 * Nb; i++)
      B(1, i) = gradn(i, 1);
    // e,z
    for (int i = 2 * Nb; i < 3 * Nb; i++)
      B(2, i) = gradn(i, 2);

    // ex,y
    for (int i = 0; i < Nb; i++)
      B(3, i) = gradn(i, 1);
    for (int i = Nb; i < 2 * Nb; i++)
      B(3, i) = gradn(i, 0);
    // ey,z
    for (int i = Nb; i < 2 * Nb; i++)
      B(4, i) = gradn(i, 2);
    for (int i = 2 * Nb; i < 3 * Nb; i++)
      B(4, i) = gradn(i, 1);
    // ex,z
    for (int i = 0; i < Nb; i++)
      B(5, i) = gradn(i, 2);
    for (int i = 2 * Nb; i < 3 * Nb; i++)
      B(5, i) = gradn(i, 0);
  }
}

void UpdatedLagrangian::calc_H_matrix(const arma::mat &gradn,
                                      arma::mat &H) const
{
  H.zeros();

  int ndofs = gradn.n_rows;    // number of dofs
  int ndime = msh.get_n_dim(); // number of dimensions
  int nnode = ndofs / ndime;   // number of nodes

  if (ndime == 1) // line
  {
    // to implement
  }
  else if ((ndime == 2 && nnode == 3) || // tri or
           (ndime == 2 && nnode == 4))   // quad
  {
    for (int i = 0; i < nnode; i++)
    {
      H(0, i) = gradn(i % nnode, 0);
      H(1, i) = gradn(i % nnode, 1);
    }

    for (int i = nnode; i < 2 * nnode; i++)
    {
      H(2, i) = gradn(i % nnode, 0);
      H(3, i) = gradn(i % nnode, 1);
    }
  }
  else if (ndime == 3)
  {
    for (int i = 0; i < nnode; i++)
    {
      H(0, i) = gradn(i % nnode, 0);
      H(1, i) = gradn(i % nnode, 1);
      H(2, i) = gradn(i % nnode, 2);
    }
    for (int i = nnode; i < 2 * nnode; i++)
    {
      H(3, i) = gradn(i % nnode, 0);
      H(4, i) = gradn(i % nnode, 1);
      H(5, i) = gradn(i % nnode, 2);
    }
    for (int i = 2 * nnode; i < 3 * nnode; i++)
    {
      H(6, i) = gradn(i % nnode, 0);
      H(7, i) = gradn(i % nnode, 1);
      H(8, i) = gradn(i % nnode, 2);
    }
  }
}

void UpdatedLagrangian::elem_resid(const int iel, const MxFE *fe,
                                   const Quadrature *qd, arma::vec &Re)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof / ndim;
  int nint = qd->get_num_ipoints();
  double detJxW, detF, press;

  arma::vec shape;
  arma::mat dshape, gradn(ndof, ndim), jacinv(ndim, ndim);
  arma::mat sigma_voigt(nvoig, 1), B(nvoig, ndof);
  arma::mat Ke(Re.size(), Re.size());
  arma::mat33 btens, *F, S, Sa, sigma;
  arma::vec3 qpt;

  std::vector<arma::vec3> elem_xe(nubf);
  std::vector<arma::vec3> elem_x0(nubf);
  get_element_x(iel, elem_xe);
  get_element_x0(iel, elem_x0);

  Mapping em = fe->get_mapping(iel, elem_xe);

  Ke.zeros();
  Re.zeros();
  press = 0;

  // Mean dilatation method
  if (material->is_incompressible())
    mean_dilatation(iel, fe, qd, press, Ke);

  // Loop over integration points
  for (int q = 0; q < qd->get_num_ipoints(); q++)
  {
    qpt = qd->get_point(q);
    fe->calc_shape_u(qpt, shape);
    fe->calc_deriv_shape_u(qpt, dshape);
    em.calc_jacobian(dshape, nubf);
    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    F = vecF[iel * nint + q];

    lcg_tensor(gradn, elem_x0, &detF, F, btens);
    calc_B_matrix(gradn, B);
    calc_B_matrix(gradn, B);

    // Compute stress and elasticity tensor
    MaterialData *md = new MaterialData(msh.get_element(iel), *F);
    if (material->is_incompressible())
    {

      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);

      // finite difference
      // im->calc_fd_stress(md, S);
      if (im->get_name() == "Guccione")
      {
        im->piola2_stress(md, S);
        if (md->get_marker() == 0)
          im->set_active_stress(iel, S, lc.load());
        im->map_to_global(md, S);
        im->push_forward(*F, S, sigma);
      }
      else
      {
        im->calc_fd_stress(iel, 
        md, S);
        // im->set_active_stress(md, S, lc.load());
        im->push_forward(*F, S, sigma);
        /*
                  arma::vec f = md->fiber();
                  cout << "Sa " << sigma <<endl;
                  im->calc_fd_stress(md, S);
                  im->push_forward(*F, S, sigma);
                  sigma = sigma + 1000.*((*F)*f)*((*F)*f).t();
                  cout << "Sigma_a " << sigma <<endl<<endl;
        */
        // exit(0);
        // im->set_active_stress(md, sigma, lc.load());
      }

      im->add_pressure(press, sigma);

      // analytic
      // im->deviatoric_stress(md, sigma);
      // im->add_pressure(press, sigma);
    }
    else
    {
      material->cauchy_stress(md, sigma);
    }

    // convert sigma stress tensor to format for matricial manipulation
    voigtvec(ndim, sigma, sigma_voigt);

    // internal nodal forces
    Re += (B.t() * sigma_voigt) * detJxW;

    // store stresses
    stressdb.tube(iel, q) = arma::vectorise(sigma, 1);
    straindb.tube(iel, q) = arma::vectorise(md->lagrangian_strain(), 1);

    // S = arma::det(*F) * arma::inv(*F) * sigma * arma::inv(*F).t();
    // stressdb.tube(iel, q) = arma::vectorise(S,1);

    delete md;
  }
}

/*
void UpdatedLagrangian::elem_stiff (const int iel, const MxFE * fe,
                                    const Quadrature * qd, arma::mat & Ke)
{
    int ndim = fe->get_ndim();
    int ndof = fe->get_ndofs_u();
    int nubf = ndof/ndim;
    int nint = qd->get_num_ipoints();
    double detJxW, detF, press;

    arma::vec shape;
    arma::mat dshape;
    arma::mat gradn(ndof,ndim);
    arma::mat jacinv(ndim,ndim);
    arma::mat sigma_voigt(nvoig,1);
    arma::mat sigma2(ndim*ndim, ndim*ndim);
    arma::mat B(nvoig,ndof), D(nvoig,nvoig), H(ndim*ndim,ndof);
    arma::mat33 btens, *F, S, sigma;
    arma::vec3 qpt;

    Tensor4 elastensorM, elastensor;

    std::vector<arma::vec3> elem_xe(nubf);
    std::vector<arma::vec3> elem_x0(nubf);
    get_element_x(iel, elem_xe);
    get_element_x0(iel, elem_x0);

    Mapping em = fe->get_mapping(iel, elem_xe);

    Ke.zeros();
    press = 0;

    // mean dilatation method
    if(material->is_incompressible())
    {
        mean_dilatation(iel, fe, qd, press, Ke);
    }

    // loop over integration points
    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
        qpt = qd->get_point(q);
        fe->calc_shape_u(qpt,shape);
        fe->calc_deriv_shape_u(qpt,dshape);
        em.calc_jacobian(dshape, nubf);
        detJxW = qd->get_weight(q) * em.get_det_jacobian();

        double jacdet = em.get_det_jacobian();
        if( !(jacdet>0) )
        {
            cerr << "Negative Jacobian determinant at element = " << iel << endl;
            exit(1);
        }

        jacinv = em.get_inv_jacobian();
        gradn  = (dshape * jacinv);

        F = vecF[iel*nint + q];

        lcg_tensor(gradn, elem_x0, &detF, F, btens);
        calc_B_matrix (gradn, B);
        calc_H_matrix (gradn, H);

        // Compute stress and elasticity tensor
        MaterialData * md = new MaterialData(msh.get_element(iel), *F);

        if( material->is_incompressible() )
        {
            IncompressibleMaterial * im;
            im = static_cast<IncompressibleMaterial*>(material);

            // stress - finite difference
            im->calc_fd_stress(md, S);
            im->push_forward(*F, S, sigma);
            im->add_pressure(press, sigma);

            // stress - analytical
            //im->deviatoric_stress(md, sigma);
            //im->add_pressure(press, sigma);



            // elasticity tensor - deviatoric part
            im->calc_fd_elastensor(md, elastensorM);
            im->push_forward(*F, elastensorM, elastensor);

            // elasticity tensor (deviatoric and volumetric) - analytical version
            //im->deviatoric_elastensor(md, elastensor);
            //im->sp_volumetric_elastensor(press, elastensor);

            // converts to Voigt notation
            elastensor.get_matrix(D);
        }
        else
        {
            material->cauchy_stress(md, sigma);
            material->sp_elastensor(md, D);
        }

        // create and convert sigma stress tensor to matrix
        sigma2 = arma::kron(arma::eye(ndim,ndim), sigma.submat(0,0,ndim-1,ndim-1));
        voigtvec(ndim, sigma, sigma_voigt);

        // material stiffness matrix - constitutive contribution (elmat_const)
        Ke += (B.t() * D * B) * detJxW;
        // geometrical stiffness matrix - stress contribution (elmat_sigma)
        Ke += (H.t() * sigma2 * H) * detJxW;

        delete md;
    }
}

*/


void UpdatedLagrangian::elem_stiff(const int iel, const MxFE *fe,
                                   const Quadrature *qd, arma::mat &Ke)
{
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof / ndim;
  int nint = qd->get_num_ipoints();
  double detJxW, detF, press;

  arma::vec shape;
  arma::mat dshape;
  arma::mat gradn(ndof, ndim);
  arma::mat jacinv(ndim, ndim);
  arma::mat sigma_voigt(nvoig, 1);
  arma::mat sigma2(ndim * ndim, ndim * ndim);
  arma::mat B(nvoig, ndof), D(nvoig, nvoig), H(ndim * ndim, ndof);
  arma::mat33 btens, *F, S, Sa, sigma;
  arma::vec3 qpt;

  Tensor4 elastensorM, elastensorA, elastensor, elastensorM_global;

  elastensorA.zero();
  elastensorM.zero();

  std::vector<arma::vec3> elem_xe(nubf);
  std::vector<arma::vec3> elem_x0(nubf);
  get_element_x(iel, elem_xe);
  get_element_x0(iel, elem_x0);

  Mapping em = fe->get_mapping(iel, elem_xe);

  Ke.zeros();
  press = 0;

  // mean dilatation method
  if (material->is_incompressible())
  {
    mean_dilatation(iel, fe, qd, press, Ke);
  }

  // cout << "Element: " << iel << endl;

  // loop over integration points
  for (int q = 0; q < qd->get_num_ipoints(); q++)
  {
    qpt = qd->get_point(q);
    fe->calc_shape_u(qpt, shape);
    fe->calc_deriv_shape_u(qpt, dshape);
    em.calc_jacobian(dshape, nubf);
    detJxW = qd->get_weight(q) * em.get_det_jacobian();

    double jacdet = em.get_det_jacobian();
    if (!(jacdet > 0))
    {
      cerr << "Negative Jacobian determinant at element = " << iel << endl;
      exit(1);
    }

    jacinv = em.get_inv_jacobian();
    gradn = (dshape * jacinv);

    F = vecF[iel * nint + q];

    lcg_tensor(gradn, elem_x0, &detF, F, btens);
    calc_B_matrix(gradn, B);
    calc_H_matrix(gradn, H);

    // Compute stress and elasticity tensor
    MaterialData *md = new MaterialData(msh.get_element(iel), *F);

    if (material->is_incompressible())
    {
      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);

      // stress - finite difference
      // im->calc_fd_stress(md, S);
      // guccione analittical
      if (im->get_name() == "Guccione")
      {
        im->piola2_stress(md, S);
        if (md->get_marker() == 0){
          im->set_active_stress(iel, S, lc.load());
        }
        im->map_to_global(md, S);
        im->push_forward(*F, S, sigma);
      }
      else
      {
        // cout << "ZZZ" << endl;
        im->calc_fd_stress(iel, md, S);
        // im->set_active_stress(md, S, lc.load());
        im->push_forward(*F, S, sigma);
        // im->set_active_stress(md, sigma, lc.load());
      }

      im->add_pressure(press, sigma);

      // stress - analytical
      // im->deviatoric_stress(md, sigma);
      // im->add_pressure(press, sigma);

      // elasticity tensor - deviatoric part
      // im->calc_fd_elastensor(md, elastensorM);
      // im->push_forward(*F, elastensorM, elastensor);

      // elasticity tensor (deviatoric and volumetric) - analytical version

      if (im->get_name() == "Guccione")
      {
        im->deviatoric_elastensor(md, elastensorM);
        im->map_elas_to_global(md, elastensorM, elastensorM_global);
        if (md->get_marker() == 0)
          im->active_stress_elastensor(iel, lc.get_nincs(), md->fiber(), elastensorM_global);
        im->push_forward(*F, elastensorM_global, elastensor);
      }
      else
      {
        im->calc_fd_elastensor(iel, md, elastensorM);
        // cout << "C: " << elastensorM(0,1,2,0) << endl;
        // im->active_stress_elastensor(lc.get_nincs(), md, elastensorM);
        // cout << "Ca: " << elastensorM(0,1,2,0) << endl << endl;
        im->push_forward(*F, elastensorM, elastensor);

        // cout << "XXX" << endl;
      }

      // analytical
      im->sp_volumetric_elastensor(press, elastensor);

      // converts to Voigt notation
      elastensor.get_matrix(D);
    }
    else
    {      
      material->cauchy_stress(md, sigma);
      material->sp_elastensor(md, D);
    }

    // create and convert sigma stress tensor to matrix
    sigma2 = arma::kron(arma::eye(ndim, ndim), sigma.submat(0, 0, ndim - 1, ndim - 1));
    voigtvec(ndim, sigma, sigma_voigt);

    // material stiffness matrix - constitutive contribution (elmat_const)
    Ke += (B.t() * D * B) * detJxW;
    // geometrical stiffness matrix - stress contribution (elmat_sigma)
    Ke += (H.t() * sigma2 * H) * detJxW;

    delete md;
  }
}

void UpdatedLagrangian::lcg_tensor(const arma::mat &gradn,
                                   const std::vector<arma::vec3> &x0,
                                   double *J, arma::mat33 *pF,
                                   arma::mat &btens)
{
  int ndim = gradn.n_cols;
  int nnod = gradn.n_rows;
  int nubf = nnod / ndim;
  arma::mat33 &F = *pF;
  arma::mat33 finvr;

  finvr.zeros();

  if (ndim == 2)
    finvr(2, 2) = 1.0;

  // Compute F^{-1}
  for (int id = 0; id < ndim; id++)
    for (int jd = 0; jd < ndim; jd++)
    {
      finvr(id, jd) = 0.0;
      for (int k = 0; k < nubf; k++)
        finvr(id, jd) += x0[k][id] * gradn(k, jd);
    }

  // Compute determinant
  *J = arma::det(finvr);

  if (*J == 0.0)
    *J = 1.0;
  *J = 1.0 / (*J);

  F = arma::inv(finvr);

  // compute b = F F^T
  btens = F * F.t();
}

void UpdatedLagrangian::mean_dilatation(const int iel, const MxFE *fe,
                                        const Quadrature *qd, double &press,
                                        arma::mat &Ke)
{
  int nint = qd->get_num_ipoints();
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof / ndim;
  double detJxW;
  arma::vec shape;
  arma::mat dshape;
  arma::mat gradn(ndof, ndim), jacinv(ndim, ndim), elacd(nubf, 4);
  arma::vec3 qpt;

  std::vector<arma::vec3> elem_xe(nubf);
  get_element_x(iel, elem_xe);

  Mapping em = fe->get_mapping(iel, elem_xe);

  elacd.zeros();

  if (material->is_incompressible())
  {
    IncompressibleMaterial *im = (IncompressibleMaterial *)material;

    double evol, Jbar, xkapp, kappa;
    evol = 0.0;

    // loops over Gauss points and finds the current volume
    // and average cartesian derivatives
    for (int q = 0; q < nint; q++)
    {
      qpt = qd->get_point(q);
      fe->calc_shape_u(qpt, shape);
      fe->calc_deriv_shape_u(qpt, dshape);
      em.calc_jacobian(dshape, nubf);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();
      jacinv = em.get_inv_jacobian();
      gradn = (dshape * jacinv);

      // compute average cartesian derivatives
      for (int i = 0; i < nubf; i++)
        for (int j = 0; j < ndim; j++)
          elacd(i, j) = elacd(i, j) + gradn(i, j) * detJxW;

      // compute element volume
      evol = evol + detJxW;
    }

    Jbar = evol / vol0(iel);

    if (!use_alg) // do NOT use AL
    {
      kappa = im->get_kappa();
      press = kappa * im->dUdJ(Jbar);
      xkapp = kappa * im->d2UdJJ(Jbar) * Jbar;
    }
    else // use augmented lagrangian
    {
      kappa = eps(iel);
      press = kappa * im->dUdJ(Jbar);
      xkapp = kappa * im->d2UdJJ(Jbar) * Jbar;
      press += eL(iel) * im->dhdp(Jbar);
      xkapp += eL(iel) * im->d2hdpp(Jbar);
    }

    ep(iel) = press;
    eJ(iel) = Jbar;

    // dilatational tangent stiffness component
    for (int i = 0; i < nubf; i++)
      for (int id = 0; id < ndim; id++)
        for (int j = 0; j < nubf; j++)
          for (int jd = 0; jd < ndim; jd++)
          {
            int ii, jj;
            double sum;
            ii = id * nubf + i;
            jj = jd * nubf + j;
            sum = xkapp * elacd(i, id) * elacd(j, jd) / evol;
            Ke(ii, jj) += sum;
          }
  }
}

void UpdatedLagrangian::pre_solve()
{

  cout << "Solving nonlinear problem (UL) using NonlinearSolver" << endl;
  // std::cout.precision(16);
  cout << "Initial inner volume: " << calc_volume(true) << endl;
  static bool start = true;

  if (start)
  {
    start = false;
    // body_forces();
    // assemble_traction();
    // assemble_const();
    int nelem = msh.get_n_elements();
    cout << "YYY " << nelem << endl;
    material->allocate_Ta(nelem);

    int iel = msh.get_n_elements();
    // arma::vec dta = arma::ones<arma::vec>(iel) * (material->get_Ta() / lc.get_nincs());

    arma::vec dta = (material->get_Ta() / lc.get_nincs());

    material->set_dTa(dta); 


    // PLACEHOLDER:
    // arma::vec ta = arma::ones<arma::vec>(iel) * (material->get_Ta());
    arma::vec ta = (material->get_Ta());
    material->set_Ta(ta);
  
    //material->set_dTa(material->get_Ta() / lc.get_nincs());
    //cout << "Active stress: " << material->get_Ta() << endl;
  }
  cout << "OKKKK" << endl;
}

void UpdatedLagrangian::solve()
{
  static NonlinearSolver *nls = NULL;
  static int cont = 1;

  int nits;
  int num_nz_prescribed = get_num_nz_prescribed();

  // create nonlinear solver
  if (nls == NULL)
  {
    nls = new NewtonLineSearch(this);
    nls->init();
  }

  // cout << "Solving nonlinear problem (UL) using NonlinearSolver" << endl;
  // std::cout.precision(16);
  // cout << "Initial volume: " << calc_volume(true) << endl;

  body_forces();
  assemble_traction();
  assemble_const();

  // use augmented lagrangian
  use_alg = false;

  if (output_step)
    output_vtk(cont, 0);

  // load increment loop
  while (lc.has_load())
  {
    lc.update();
    // augmented Lagrangian
    int al_iter = 0;

    double al_tol = 0.01;
    // double al_tol = 0.05;

    double al_norm = 0;
    double al_norm0 = 0;
    double al_reln;

    al_norm0 = arma::norm(eL, 2);

    do
    {
      al_iter++;
      evaluate_forces(((Newton *)nls)->residual());
      ////assemble_pressure();
      // cout<<"Cavity volume: "<<total_volume_cavity()<<endl;

      if (num_nz_prescribed > 0)
      {
        cout << "Prescribing non-zero displacements" << endl;
        prescribe_displacements();
        react = 0.0;
      }

      // call the nonlinear solver
      nits = nls->solve();

      // output deformation from this load increment
      if (output_step)
        output_vtk(cont, lc.increment());

      first_step = true;
      log << nits << " ";

      // augmented Lagrangian
      al_norm0 = al_norm;
      al_norm = al_augment(al_tol);
      al_reln = fabs((al_norm - al_norm0) / (al_norm));

      if (!use_alg)
        al_reln = 0;
      else
        cout << "   ALG step " << al_iter << " - lagnorm " << al_reln << endl;

    } while (al_reln > al_tol);

    // update penalty parameter
    al_update(lc.increment());
  }

  // cout << "Cavity volume: " << total_volume_cavity() << endl;

  // prepare to leave
  cont += 1;     // counter for output
  fext0 += fext; // save nodal loads
  log << calc_volume() << endl;
  
  // std::cout.precision(16);

  cout << "End inner volume: " << calc_volume() << endl;
  nls->timer.summary();

  cout << " Done" << endl;
}
