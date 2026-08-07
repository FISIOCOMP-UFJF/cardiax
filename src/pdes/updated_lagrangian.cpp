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
        // adaptative kappa
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
    md->set_active_stress(0.0); 

    if (material->is_incompressible())
    {

      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);
      im->calc_fd_elastensor(md, elastensorM);
      im->push_forward(F, elastensorM, elastensor);
      
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
  std::vector<int> pnums(nubf);
  get_element_x(iel, elem_xe);
  get_element_x0(iel, elem_x0);

  Mapping em = fe->get_mapping(iel, elem_xe);

  Ke.zeros();
  Re.zeros();
  press = 0;

  arma::vec Ta_e(nubf);
  msh.get_element_pt_nums(iel, pnums);
  const arma::vec &Ta = material->get_Ta(); 
  for(int i=0;i<nubf;i++)
    Ta_e(i) = Ta(pnums[i]);


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

    double Ta_ip  =0.0; 
    for (int j = 0; j < nubf; ++j)
    {
        Ta_ip += shape(j) * Ta_e(j);
    }

    // compute stress and elasticity tensor
    MaterialData *md = new MaterialData(msh.get_element(iel), *F);
    if(md->get_marker() == 0)
      md->set_active_stress(Ta_ip*lc.load());
    else
      md->set_active_stress(0.0);

    if (material->is_incompressible())
    {
      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);

      im->calc_fd_stress(iel, md, S);
      im->push_forward(*F, S, sigma);
      

      im->add_pressure(press, sigma);
      im->calc_fd_elastensor(md, elastensorM);
      im->push_forward(*F, elastensorM, elastensor);

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
  std::vector<int> pnums(nubf);

  get_element_x(iel, elem_xe);
  get_element_x0(iel, elem_x0);

  Mapping em = fe->get_mapping(iel, elem_xe);

  Ke.zeros();
  Re.zeros();
  press = 0;

  arma::vec Ta_e(nubf);
  msh.get_element_pt_nums(iel, pnums);
  const arma::vec &Ta = material->get_Ta(); 
  for(int i=0;i<nubf;i++)
    Ta_e(i) = Ta(pnums[i]);

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

    
    double Ta_ip  =0.0; 
    for (int j = 0; j < nubf; ++j)
    {
        Ta_ip += shape(j) * Ta_e(j);
    }


    // Compute stress and elasticity tensor
    MaterialData *md = new MaterialData(msh.get_element(iel), *F);
    if(md->get_marker() == 0)
      md->set_active_stress(Ta_ip*lc.load());
    else
      md->set_active_stress(0.0);

    if (material->is_incompressible())
    {
      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);
      
      im->calc_fd_stress(iel, md, S);
      im->push_forward(*F, S, sigma);

      im->add_pressure(press, sigma);
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

    delete md;
  }
}

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
  std::vector<int> pnums(nubf);

  get_element_x(iel, elem_xe);
  get_element_x0(iel, elem_x0);

  Mapping em = fe->get_mapping(iel, elem_xe);

  Ke.zeros();
  press = 0;

  arma::vec Ta_e(nubf);
  msh.get_element_pt_nums(iel, pnums);
  const arma::vec &Ta = material->get_Ta(); 
  for(int i=0;i<nubf;i++){
    Ta_e(i) = Ta(pnums[i]);
  }
  // mean dilatation method
  if (material->is_incompressible())
  {
    mean_dilatation(iel, fe, qd, press, Ke);
  }

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
   
    double Ta_ip  =0.0; 
    for (int j = 0; j < nubf; ++j)
    {
        Ta_ip += shape(j) * Ta_e(j);
    }

    // Compute stress and elasticity tensor
    MaterialData *md = new MaterialData(msh.get_element(iel), *F);
    if(md->get_marker() == 0)
      md->set_active_stress(Ta_ip*lc.load());
    else
      md->set_active_stress(0.0);
    
    if (material->is_incompressible())
    {
      IncompressibleMaterial *im;
      im = static_cast<IncompressibleMaterial *>(material);

      im->calc_fd_stress(iel, md, S);
      im->push_forward(*F, S, sigma);

      im->add_pressure(press, sigma);
      im->calc_fd_elastensor(md, elastensorM);
      im->push_forward(*F, elastensorM, elastensor);

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
  cout << "Initial inner volume: " << calc_volume(true) << endl;
  static bool start = true;

  if (start)
  {
    start = false;
    
    int nelem = msh.get_n_elements();
    cout << "Number of Elements " << nelem << endl;

    arma::vec dta = (material->get_Ta() / lc.get_nincs());
    material->set_dTa(dta);  //TODO; probably this is unecessary
  }
  cout << "End of UL pre solve" << endl;
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
  body_forces();
  assemble_traction();
  assemble_const();

  // use augmented lagrangian
  use_alg = false;

  // per-increment PV history: only meaningful for a monotonic load ramp
  bool first_call = !pv_history_is_open();
  if (pv_record_increments)
  {
    open_pv_history(this->basename);
    if (first_call)
      record_pv_history(0);
  }

  if (output_step)
    output_vtk(cont, 0);

  // load increment loop
  while (lc.has_load())
  {
    lc.update();
    // augmented Lagrangian
    int al_iter = 0;

    double al_tol = 0.01;
    double al_norm = 0;
    double al_norm0 = 0;
    double al_reln;

    al_norm0 = arma::norm(eL, 2);

    do
    {
      al_iter++;
      evaluate_forces(((Newton *)nls)->residual());
      
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

    // record pressure and volume of both cavities for this increment
    if (pv_record_increments)
      record_pv_history(lc.increment());
  }

  // NOTE: the PV history file is intentionally left open here. solve() may
  // be called once per cardiac-cycle step by an outer driver, and closing
  // it now would truncate the history on the next call. It is closed in
  // run() (standalone path) or by the driver via close_pv_history().

  // prepare to leave
  cont += 1;     // counter for output
  fext0 += fext; // save nodal loads
  log << calc_volume() << endl;
  
  cout << "End inner volume: " << calc_volume() << endl;
  cout << "End LV cavity volume: " << volume_LV() << endl;
  cout << "End RV cavity volume: " << volume_RV() << endl;
  //nls->timer.summary();

  cout << " Done" << endl;
}
