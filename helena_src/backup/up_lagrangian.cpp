// **************** BACKUP Routines for Coupled electromechanics **************

void UpdatedLagrangian::calc_deformation_grad(const arma::mat & gradn,
    const std::vector<arma::vec3> & x0, Real * detf, arma::mat & ftens)
{
  int ndim = gradn.n_cols;
  int nnod = gradn.n_rows;
  int nubf = nnod/ndim;
  arma::mat finvr(3,3);//, ftens(3,3);
  finvr.zeros();

  if (ndim == 2) finvr(2,2) = 1.0;

  // Computes F^{-1}
  for(int id=0; id<ndim; id++)
    for(int jd=0; jd<ndim; jd++)
    {
      finvr(id,jd)=0.0;
      for(int k=0; k<nubf; k++)
        finvr(id,jd) += x0[k][id] * gradn(k,jd);
    }

  // Computes determinant
  *detf = finvr(0, 0) * finvr(1, 1) * finvr(2, 2)
        + finvr(0, 1) * finvr(1, 2) * finvr(2, 0)
        + finvr(0, 2) * finvr(1, 0) * finvr(2, 1)
        - finvr(0, 2) * finvr(1, 1) * finvr(2, 0)
        - finvr(0, 1) * finvr(1, 0) * finvr(2, 2)
        - finvr(0, 0) * finvr(1, 2) * finvr(2, 1);

  if (*detf == 0.0) *detf = 1.0;
  *detf = 1.0/(*detf);

  // Computes F
  ftens(0,0) = (finvr(1,1)*finvr(2,2) - finvr(1,2)*finvr(2,1)) * (*detf);
  ftens(0,1) = (finvr(0,2)*finvr(2,1) - finvr(0,1)*finvr(2,2)) * (*detf);
  ftens(0,2) = (finvr(0,1)*finvr(1,2) - finvr(0,2)*finvr(1,1)) * (*detf);

  ftens(1,0) = (finvr(1,2)*finvr(2,0) - finvr(1,0)*finvr(2,2)) * (*detf);
  ftens(1,1) = (finvr(0,0)*finvr(2,2) - finvr(0,2)*finvr(2,0)) * (*detf);
  ftens(1,2) = (finvr(0,2)*finvr(1,0) - finvr(0,0)*finvr(1,2)) * (*detf);

  ftens(2,0) = (finvr(1,0)*finvr(2,1) - finvr(1,1)*finvr(2,0)) * (*detf);
  ftens(2,1) = (finvr(0,1)*finvr(2,0) - finvr(0,0)*finvr(2,1)) * (*detf);
  ftens(2,2) = (finvr(0,0)*finvr(1,1) - finvr(0,1)*finvr(1,0)) * (*detf);
}

void UpdatedLagrangian::output_F_tensor(std::vector<arma::mat*> & v_ftens)
{
  std::vector<int> dnums;
  MixedFiniteElement * fe = fespace.createFE();

  v_ftens.clear();

  for(int i=0; i<msh.get_n_elements(); i++)
  {
    // BEGIN CALC ELMAT VEC
    int ndim = fe->get_ndim();
    int ndof = fe->get_ndofs_u();
    int nubf = ndof/ndim;

    Real detF;
    arma::vec shape;
    arma::mat dshape;
    arma::mat gradn(ndof,ndim);
    arma::mat jacinv(ndim,ndim);
    arma::mat * F = new arma::mat(3,3);

    arma::vec3 pt, qpt;

    std::vector<arma::vec3> xe(ndof);
    std::vector<arma::vec3> x0(ndof);
    this->get_element_x(i, xe);
    this->get_element_x0(i, x0);
    Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type());
    Mapping em = fe->get_mapping(i, xe);

    fespace.get_element_dofs_u(i,dnums);

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      //qd.get_point(q, qpt);
      qpt = qd->get_point(q);
      fe->calc_shape_u(qpt,shape);
      fe->calc_deriv_shape_u(qpt,dshape);
      em.calc_jacobian(dshape, nubf);
      jacinv = em.get_inv_jacobian();
      gradn = (dshape * jacinv);

      calc_deformation_grad(gradn, x0, &detF, *F);

      v_ftens.push_back(F);
    }
    delete qd;
  }
  delete fe;
}

void UpdatedLagrangian::calc_right_cauchy_green(const arma::mat & gradn,
    const std::vector<arma::vec3> & x0, Real * detf, arma::mat & ctens)
{
  int ndim = gradn.n_cols;
  int nnod = gradn.n_rows;
  arma::mat finvr(3,3), ftens(3,3);
  finvr.zeros();

  if (ndim == 2) finvr(2,2) = 1.0;

  // Computes F^{-1}
  for(int id=0; id<ndim; id++)
    for(int jd=0; jd<ndim; jd++)
    {
      finvr(id,jd)=0.0;
      for(int k=0; k<nnod; k++)
        finvr(id,jd) += x0[k][id] * gradn(k,jd);
    }

  // Computes determinant
  *detf = finvr(0, 0) * finvr(1, 1) * finvr(2, 2)
        + finvr(0, 1) * finvr(1, 2) * finvr(2, 0)
        + finvr(0, 2) * finvr(1, 0) * finvr(2, 1)
        - finvr(0, 2) * finvr(1, 1) * finvr(2, 0)
        - finvr(0, 1) * finvr(1, 0) * finvr(2, 2)
        - finvr(0, 0) * finvr(1, 2) * finvr(2, 1);

  if (*detf == 0.0) *detf = 1.0;
  *detf = 1.0/(*detf);

  // Computes F
  ftens(0,0) = (finvr(1,1)*finvr(2,2) - finvr(1,2)*finvr(2,1)) * (*detf);
  ftens(0,1) = (finvr(0,2)*finvr(2,1) - finvr(0,1)*finvr(2,2)) * (*detf);
  ftens(0,2) = (finvr(0,1)*finvr(1,2) - finvr(0,2)*finvr(1,1)) * (*detf);

  ftens(1,0) = (finvr(1,2)*finvr(2,0) - finvr(1,0)*finvr(2,2)) * (*detf);
  ftens(1,1) = (finvr(0,0)*finvr(2,2) - finvr(0,2)*finvr(2,0)) * (*detf);
  ftens(1,2) = (finvr(0,2)*finvr(1,0) - finvr(0,0)*finvr(1,2)) * (*detf);

  ftens(2,0) = (finvr(1,0)*finvr(2,1) - finvr(1,1)*finvr(2,0)) * (*detf);
  ftens(2,1) = (finvr(0,1)*finvr(2,0) - finvr(0,0)*finvr(2,1)) * (*detf);
  ftens(2,2) = (finvr(0,0)*finvr(1,1) - finvr(0,1)*finvr(1,0)) * (*detf);

  // Obtains C = F^T F
  ctens.zeros();
  //ctens = ftens.t() * ftens;

  if (ndim == 2) ctens(2, 2) = 1.0;

  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      ctens(i,j) = 0.0;
      for(int k=0; k<ndim; k++){
        ctens(i,j) += ftens (k,i) * ftens(k,j);
      }
    }
  }
  
}

void UpdatedLagrangian::output_C_tensor(std::vector<arma::mat*> & v_ctens)
{
  std::vector<int> dnums;
  MixedFiniteElement * fe = fespace.createFE();

  v_ctens.clear();

  for(int i=0; i<msh.get_n_elements(); i++)
  {
    // BEGIN CALC ELMAT VEC
    int ndim = fe->get_ndim();
    int ndof = fe->get_ndofs_u();
    int nubf = ndof/ndim;

    Real detF;
    arma::vec shape;
    arma::mat dshape;
    arma::mat gradn(ndof,ndim);
    arma::mat jacinv(ndim,ndim);
    arma::mat * Ctens = new arma::mat(3,3);

    arma::vec3 pt, qpt;

    std::vector<arma::vec3> xe(ndof);
    std::vector<arma::vec3> x0(ndof);
    this->get_element_x(i, xe);
    this->get_element_x0(i, x0);
    Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type());
    Mapping em = fe->get_mapping(i, xe);

    fespace.get_element_dofs_u(i,dnums);

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      //qd.get_point(q, qpt);
      qpt = qd->get_point(q);

      fe->calc_shape_u(qpt,shape);
      fe->calc_deriv_shape_u(qpt,dshape);
      em.calc_jacobian(dshape, nubf);
      jacinv = em.get_inv_jacobian();
      gradn = (dshape * jacinv);

      calc_right_cauchy_green(gradn, x0, &detF, *Ctens);

      v_ctens.push_back( Ctens );
    }
    delete &qd;
  }
  delete fe;
}
*/

// ****************************************************************************


// ************** BACKUP Calcelmatvec in indicial form ************************
/*
calc_elmatvec(){

   // Computations in indicial form

    // Internal equivalent nodal forces
    //for(int i=0; i<nubf; i++)
    //  for(int id=0; id<ndim; id++)
    //    for(int jd=0; jd<ndim; jd++)
    //      elvec(i+id*nubf) += sigma(id,jd) * gradn(i,jd) * detJxW;

    // Constitutive part 
    // *** Need to check
    for(int i=0; i<nubf; i++)
      for(int id=0; id<ndim; id++)
      {
        for(int j=0; j<nubf; j++)
          for(int jd=0; jd<ndim; jd++)
          {
            Real sum = 0.0;

            for(int kd=0; kd<ndim; kd++)
              for(int ld=0; ld<ndim; ld++)
                sum += gradn(i,kd) * ctens(id,jd,kd,ld) * gradn(j,ld) * detJxW;

            elmat_const(i+id*nubf,j+id*nubf) += sum;
          }
      }

    // Stress component using indicial notation
    for(int i=0; i<nubf; i++)
      for(int j=0; j<nubf; j++)
      {
        Real sum = 0.0;
        for(int kd=0; kd<ndim; kd++)
          for(int ld=0; ld<ndim; ld++)
            sum += gradn(i,kd) * sigma(kd,ld) * gradn(j,ld) * detJxW;
	    //elmat_sigma(i+kd*nubf,j+kd*nubf) 
	    //+=gradn(i,kd) * sigma(kd,ld) * gradn(j,ld) * detJxW;

        for(int id=0; id<ndim; id++)
          elmat_sigma(i+id*nubf,j+id*nubf) += sum;
      }
}
*/


// ****************************************************************************


/*
// Calculate Von Mises Stress of Mixed FE Variable
void UpdatedLagrangian::von_mises(const arma::vec & u, arma::vec & sigma)
{
  int nels = msh.get_n_elements();
  vector<int> dnums,sdnums;

  // Scalar field associated with Von Mises Stress in each element
  sigma.resize(nels);
  sigma = 0.0;

  for(int i=0;i<nels; i++)
  {
    MixedFiniteElement * fe = fespace.createFE();
    Mapping em = msh.get_element_mapping(i);
    arma::vec elvecu, elsigma;

    fespace.get_element_dofs_u(i,dnums);

    this->get_elvec(dnums, udisp, elvecu);
    this->calc_elstress(em, fe, elvecu, msh.get_element_centroid(i), elsigma);

    Real v = 0.3;
    Real sigma11 = elsigma(0);
    Real sigma22 = elsigma(1);
    Real sigma12 = elsigma(2);
    Real sigma33 = v*(sigma11 + sigma22);

    Real s1 = sigma11*sigma11;
    Real s2 = sigma22*sigma22;
    Real s3 = sigma33*sigma33;
    Real s4 = -(sigma11*sigma22);
    Real s5 = -(sigma22*sigma33);
    Real s6 = -(sigma33*sigma11);
    Real s7 = (3.0*sigma12*sigma12);

    sigma(i) = sqrt(s1 + s2 + s3 + s4 + s5 + s6 + s7);

    delete fe;
  }
}

void UpdatedLagrangian::calc_elstress(Mapping & em, const MxFE * fe,
				      const arma::vec & u, const arma::vec3 & X, 
				      arma::vec & sigma)
{
  int index;
  int numc = num_components;
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  arma::mat dshape;
  arma::mat gradT_jac(ndof,ndim);
  arma::mat jacinv(ndim,ndim);
  arma::mat B(numc,ndof), D(numc,numc);

  sigma.resize(numc);

  index = em.get_index();

  fe->calc_deriv_shape_u(X,dshape);
  em.calc_jacobian(dshape, nubf);
  jacinv = em.get_inv_jacobian();

  gradT_jac = dshape * jacinv;

  // Material properties
  Real xmu=3.0;
  Real xla=2.0;

  calc_B_matrix (gradT_jac, B);
  calc_D_matrix (X, index, D, xmu, xla);

  sigma = D * B * u;
}

