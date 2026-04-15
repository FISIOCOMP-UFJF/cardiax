#include "total_lagrangian.hpp"

TotalLagrangian::TotalLagrangian() : NonlinearElasticity()
{
  // do nothing
}

TotalLagrangian::~TotalLagrangian()
{
  for(uint i=0; i<fedata.size(); i++) 
    delete fedata[i];
}

void TotalLagrangian::assemble_const()
{
  int n  = msh.get_n_dim() * msh.get_nen();
  int ne = msh.get_n_elements();
  int * pidx;
  arma::mat Ke(n,n);
  std::vector<int> dnums;
  
  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type());

  for(int i=0; i<ne; i++)
  {
    calc_elmat_const(i, fe, qd, Ke);
    fespace.get_element_dofs_u (i,dnums);
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

void TotalLagrangian::assemble_stiff()
{
  int n  = msh.get_n_dim() * msh.get_nen();
  int ne = msh.get_n_elements();
  int * pidx;
  arma::vec Re(n);
  arma::mat Ke(n,n);
  std::vector<int> dnums; 

  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());
  
  for(int i=0; i<ne; i++)
  {
    calc_elmatvec (i, fe, qd, Ke, Re);
    fespace.get_element_dofs_u (i, dnums);

    // Fast assembling
    pidx = &dnums[0];
    Ke = Ke.t();
    K.add(n, n, pidx, pidx, Ke.memptr());    

    for(int k=0; k<n; k++)
    {
      if (ldgof[dnums[k]]) 
        r.add(dnums[k], -Re(k)); // add -R
      else
        react.add(dnums[k], Re(k));
    }    
  }

  K.assemble();
  r.assemble();
  apply_boundary(K);

  delete qd;
  delete fe;
}

void TotalLagrangian::assemble_active(const arma::vec & is,
                                      std::vector<arma::mat33*> & vstrs,
                                      std::vector<arma::mat33*> & vfibcoords)
{
  MixedFiniteElement * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(0, fe->get_type());

  int n  = msh.get_n_dim() * msh.get_nen();
  int ne = msh.get_n_elements();
  int ndim  = fe->get_ndim();
  int ndof  = fe->get_ndofs_u();
  int nnode = fe->get_nnode();
  int nint  = qd->get_num_ipoints();

  arma::vec tmp_fext(fespace.get_ndofs());
  arma::vec elvec(n);
  std::vector<int> dnums;

  tmp_fext.zeros();

  for(int i=0; i<ne; i++)
  {
    double detJxW, J;
    arma::vec shape;
    arma::vec3 pt, qpt;
    arma::mat dshape, gradn(ndof,ndim), jacinv(ndim,ndim);
    arma::mat33 *F, C, S;
    arma::mat Sv(nvoig,1);
    arma::mat Bl(nvoig,ndof);
    std::vector<arma::vec3> x0(ndof);    
    get_element_x0(i, x0);
      
    Mapping em = fe->get_mapping(i, x0);

    elvec.zeros();  
    fespace.get_element_dofs_u(i,dnums);
//cout << "passou 2\n";

    for(int q=0; q<nint; q++)
    {
      qpt = qd->get_point(q);
      fe->calc_shape_u(qpt, shape);
      fe->calc_deriv_shape_u(qpt, dshape);
      em.calc_jacobian(dshape, nnode);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();
      jacinv = em.get_inv_jacobian();      
      gradn = (dshape * jacinv);

      // deformation gradient F and right cauchy green C=F^T F
      F = vecF[i*nint + q];
      rcg_tensor(gradn, x0, J, F, C);
      Bl_matrix (gradn, F, Bl);            
      // interpolation of active stress from nodal to integration points
      arma::mat33 T;                        // active stress tensor 
      arma::mat33 Tl(arma::fill::zeros);    // active local stress tensor
      arma::mat33 & M  = *(vfibcoords[i]);  // change of basis matrix
      arma::mat33 & Ta = *(vstrs[i]);       // anisotropic Cauchy active stress
      arma::mat33 Fi   = arma::inv(*F);     // inverse of F
      arma::vec3 f0 = M.col(0);
      arma::vec3 f  = (*F) * f0;
      double lambda2 = pow(arma::norm(f,2),2);

      Ta.zeros();
      for(int j=0; j<nnode; j++)
      {	
        // int node = ( dnums[j] - 2) / 3 ; // global node number
        int node = dnums[j];

//cout << "j: " << j << " node: " << node <<endl;

        Tl(0,0) = is(node);  // local active tension in fiber direction	
        //T = M * Tl * M.t();  // global active tension
        //T = (Tl(0,0)/lambda2) * (f0 * f0.t());
        T = is(node)*(f * f.t());
        
	      // interpolate from node to integration point
	      for(int ii=0; ii<3; ii++)
          for(int jj=0; jj<3; jj++)
	          Ta(ii,jj) += shape(j) * T(ii,jj); //Ta(1,1) += shape(j) * is(k);
      }
      cout << Ta << endl;

      // Total Lagrangian formulation, thus
      // since Ta is Cauchy stress, we need to convert to PK2
      //S = Ta;
      S = J * (Fi * Ta * Fi.t());
      voigtvec(ndim, S, Sv);
      
      // indicial notation
      //for(int i=0; i<nnode; i++)
	    //  for(int id=0; id<ndim; id++)
	    //    for(int jd=0; jd<ndim; jd++)
	    //      elvec(i+id*nnode) += S(id,jd) * gradn(i,jd) * detJxW;     
      
      // matrix notation
      elvec += (Bl.t() * Sv) * detJxW;
    }

    // assembly initial stress load vector
    for(int k=0; k<n; k++)
      if (ldgof[dnums[k]])
        tmp_fext(dnums[k]) -= elvec(k);
  }

  // copy
  fext = tmp_fext;

  delete qd;
  delete fe; 
}

void TotalLagrangian::calc_elmat_const (const int iel, const MxFE * fe, 
                                        const Quadrature * qd,
                                        arma::mat & elmat)
{
  int ndim  = fe->get_ndim();
  int ndof  = fe->get_ndofs_u();
  int nnode = ndof/ndim;
  int nint  = qd->get_num_ipoints();
  double detJxW, detF;
  arma::mat Bl(nvoig,ndof), D(nvoig,nvoig);
  arma::mat33 C, *F;

  std::vector<arma::vec3> elem_xe(nnode);
  get_element_x (iel, elem_xe);
  
  //std::vector<arma::vec3> elem_x0(nnode);
  //get_element_x0(iel, elem_x0);

  // map to ref element
  //Mapping em = fe->get_mapping(iel, elem_x0); 

  elmat.zeros();
  for(int q=0; q<nint; q++)
  {
    //qpt = qd->get_point(q);
    //fe->calc_deriv_shape_u(qpt,dshape);
    //em.calc_jacobian(dshape, nnode);
    //detJxW = qd->get_weight(q) * em.get_det_jacobian();
    //jacinv = em.get_inv_jacobian();
    //gradn  = (dshape * jacinv);
    
    // Using MxFEData
    const arma::mat & gradn = fedata[iel]->get_gradn(q);
    double d = fedata[iel]->get_jacobian_det(q);
    detJxW = qd->get_weight(q) * d;
    F = vecF[iel*nint + q];

    rcg_tensor(gradn, elem_xe, detF, F, C);
    Bl_matrix(gradn, F, Bl);

    // compute elasticity tensor using FD
    MaterialData * md = new MaterialData(msh.get_element(iel), *F);
    material->calc_fd_elastensor(iel, md, D);
    
    // assemble element matrix
    elmat += (Bl.t() * D * Bl) * detJxW;

    delete md;
  }
}

void TotalLagrangian::calc_elmatvec (const int iel, const MxFE * fe,
                                     const Quadrature * qd,
                                     arma::mat & Ke, arma::vec & Re)
{
  int ndim  = fe->get_ndim();
  int ndof  = fe->get_ndofs_u();
  int nnode = fe->get_nnode();
  int nint, ndim2 = ndim*ndim;
  double detJxW, detF;

  //arma::vec shape;
  //arma::mat dshape;
  //arma::mat gradn(ndof,ndim);
  //arma::mat jacinv(ndim,ndim);

  arma::mat Sb(ndim2, ndim2); // block
  arma::mat Sv(nvoig,1);      // voigt
  arma::mat Bl(nvoig,ndof), Bnl(ndim2,ndof), D(nvoig,nvoig);
  arma::mat33 *F, C, S, T;
  arma::vec3 qpt;
  static const arma::mat I = arma::eye(ndim,ndim);

  std::vector<arma::vec3> elem_xe(nnode);
  get_element_x (iel, elem_xe);

  //Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type());
  //Mapping em = fe->get_mapping(iel, elem_x0);

  nint = qd->get_num_ipoints();

  Ke.zeros();
  Re.zeros();

  for(int q=0; q<nint; q++)
  {
    qpt = qd->get_point(q);
    
    // old scheme
    //fe->calc_shape_u(qpt,shape);
    //fe->calc_deriv_shape_u(qpt,dshape);
    //em.calc_jacobian(dshape, nnode);
    //detJxW = qd->get_weight(q) * em.get_det_jacobian();
    //jacinv = em.get_inv_jacobian();

    // Compute dN/dX (with respect to X reference coords )
    //gradn = (dshape * jacinv);
    
    // Using MxFEData
    const arma::mat & gradn = fedata[iel]->get_gradn(q);
    double d = fedata[iel]->get_jacobian_det(q);
    detJxW = qd->get_weight(q) * d;

    F = vecF[iel*nint + q];
    
    rcg_tensor (gradn, elem_xe, detF, F, C);
    Bl_matrix  (gradn, F, Bl);
    Bnl_matrix (gradn, Bnl);

    // Compute stress and stress derivatives using FD
    MaterialData * md = new MaterialData(msh.get_element(iel), *F);
    material->calc_fd_stress(iel, md, S);
    material->calc_fd_elastensor(iel, md, D);

    // stress
    voigtvec(ndim, S, Sv); 
    Sb = arma::kron(I, S.submat(0, 0, ndim-1, ndim-1));

    // computations in matricial form
    Re += (Bl.t()  * Sv)       * detJxW;    
    Ke += (Bl.t()  * D  * Bl)  * detJxW;
    Ke += (Bnl.t() * Sb * Bnl) * detJxW;

    // store Cauchy stress
    T = (1.0/detF) * (*F) * S * (*F).t();
    stressdb.tube(iel, q) = arma::vectorise(T,1);
    //stressdb.tube(iel, q) = arma::vectorise(S,1);
    straindb.tube(iel, q) = arma::vectorise(md->lagrangian_strain(),1);
    delete md;
  }
}

void TotalLagrangian::elem_stiff(const int iel, const MxFE * fe,
                                 const Quadrature * qd, arma::mat & Ke)
{
  int ndim  = fe->get_ndim();
  int ndof  = fe->get_ndofs_u();
  int nnode = fe->get_nnode();
  int nint  = qd->get_num_ipoints();
  int ndim2 = ndim*ndim;
  double detJxW, detF;

  arma::mat Sb(ndim2, ndim2); // block
  arma::mat Bl(nvoig,ndof), Bnl(ndim2,ndof), D(nvoig,nvoig);
  arma::mat33 *F, C, S, T;
  static const arma::mat I = arma::eye(ndim,ndim);

  std::vector<arma::vec3> elem_xe(nnode);
  get_element_x (iel, elem_xe);

  Ke.zeros();

  for(int q=0; q<nint; q++)
  {
    // using MxFEData
    const arma::mat & gradn = fedata[iel]->get_gradn(q);
    double d = fedata[iel]->get_jacobian_det(q);
    detJxW = qd->get_weight(q) * d;

    F = vecF[iel*nint + q];
    
    rcg_tensor (gradn, elem_xe, detF, F, C);
    Bl_matrix  (gradn, F, Bl);
    Bnl_matrix (gradn, Bnl);

    //double lamb=1.05;
    //(*F).fill(0.0);
    //(*F)(0,0)=lamb;
    //(*F)(1,1)=1.0/sqrt(lamb);
    //(*F)(2,2)=1.0/sqrt(lamb);
    //cout << (*F) << endl;

    // Compute stress and stress derivatives using FD
    MaterialData * md = new MaterialData(msh.get_element(iel), *F);
    material->calc_fd_stress(iel, md, S);
    material->calc_fd_elastensor(iel, md, D);
    Sb = arma::kron(I, S.submat(0, 0, ndim-1, ndim-1));

    // computations in matricial form
    Ke += (Bl.t()  * D  * Bl)  * detJxW;
    Ke += (Bnl.t() * Sb * Bnl) * detJxW;

    // store Cauchy stress
    T = (1.0/detF) * (*F) * S * (*F).t();
    stressdb.tube(iel, q) = arma::vectorise(T,1);
    //stressdb.tube(iel, q) = arma::vectorise(S,1);
    straindb.tube(iel, q) = arma::vectorise(md->lagrangian_strain(),1);

    delete md;
  }

}

void TotalLagrangian::elem_resid (const int iel, const MxFE * fe, 
                                  const Quadrature * qd, arma::vec & Re)
{
  int ndim  = fe->get_ndim();
  int ndof  = fe->get_ndofs_u();
  int nnode = fe->get_nnode();
  int nint  = qd->get_num_ipoints();
  double detJxW, detF; 
  arma::mat Sv(nvoig,1);
  arma::mat Bl(nvoig,ndof);
  arma::mat33 *F, C, S;

  std::vector<arma::vec3> elem_xe(nnode);
  get_element_x (iel, elem_xe);

  Re.zeros();   

  for(int q=0; q<nint; q++)
  {   
    // using MxFEData
    const arma::mat & gradn = fedata[iel]->get_gradn(q);
    double d = fedata[iel]->get_jacobian_det(q);
    detJxW = qd->get_weight(q) * d;

    F = vecF[iel*nint + q];
    
    // kinematics and FE matrices
    rcg_tensor (gradn, elem_xe, detF, F, C);
    Bl_matrix  (gradn, F, Bl);

    // compute stress using FD
    MaterialData * md = new MaterialData(msh.get_element(iel), *F);
    material->calc_fd_stress(iel, md, S);
    voigtvec(ndim, S, Sv); 
    delete md;
   
    // compute residual
    Re += (Bl.t() * Sv) * detJxW;
  }
}

void TotalLagrangian::Bl_matrix (const arma::mat & gradn, const arma::mat * F, 
                                 arma::mat & Bl) const
{
  Bl.zeros();

  int ndof  = gradn.n_rows;    // ndof
  int ndim  = msh.get_n_dim(); // ndim
  int nnode = ndof/ndim;       // nnode

  // cout << "GRADN " << gradn.n_rows << " " << gradn.n_cols << endl;
  // cout << gradn << endl << endl;

  if (ndim == 2)
  {
    // row 1
    for(int i=0; i<nnode; i++) 
      Bl(0,i) = (*F)(0,0) * gradn(i,0);
    for(int i=nnode; i<2*nnode; i++) 
      Bl(0,i) = (*F)(1,0) * gradn(i,0);

    // row 2
    for(int i=0; i<nnode; i++) 
      Bl(1,i) = (*F)(0,1) * gradn(i,1);
    for(int i=nnode; i<2*nnode; i++) 
      Bl(1,i) = (*F)(1,1) * gradn(i,1);

    // row 3
    for(int i=0; i<nnode; i++)
      Bl(2,i) = (*F)(0,0) * gradn(i,1) + (*F)(0,1) * gradn(i,0);
    for(int i=nnode; i<2*nnode; i++)
      Bl(2,i) = (*F)(1,0) * gradn(i,1) + (*F)(1,1) * gradn(i,0);
  }
  else if (ndim == 3)
  {
    // row 1
    for(int i=0; i<nnode; i++)
      Bl(0,i) = (*F)(0,0) * gradn(i,0);
    for(int i=nnode; i<2*nnode; i++)
      Bl(0,i) = (*F)(1,0) * gradn(i,0);
    for(int i=2*nnode; i<3*nnode; i++)
      Bl(0,i) = (*F)(2,0) * gradn(i,0);

    // row 2
    for(int i=0; i<nnode; i++)
      Bl(1,i) = (*F)(0,1) * gradn(i,1);
    for(int i=nnode; i<2*nnode; i++)
      Bl(1,i) = (*F)(1,1) * gradn(i,1);
    for(int i=2*nnode; i<3*nnode; i++)
      Bl(1,i) = (*F)(2,1) * gradn(i,1);

    // row 3
    for(int i=0; i<nnode; i++)
      Bl(2,i) = (*F)(0,2) * gradn(i,2);
    for(int i=nnode; i<2*nnode; i++)
      Bl(2,i) = (*F)(1,2) * gradn(i,2);
    for(int i=2*nnode; i<3*nnode; i++)
      Bl(2,i) = (*F)(2,2) * gradn(i,2);

    // row 4
    for(int i=0; i<nnode; i++)
      Bl(3,i) = (*F)(0,0) * gradn(i,1) + (*F)(0,1) * gradn(i,0);
    for(int i=nnode; i<2*nnode; i++)
      Bl(3,i) = (*F)(1,0) * gradn(i,1) + (*F)(1,1) * gradn(i,0);
    for(int i=2*nnode; i<3*nnode; i++)
      Bl(3,i) = (*F)(2,0) * gradn(i,1) + (*F)(2,1) * gradn(i,0);

    // row 5
    for(int i=0; i<nnode; i++)
      Bl(4,i) = (*F)(0,1) * gradn(i,2) + (*F)(0,2) * gradn(i,1);
    for(int i=nnode; i<2*nnode; i++)
      Bl(4,i) = (*F)(1,1) * gradn(i,2) + (*F)(1,2) * gradn(i,1);
    for(int i=2*nnode; i<3*nnode; i++)
      Bl(4,i) = (*F)(2,1) * gradn(i,2) + (*F)(2,2) * gradn(i,1);

    // row 6
    for(int i=0; i<nnode; i++)
      Bl(5,i) = (*F)(0,0) * gradn(i,2) + (*F)(0,2) * gradn(i,0);
    for(int i=nnode; i<2*nnode; i++)
      Bl(5,i) = (*F)(1,0) * gradn(i,2) + (*F)(1,2) * gradn(i,0);
    for(int i=2*nnode; i<3*nnode; i++)
      Bl(5,i) = (*F)(2,0) * gradn(i,2) + (*F)(2,2) * gradn(i,0);
  }  
}

// void TotalLagrangian::Bl_matrix (const arma::mat & gradn, const arma::mat * F, 
//                                  arma::mat & Bl) const
// {
//   Bl.zeros();

//   int ndim  = msh.get_n_dim();
//   int nnode = gradn.n_rows / ndim;

//   if (ndim == 2)
//   {
//     // row 1 (E11)
//     for(int i=0; i<nnode; i++) Bl(0, i)       = (*F)(0,0) * gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(0, i+nnode) = (*F)(1,0) * gradn(i,0);

//     // row 2 (E22)
//     for(int i=0; i<nnode; i++) Bl(1, i)       = (*F)(0,1) * gradn(i,1);
//     for(int i=0; i<nnode; i++) Bl(1, i+nnode) = (*F)(1,1) * gradn(i,1);

//     // row 3 (2E12)
//     for(int i=0; i<nnode; i++) Bl(2, i)       = (*F)(0,0)*gradn(i,1) + (*F)(0,1)*gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(2, i+nnode) = (*F)(1,0)*gradn(i,1) + (*F)(1,1)*gradn(i,0);
//   }
//   else if (ndim == 3)
//   {
//     // row 1 (E11)
//     for(int i=0; i<nnode; i++) Bl(0, i)         = (*F)(0,0) * gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(0, i+nnode)   = (*F)(1,0) * gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(0, i+2*nnode) = (*F)(2,0) * gradn(i,0);

//     // row 2 (E22)
//     for(int i=0; i<nnode; i++) Bl(1, i)         = (*F)(0,1) * gradn(i,1);
//     for(int i=0; i<nnode; i++) Bl(1, i+nnode)   = (*F)(1,1) * gradn(i,1);
//     for(int i=0; i<nnode; i++) Bl(1, i+2*nnode) = (*F)(2,1) * gradn(i,1);

//     // row 3 (E33)
//     for(int i=0; i<nnode; i++) Bl(2, i)         = (*F)(0,2) * gradn(i,2);
//     for(int i=0; i<nnode; i++) Bl(2, i+nnode)   = (*F)(1,2) * gradn(i,2);
//     for(int i=0; i<nnode; i++) Bl(2, i+2*nnode) = (*F)(2,2) * gradn(i,2);

//     // row 4 (2E12)
//     for(int i=0; i<nnode; i++) Bl(3, i)         = (*F)(0,0)*gradn(i,1) + (*F)(0,1)*gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(3, i+nnode)   = (*F)(1,0)*gradn(i,1) + (*F)(1,1)*gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(3, i+2*nnode) = (*F)(2,0)*gradn(i,1) + (*F)(2,1)*gradn(i,0);

//     // row 5 (2E23)
//     for(int i=0; i<nnode; i++) Bl(4, i)         = (*F)(0,1)*gradn(i,2) + (*F)(0,2)*gradn(i,1);
//     for(int i=0; i<nnode; i++) Bl(4, i+nnode)   = (*F)(1,1)*gradn(i,2) + (*F)(1,2)*gradn(i,1);
//     for(int i=0; i<nnode; i++) Bl(4, i+2*nnode) = (*F)(2,1)*gradn(i,2) + (*F)(2,2)*gradn(i,1);

//     // row 6 (2E13)
//     for(int i=0; i<nnode; i++) Bl(5, i)         = (*F)(0,0)*gradn(i,2) + (*F)(0,2)*gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(5, i+nnode)   = (*F)(1,0)*gradn(i,2) + (*F)(1,2)*gradn(i,0);
//     for(int i=0; i<nnode; i++) Bl(5, i+2*nnode) = (*F)(2,0)*gradn(i,2) + (*F)(2,2)*gradn(i,0);
//   }
// }

void TotalLagrangian::Bnl_matrix (const arma::mat & gradn, arma::mat & B) const
{
  int nd = msh.get_n_dim();  // n_dim
  int nb = gradn.n_rows/nd; // n_node
  arma::mat I(arma::eye(nd,nd));
  B.zeros();
  B = arma::kron(I, gradn.submat(0,0,nb-1,nd-1).t() );  
}

void TotalLagrangian::pre_computation()
{  
  MixedFiniteElement * fe = fespace.createFE();

  int ndim  = fe->get_ndim();
  int ndof  = fe->get_ndofs_u();
  int nnode  = ndof/ndim;
  int nelem = msh.get_n_elements();

  cout << " Precomputing FE data" << endl;
  
  for(int i=0; i<nelem; i++)
  {
    std::vector<arma::vec3> elem_x0(nnode);
    get_element_x0(i, elem_x0);
    
    MxFEData * fd = new MxFEData(i, fe, elem_x0);
    fedata.push_back( fd );
  }

  delete fe;
}

void TotalLagrangian::rcg_tensor(const arma::mat & gradn, 
                                 const std::vector<arma::vec3> & elem_x, 
                                 double & J, arma::mat33 * pF, arma::mat & C)
{
  int ndim = gradn.n_cols;
  int nnod = gradn.n_rows;
  int nnode = nnod/ndim;
  arma::mat33 & F = *pF;

  F.zeros();
  if (ndim == 2) F(2,2) = 1.0;

  // computes deformation gradient tensor: F = dx/dX
  for(int id=0; id<ndim; id++)
    for(int jd=0; jd<ndim; jd++)
      for(int k=0; k<nnode; k++)
        F(id,jd) += elem_x[k][id] * gradn(k,jd);

  J = arma::det(F);

  // trick to work in 2D
  if (ndim==2 && J == 0.0) J = 1.0;

  // right Cauchy Green tensor
  C = F.t() * F;
}

void TotalLagrangian::pre_solve(){
  cout << "Solving nonlinear problem (Total Lagrangian)" << endl;
  cout << "Initial volume: " << calc_volume(true) << endl;
}

void TotalLagrangian::solve()
{
  static int cont = 1;
  static bool precomp = true;
  int nits;
  int num_prescribed = get_num_nz_prescribed();
  
  static NonlinearSolver * nls = NULL;

#ifdef USE_BFGS
  static PardisoSolver * parsolver = NULL;
#endif

  ofstream sfile;
  sfile.open("stress.dat");
  
  if(precomp)
  {
    // perform some FE pre computations for TL
    pre_computation();
    precomp = false;
    
    // create Newton nonlinear solver
    nls = new NewtonLineSearch(this);
    nls->init();    

#ifdef USE_BFGS    
    // create BFGS nonlinear solver
    parsolver = new PardisoSolver();    
    nls = new BFGSSolver(this);    
    ((BFGSSolver*)nls)->init(num_dofs, parsolver);
    
    cout << "BFGS setup done" << endl;
#endif    
  }
 
  if(output_step)
    output_vtk(cont,0);

  body_forces();
  assemble_traction();
  assemble_const();
  
  // load increment loop
  while( lc.has_load() )
  {
    lc.update();

#ifndef USE_BFGS    
    evaluate_forces( ((Newton*)nls)->residual() );
#endif

#ifdef USE_BFGS    
    evaluate_forces( ((BFGSSolver*)nls)->residual0() );    
#endif
    
    if (num_prescribed > 0 || neumann_map.size() > 0)
    {     
      prescribe_displacements();
      react = 0.0;
    }       

    // call the nonlinear solver
    timer.enter("Nonlinear solver");
    nits = nls->solve();
    timer.leave();
    
    // output deformation from this load increment
    if(output_step)
    {
      output_vtk(cont, lc.increment());
    }
/*
    if(output_step)
    {
      cout << "Writing stress to file" << endl;
      for(int kk=0; kk<9; kk++) // stress component
      {
        double mediastress = 0;
        for (int ii = 0; ii < 8; ii++) // integration point
          mediastress += stressdb(0, ii, kk);
        mediastress = mediastress/8.0;
        sfile.setf(ios::scientific);
        sfile << mediastress << "\t";
      }
      sfile << endl;
    }
*/
    first_step = true;
            
    // adapt load step size
    //lc.adapt(nits);
  }
 
  // prepare to leave
  cont  += 1;     // counter for output
  fext0 += fext ; // save nodal loads
  log << calc_volume() << endl;
  cout << " Done" << endl;
    
  sfile.close();
}



// --- BACKUP -----------------------------------------------------------------

/*
 void TotalLagrangian::solve()
{
  static int cont = 1;
  static bool precomp = true;
  int num_prescribed = get_num_nz_prescribed();

  // --------------------------------------------------------------------------
  // Setup linear solver for Newton step
  // --------------------------------------------------------------------------

  std::pair<int,double> ir;
  petsc::LinearSolver kspsolver;
  kspsolver.init();  
  kspsolver.use_mumps();
  kspsolver.set_ordering("rcm");

  // --------------------------------------------------------------------------
  // Setup and pre-computations
  // --------------------------------------------------------------------------
  
  cout << "Solving nonlinear problem (Total Lagrangian)" << endl;
  cout << "Initial volume: " << calc_mesh_volume() << endl;

  if(output_step) 
    output_vtk(cont,0);

  if(precomp)
  {
    pre_computation();
    precomp = false;
  }
  
  body_forces();
  assemble_traction();
  assemble_const();
  
  // --------------------------------------------------------------------------
  // Load increment loop
  // --------------------------------------------------------------------------

  while( lc.has_load() )
  {
    lc.update();
    evaluate_forces(r);
    assemble_pressure();
    
    if (num_prescribed > 0 || neumann_map.size() > 0)
    {
      prescribe_displacements();
      init_resid_stiff();
      assemble_pressure();
      assemble_stiff();
    }

    // ------------------------------------------------------------------------
    // Newton iteration
    // ------------------------------------------------------------------------
    int niter    = 0;
    int miter    = 50;
    double rtol  = 1.0e-6;
    double atol  = 1.0e-5;
    double rnorm = 2 * rtol;
    double rtoli = 1.0e-10;
    bool newtonconv = false;
    double enorm=0, enorm0=0;

    rnorm0 = react.l2norm() + arma::norm(lc.load()*fext,2);        
    
    while (!newtonconv && (niter < miter))
    {
      niter++;
        
      // Solve Ku = -r
      timer.enter("Solve");
      ir = kspsolver.solve(K, u, r, rtoli);
      timer.leave();

      if(niter == 1) enorm0 = calc_energy();

      // -----------------------------------------------------------------------
      // Line Search iteration
      //  Enable search -> 0.5 or disable -> search = 1e+6
      // -----------------------------------------------------------------------
      int nsearch = 0;
      int msearch = 5;
      double eta0 = 0.0, eta = 1.0;
      double rtu0 = calc_energy();
      double search = 0.5; 
      double rtu = rtu0 * search * 2;
      
      while( (abs(rtu) > abs(rtu0 * search)) && (nsearch < msearch) )
      {
        // increment number of line search steps
	      nsearch++;
               
	      // update geometry
	      update_geometry(eta-eta0);

	      // compute new residual forces and stiffness matrix
	      timer.enter("Assemble");	
        init_resid_stiff();    
	      assemble_pressure();
	      assemble_stiff();
	      timer.leave();

      	// perform line search
	      line_search(eta0, eta, rtu0, rtu);
      }

      // check for equilibrium convergence
      newtonconv = newton_converged(rnorm0,enorm0,rtol,atol,rnorm,enorm);

      // output information
      cout << fixed << scientific;
      cout << "  Newton step" << setprecision(4) << niter;
      cout << " rnorm " << rnorm << "  enorm " << enorm << endl;
    }

    if (niter >= miter)
    {
      kspsolver.converged_reason();
      ERROR("newton solver did not converge");
    }
    else
    {
      if(output_step) 
	      output_vtk(cont, lc.increment());

      log << niter << "\t" << ir.first << "\t" << endl;
    }   
  }
 
  // prepare to leave
  cont  += 1;     // counter for output
  fext0 += fext ; // save nodal loads
  log << calc_mesh_volume() << endl;
  cout << " Done" << endl;
}
*/
