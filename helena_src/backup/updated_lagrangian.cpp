#include "updated_lagrangian.hpp"

#ifdef USE_SNES
extern PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal norm, void* mctx);
#endif

/*
// MOVED
void UpdatedLagrangian::run(const string & name)
{
  set_mesh_filename(name);
  init();
  assemble_const();

#ifndef USE_SNES
  solve();
#else
  solve_snes();
#endif

}

// MOVED
void UpdatedLagrangian::init_mesh()
{
  // Read mesh
  std::cout << ":: Reading mesh" << endl;

  GmshIO gmshreader(msh);
  gmshreader.read(filename);

  //std::cout << "Done." << std::endl;
  std::cout << msh << std::endl;

  fespace.set_mesh(&msh);
}

void UpdatedLagrangian::init_coords()
{
  const std::vector<Point3d> & pts = msh.get_points();
  //x.reserve(pts.size());
  //x0.reserve(pts.size());
  for(uint i=0; i<pts.size(); i++)
  {
    const Point3d & p = pts[i];
    x.push_back(p);
    x0.push_back(p);
  }
}

void UpdatedLagrangian::init_control()
{
  // TODO: Change xlmax. User input.
  xlamb = 0.0;
  xlmax = 1.0;
  dlamb = 1.0/num_increments;
}

void UpdatedLagrangian::init_matvecs()
{
  int nod, dir, idx;
  int N = fespace.get_ndofs();
  NodalData inf;
  std::map<int,NodalData>::const_iterator it;

  num_dofs = N;

  udisp.resize(N);
  eload.resize(N);
  eload0.resize(N);

  udisp = 0.0;
  eload = 0.0;
  eload0 = 0.0;

  for(int i=0; i<N; i++)
    ldgof.push_back(true);

  //
  // Initialize ldgof array
  //  --> ldgof[i] = true, if dof i is active
  //  --> ldgof[i] = false, if dof i is prescribed
  //

  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    idx = nod + (dir * msh.get_n_points());
    ldgof[idx] = false;
  }

  // Initialize nodal loads vector
  for(it=nodal_loads_map.begin(); it!=nodal_loads_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    idx = nod + (dir * msh.get_n_points());
    eload[idx] = inf.second;
  }

  // Allocate matrices and vectors
  if (first)
  {
    K.create(num_dofs,num_dofs,120);
    u.create(num_dofs);
    //f.create(ndofs);
    r.create(num_dofs);
    t.create(num_dofs);
    react.create(num_dofs);
    ubar.create(num_dofs);
  }

  // Initialize residual and reactions
  K = 0.0;
  r = 0.0;
  react = 0.0;
}


void UpdatedLagrangian::init_resid_stiff()
{
  K = 0.0;
  u = 0.0;
  react = 0.0;
  //r = 0.0;

  for(int i=0; i<fespace.get_ndofs(); i++)
    r.set(i, eload0(i) + xlamb * eload(i) );
    //r.set(i, xlamb * eload[i]);

#ifdef USE_SNES
  //form_jacobian();
  //apply_boundary();
#endif
}


// MOVED

void UpdatedLagrangian::save_mesh()
{
  std::stringstream ss;
  ss << this->output;
  std::string output_file = basename + "_deformed_" + ss.str();
  this->output_vtk(output_file);
  output++;
}

void UpdatedLagrangian::init()
{
  if(!restart)
  {
    init_mesh();
    init_coords();
    init_control();
    init_matvecs();
    save_mesh();

    // Write undeformed mesh
    //std::string output_file = basename + "_0";
    //output_vtk(output_file);

    //std::stringstream ssi;
    //ssi << cont;
    //ssi << "_";
    //ssi << increment;
    //std::string output_file = basename + "_" + ssi.str();
    //this->output_vtk(output_file);

#ifdef USE_SNES
    int N = fespace.get_ndofs();
    int npts = msh.get_n_points();
    const std::vector<Point3d> & pts = msh.get_points();
    X.create(N);
    X0.create(N);
    xaux.resize(N);

    // Initial coordinates
    for(int i=0; i<npts; i++)
      for(int j=0; j<msh.get_ndim(); j++)
      {
        int k = i + j * npts;
        x[i][j]  = pts[i][j];
        X.set(k,x[i][j]);
        X0.set(k,x[i][j]);
      }
#endif

  }
  else
  {
    // TODO: Implement restar1, restar2, dump functions!
  }

}

void UpdatedLagrangian::reset()
{
  // PETSc
  K = 0.0;
  u = 0.0;
  udisp = 0.0;

  /////////////////r = 0.0;
  //react = 0.0;
  /////////////////eload = 0.0;

  init_control();
}
*/

void UpdatedLagrangian::initial_stress(const Vector & is, std::vector<Matrix*> & vstrs)
{
  std::cout << " Loading initial stress\n";

  int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();
  Matrix elmat(n,n);
  Matrix elmat_const(n,n);
  Matrix elmat_sigma(n,n);
  Vector tmp_eload(fespace.get_ndofs());
  Vector elvec(n);
  std::vector<int> dnums;
  MixedFiniteElement * fe = fespace.createFE();

  tmp_eload = 0.0;

  // Loop over all the elements
  for(int i=0; i<nelem; i++)
  {
    //calc_elmatvec (i, fe, elmat_const, elmat_sigma, elvec, dnums);

    // BEGIN CALC ELMAT VEC
    int numc = num_components;
    int ndim = fe->get_ndim();
    int ndof = fe->get_ndofs_u();
    int nubf = ndof/ndim;
    int vsize = (ndim == 2 ? 3 : 6);

    Real detJxW, detF;
    Vector shape;
    Matrix dshape;
    Matrix gradn(ndof,ndim);
    Matrix jacinv(ndim,ndim);
    Matrix btens(3,3), sigma(3, 3), sigma_voigt(vsize,1), sigma2(2,2);
    Matrix B(numc,ndof), D(numc,numc), H(ndim,ndof);
    Point3d pt, qpt;

    std::vector<Point3d> xe(ndof);
    std::vector<Point3d> x0(ndof);

    this->get_element_x(i, xe);
    this->get_element_x0(i, x0);

    Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type()); // CHECK
    Mapping em = fe->get_mapping(i, xe);

    elmat_const = 0.0;
    elmat_sigma = 0.0;
    elvec = 0.0;

    fespace.get_element_dofs_u(i,dnums);

    for(int q=0; q<qd->get_num_ipoints(); q++)
    {
      //d.get_point(q, qpt);
      qpt = qd->get_point(q);

      fe->calc_shape_u(qpt,shape);
      fe->calc_deriv_shape_u(qpt,dshape);

      em.calc_jacobian(dshape, nubf);
      detJxW = qd->get_weight(q) * em.get_det_jacobian();
      jacinv = em.get_inv_jacobian();

      gradn = (dshape * jacinv);

      // Compute left Cauchy-Green tensor b=FF^{T}
      calc_left_cauchy_green(gradn, x0, &detF, btens);

      calc_B_matrix (gradn, B);
      //calc_H_matrix (shape, H);

      material->cauchy_stress(ndim, detF, btens, sigma);
      material->elas_tensor(elastype, detF, D);

      // Mount sigma in Voigt notation
      //material->cauchy_stress_voigt(ndim, detF, btens, sigma_voigt);
      if(ndim==2)
      {       
	sigma_voigt(0,0) = sigma(0,0);
	sigma_voigt(1,0) = sigma(1,1);
	sigma_voigt(2,0) = sigma(0,1);
      }
      else if (ndim==3)
      {
	sigma_voigt(0,0) = sigma(0,0);
	sigma_voigt(1,0) = sigma(1,1);
	sigma_voigt(2,0) = sigma(2,2);
	sigma_voigt(3,0) = sigma(0,1);
	sigma_voigt(4,0) = sigma(1,2);
	sigma_voigt(5,0) = sigma(0,2);
      }

      // *********************************************************************
      // INTERPOLATE Ta to integration POINTS
      // *********************************************************************
      Matrix *lstress = vstrs[i];
      *lstress=0.0;

      // *********************************************************************
      // HARDCODED FOR QUADS *************************************************
      // *********************************************************************
      for(int j=0; j<4; j++)
      {
        int jj = dnums[j]; // no global

        (*lstress)(0,0) += 0.0; //shape(j) * is(jj);
        (*lstress)(0,1) += 0.0;//shape(j) * is(jj);
        //(*lstress)(1,1) += shape(j) * is(jj);
        (*lstress)(1,1) += shape(j) * is(jj);
        //
        // TODO: colocar os vetores das fibras
        //+ eledb(1,jn,in) * (xvecf(jp,1)*xta(jp)*xvecf(jp,1))
        //
      }

      //////index++;
      //
      // *********************************************************************

      // Indicial notation
      for(int i=0; i<nubf; i++)
        for(int id=0; id<ndim; id++)
          for(int jd=0; jd<ndim; jd++)
            elvec(i+id*nubf) += (*lstress)(id,jd) * gradn(i,jd) * detJxW;

    }

    delete qd;
    // END CALC ELMAT VEC

    // Residual vector and reaction vector
    for(int k=0;k<n;k++){

      if (ldgof[dnums[k]])
        tmp_eload(dnums[k]) -= elvec(k);
      //r.add(dnums[k], -elvec(k)); // -R
      //else
      //  react.add(dnums[k], elvec(k));
    }
  }

  for(int j=0; j<eload.size(); j++)
    eload(j) = tmp_eload(j);

  delete fe;
}

void UpdatedLagrangian::output_F_tensor(std::vector<Matrix*> & v_ftens)
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
    Vector shape;
    Matrix dshape;
    Matrix gradn(ndof,ndim);
    Matrix jacinv(ndim,ndim);
    Matrix * F = new Matrix(3,3);

    Point3d pt, qpt;

    std::vector<Point3d> xe(ndof);
    std::vector<Point3d> x0(ndof);
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

void UpdatedLagrangian::output_C_tensor(std::vector<Matrix*> & v_ctens)
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
    Vector shape;
    Matrix dshape;
    Matrix gradn(ndof,ndim);
    Matrix jacinv(ndim,ndim);
    Matrix * Ctens = new Matrix(3,3);

    Point3d pt, qpt;

    std::vector<Point3d> xe(ndof);
    std::vector<Point3d> x0(ndof);
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

/*
void UpdatedLagrangian::apply_boundary()
{
  // Impose prescribed displacements
  if (fixed_nodes_map.size())
  {
    int nod, idx, dir;
    Real val;
    NodalData inf;

    std::map<int,Real> boundary_values;
    std::map<int,NodalData>::iterator it;

    for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
    {
      nod = it->first;
      inf = it->second;
      dir = inf.first;
      val = dlamb * inf.second;
      idx = (dir * msh.get_n_points()) + nod;

      //cout << "idx = " << idx << endl;

      boundary_values.insert( std::pair<int,Real>(idx,val) );
    }
    FETools::apply_boundary_values(boundary_values, K, r, ubar, false);
  }
}
*/

void UpdatedLagrangian::assemble_const()
{
  int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();
  std::vector<int> dnums;
  Matrix elmat(n,n);
  MixedFiniteElement * fe = fespace.createFE();

  for(int i=0; i<nelem; i++)
  {
    calc_elmat_const (i, fe, elmat);
    fespace.get_element_dofs_u(i,dnums);

    for(int k=0;k<n;k++)
      for(int l=0;l<n;l++)
        K.add(dnums[k],dnums[l],elmat(k,l));
  }

  K.assemble();
  apply_boundary();

  delete fe;
}

void UpdatedLagrangian::assemble_stiff()
{
  int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();
  std::vector<int> dnums;
  Matrix elmat(n,n);
  Matrix elmat_const(n,n);
  Matrix elmat_sigma(n,n);
  Vector elvec(n);

  MixedFiniteElement * fe = fespace.createFE();

  // Loop over all the elements
  for(int i=0; i<nelem; i++)
  {
    // calc_elmatvec substitutes subroutine: internal, kconst, ksigma
    calc_elmatvec (i, fe, elmat_const, elmat_sigma, elvec, dnums);

    fespace.get_element_dofs_u(i,dnums);

    //cout << "elmat-const  \n" << elmat_const << endl;
    //cout << "elmat-sigma  \n" << elmat_sigma << endl;
    //cout << "elvec        \n" << elvec << endl;

    // K_const + K_sigma
    for(int k=0;k<n;k++)
      for(int l=0;l<n;l++)
        K.add(dnums[k],dnums[l], elmat_const(k,l)+elmat_sigma(k,l) );

    // Residual vector and reaction vector
    for(int k=0;k<n;k++)
    {
      if (ldgof[dnums[k]])
        r.add(dnums[k], -elvec(k)); // -R
      else
        react.add(dnums[k], elvec(k));
    }
  }

  // Assemble matrix and vector and
  // impose prescribed displacements
  K.assemble();
  r.assemble();
  apply_boundary();

  delete fe;
}

void UpdatedLagrangian::calc_left_cauchy_green(const Matrix & gradn,
    const std::vector<Point3d> & x0, Real * detf, Matrix & btens)
{
  int ndim = gradn.width();
  int nnod = gradn.height();
  int nubf = nnod/ndim;
  Matrix finvr(3,3), ftens(3,3);
  finvr = 0.0;

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

  // Obtains b = F F^T
  btens = 0.0;
  if (ndim == 2) btens(2, 2) = 1.0;

  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      btens(i,j) = 0.0;
      for(int k=0; k<ndim; k++){
        btens(i,j) += ftens (i,k) * ftens(j,k);
      }
    }
  }

#ifdef DEBUG
  cout << "\nCartesian derivatives of shape functions\n";
  cout << gradn << endl;
  cout << "F^{-1}\n";
  cout << finvr << endl;
  cout << "Det F\n";
  cout << *detf << endl;
  cout << "\nDeformation Gradient tensor F\n";
  cout << ftens << endl;
  cout << "Left Cauchy-Green tensor\n";
  cout << btens << endl;
#endif

}

void UpdatedLagrangian::calc_deformation_grad(const Matrix & gradn,
    const std::vector<Point3d> & x0, Real * detf, Matrix & ftens)
{
  int ndim = gradn.width();
  int nnod = gradn.height();
  int nubf = nnod/ndim;
  Matrix finvr(3,3);//, ftens(3,3);
  finvr = 0.0;

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

void UpdatedLagrangian::calc_right_cauchy_green(const Matrix & gradn,
    const std::vector<Point3d> & x0, Real * detf, Matrix & ctens)
{
  int ndim = gradn.width();
  int nnod = gradn.height();
  Matrix finvr(3,3), ftens(3,3);
  finvr = 0.0;

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
  ctens = 0.0;
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

void UpdatedLagrangian::calc_elmat_const (const int eindex,
    const MixedFiniteElement * fe, Matrix & elmat)
{
  int numc = num_components;
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;

  Real detJxW;
  Matrix dshape;
  Matrix gradn(ndof,ndim);
  Matrix jacinv(ndim,ndim);
  Matrix B(numc,ndof), D(numc,numc);
  std::vector<Point3d> xe(nubf);
  Point3d qpt;

  this->get_element_x(eindex, xe);

  Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type()); // CHECK
  Mapping em = fe->get_mapping(eindex,xe);

  elmat = 0.0;

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    // *******************************************************
    // IMPORTANT !!! Correct the BUG in all quadrature rules *
    // *******************************************************

    //qd.get_point(q, qpt);
    qpt = qd->get_point(q);

    fe->calc_deriv_shape_u(qpt,dshape);
    em.calc_jacobian(dshape, nubf);

    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn  = (dshape * jacinv);

    if(elastype == PLANE_STRESS)
    {
      PlaneStressNeoHookean * pm = (PlaneStressNeoHookean*) material;
      detJxW = detJxW * pm->get_thick();
    }

    calc_B_matrix (gradn, B);
    material->elas_tensor(elastype, 1.0, D);

    elmat += (transpose(B) * D * B) * detJxW;

  }

  delete qd;
}

void UpdatedLagrangian::calc_elmatvec (const int eindex,
    const MixedFiniteElement * fe, Matrix & elmat_const,
    Matrix & elmat_sigma, Vector & elvec, std::vector<int> lnods)
{
  //int index;
  int numc = num_components;
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;
  int vsize = (ndim == 2 ? 3 : 6);
  int ssize = (ndim == 2 ? 4 : 9); // stress matrix - voigt size

  Real detJxW, detF;
  Vector shape;
  Matrix dshape;
  Matrix gradn(ndof,ndim);
  Matrix jacinv(ndim,ndim);
  Matrix btens(3,3), sigma(3, 3), sigma_voigt(vsize,1), sigma2(ssize,ssize);
  Matrix B(numc,ndof), D(numc,numc);
  Matrix H;

  if (ndim == 2)
    H.resize(2*ndim,ndof);
  else if (ndim == 3)
    H.resize(3*ndim,ndof);

  Point3d pt, qpt;

  std::vector<Point3d> elem_xe(nubf);
  std::vector<Point3d> elem_x0(nubf);

  this->get_element_x(eindex, elem_xe);
  this->get_element_x0(eindex, elem_x0);

  Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type()); // CHECK
  Mapping em = fe->get_mapping(eindex, elem_xe);

  elmat_const = 0.0;
  elmat_sigma = 0.0;
  elvec = 0.0;

  for(int q=0; q<qd->get_num_ipoints(); q++)
  {
    //qd.get_point(q, qpt);
    qpt = qd->get_point(q);

    fe->calc_shape_u(qpt,shape);
    fe->calc_deriv_shape_u(qpt,dshape);

    em.calc_jacobian(dshape, nubf);
    detJxW = qd->get_weight(q) * em.get_det_jacobian();
    jacinv = em.get_inv_jacobian();
    gradn  = (dshape * jacinv);

    calc_left_cauchy_green(gradn, elem_x0, &detF, btens);
    calc_B_matrix (gradn, B);
    calc_H_matrix (gradn, H);

    //
    // ******* BUG *****
    // TOMAR CUIDADO COM AS IMPLEMENTACOES DO SIGMA_VOIGT!
    // NAO ESQUECER!
    //
    material->cauchy_stress(ndim, detF, btens, sigma);
    material->elas_tensor(elastype, detF, D);

    sigma2 = 0.0;

    if (ndim==2)
    {
      sigma2(0,0)=sigma(0,0);
      sigma2(0,1)=sigma(0,1);
      sigma2(1,0)=sigma(1,0);
      sigma2(1,1)=sigma(1,1);
      sigma2(2,2)=sigma(0,0);
      sigma2(2,3)=sigma(0,1);
      sigma2(3,2)=sigma(1,0);
      sigma2(3,3)=sigma(1,1);

      sigma_voigt(0,0) = sigma(0,0);
      sigma_voigt(1,0) = sigma(1,1);
      sigma_voigt(2,0) = sigma(0,1);
    }
    else if (ndim == 3)
    {
      // block 1
      sigma2(0,0) = sigma(0,0);
      sigma2(0,1) = sigma(0,1);
      sigma2(0,2) = sigma(0,2);
      sigma2(1,0) = sigma(1,0);
      sigma2(1,1) = sigma(1,1);
      sigma2(1,2) = sigma(1,2);
      sigma2(2,0) = sigma(2,0);
      sigma2(2,1) = sigma(2,1);
      sigma2(2,2) = sigma(2,2);
      // block 2
      sigma2(3,3) = sigma(0,0);
      sigma2(3,4) = sigma(0,1);
      sigma2(3,5) = sigma(0,2);
      sigma2(4,3) = sigma(1,0);
      sigma2(4,4) = sigma(1,1);
      sigma2(4,5) = sigma(1,2);
      sigma2(5,3) = sigma(2,0);
      sigma2(5,4) = sigma(2,1);
      sigma2(5,5) = sigma(2,2);
      // block 3
      sigma2(6,6) = sigma(0,0);
      sigma2(6,7) = sigma(0,1);
      sigma2(6,8) = sigma(0,2);
      sigma2(7,6) = sigma(1,0);
      sigma2(7,7) = sigma(1,1);
      sigma2(7,8) = sigma(1,2);
      sigma2(8,6) = sigma(2,0);
      sigma2(8,7) = sigma(2,1);
      sigma2(8,8) = sigma(2,2);

      sigma_voigt(0,0) = sigma(0,0);
      sigma_voigt(1,0) = sigma(1,1);
      sigma_voigt(2,0) = sigma(2,2);
      sigma_voigt(3,0) = sigma(0,1);
      sigma_voigt(4,0) = sigma(1,2);
      sigma_voigt(5,0) = sigma(0,2);
    }

    if(elastype == PLANE_STRESS)
    {
      PlaneStressNeoHookean * pm = (PlaneStressNeoHookean*) material;
      detJxW = detJxW * (pm->get_thick()/detF);
    }

    // Computations in matricial form

    elvec       +=  transpose(B) * sigma_voigt * detJxW;
    elmat_const += (transpose(B) * D * B) * detJxW;
    elmat_sigma += (transpose(H) * sigma2 * H) * detJxW;

    // Computations in indicial form

    // Internal equivalent nodal forces
    //for(int i=0; i<nubf; i++)
    //  for(int id=0; id<ndim; id++)
    //    for(int jd=0; jd<ndim; jd++)
    //      elvec(i+id*nubf) += sigma(id,jd) * gradn(i,jd) * detJxW;

     /*
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
	    //elmat_sigma(i+kd*nubf,j+kd*nubf) += gradn(i,kd) * sigma(kd,ld) * gradn(j,ld) * detJxW;

        for(int id=0; id<ndim; id++)
          elmat_sigma(i+id*nubf,j+id*nubf) += sum;
      }
    */

  }

  delete qd;

}

void UpdatedLagrangian::calc_elstress(Mapping & em, const MixedFiniteElement * fe,
    const Vector & u, const Point3d & X, Vector & sigma)
{
  int index;
  int numc = num_components;
  int ndim = fe->get_ndim();
  int ndof = fe->get_ndofs_u();
  int nubf = ndof/ndim;

  sigma.resize(numc);

  Matrix dshape;
  Matrix gradT_jac(ndof,ndim);
  Matrix jacinv(ndim,ndim);
  Matrix B(numc,ndof), D(numc,numc);

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

inline void UpdatedLagrangian::calc_B_matrix (const Matrix & gradn, Matrix & B) const
{
  B=0.0;

  int Nr = gradn.height(); // ndof
  int Nd = msh.get_ndim(); // ndim
  int Nb = Nr/Nd;          // nubf

  if (Nd == 2)
  {
    // e,x
    for(int i=0; i<Nb; i++)
      B(0,i) = gradn(i,0);

    // e,y
    for(int i=Nb; i<2*Nb; i++)
      B(1,i) = gradn(i,1);

    // ex,y
    for(int i=0; i<Nb; i++)
      B(2,i) = gradn(i,1);
    for(int i=Nb; i<2*Nb; i++)
      B(2,i) = gradn(i,0);
  }
  else if (Nd == 3)
  {
    // e,x
    for(int i=0; i<Nb; i++)
      B(0,i) = gradn(i,0);
    // e,y
    for(int i=Nb; i<2*Nb; i++)
      B(1,i) = gradn(i,1);
    // e,z
    for(int i=2*Nb; i<3*Nb; i++)
      B(2,i) = gradn(i,2);

    // ex,y
    for(int i=0; i<Nb; i++)
      B(3,i) = gradn(i,1);
    for(int i=Nb; i<2*Nb; i++)
      B(3,i) = gradn(i,0);
    // ey,z
    for(int i=Nb; i<2*Nb; i++)
      B(4,i) = gradn(i,2);
    for(int i=2*Nb; i<3*Nb; i++)
      B(4,i) = gradn(i,1);
    // ex,z
    for(int i=0; i<Nb; i++)
      B(5,i) = gradn(i,2);
    for(int i=2*Nb; i<3*Nb; i++)
      B(5,i) = gradn(i,0);
  }
}

inline void UpdatedLagrangian::calc_D_matrix (const Point3d & X, int index, Matrix & D,
    Real xmu, Real xlamb) const
{
  D=0.0;
  if (elastype == PLANE_STRAIN || elastype == PLANE_STRESS)
  {
    const Real c = xlamb + 2.0*xmu;
    D(0,0) = c;
    D(0,1) = xlamb;
    D(1,0) = xlamb;
    D(1,1) = c;
    D(2,2) = xmu;
  }
  else if (elastype == THREE_DIM)
  {
  }
}

inline void UpdatedLagrangian::calc_H_matrix (const Matrix & gradn, Matrix & H) const
{
  H=0.0;

  int Nr = gradn.height();    // ndof
  int Nd = msh.get_ndim();    // ndim
  int Nb = Nr/Nd;             // nubf

  if (Nd == 1) // line
  {
    // to implement
  }
  else if ((Nd == 2 && Nb == 3) || (Nd == 2 && Nb == 4)) // tri or quad
  {
    for(int i=0; i<Nb; i++)
    {
      H(0,i) = gradn(i%Nb,0);
      H(1,i) = gradn(i%Nb,1);
    }

    for(int i=Nb; i<2*Nb; i++)
    {
      H(2,i) = gradn(i%Nb,0);
      H(3,i) = gradn(i%Nb,1);
    }
  }
  else if (Nd == 3)
  {
    for(int i=0; i<Nb; i++)
    {
      H(0,i) = gradn(i%Nb,0);
      H(1,i) = gradn(i%Nb,1);
      H(2,i) = gradn(i%Nb,2);
    }
    for(int i=Nb; i<2*Nb; i++)
    {
      H(3,i) = gradn(i%Nb,0);
      H(4,i) = gradn(i%Nb,1);
      H(5,i) = gradn(i%Nb,2);
    }
    for(int i=2*Nb; i<3*Nb; i++)
    {
      H(6,i) = gradn(i%Nb,0);
      H(7,i) = gradn(i%Nb,1);
      H(8,i) = gradn(i%Nb,2);
    }
  }

}

void UpdatedLagrangian::solve()
{
  // Starts the load increment loop
  // TODO: Move to a member function later ---> load_increment_loop()

  int increment = 0;
  int num_nz_prescribed = get_num_nz_prescribed();
  static int cont = 1;
  std::string textmsg;
  std::pair<PetscInt,PetscReal> ir;
  PETSc::LinearSolver kspsolver;

  if(first)
  {
    kspsolver.init();
  }

  kspsolver.set_solver_type(1);    // CG
  kspsolver.set_preconditioner(3); // JACOBI

  cout << "-----------------------------------------";
  cout << "-----------------------------------------";
  cout << endl;
  cout << "Solving nonlinear problem\n";

  // Save undeformed mesh at each step
  std::stringstream ss;
  ss << cont;
  ss << "_";
  ss << 0;
  std::string output_fileu = basename + "_" + ss.str();
  this->output_vtk(output_fileu);

  // **************************************************************************
  // Load increment loop
  // **************************************************************************

  while ( (xlamb < xlmax) && (increment < num_increments))
  {
    xlamb += dlamb;
    increment++;

    evaluate_forces();
    boundary_normal_pressure();

    cout << "-----------------------------------------";
    cout << "-----------------------------------------";
    cout << endl;
    cout << "Load increment " << increment << endl;

    if (num_nz_prescribed > 0)
    {
      prescribe_displacements();
      init_resid_stiff();
      boundary_normal_pressure();
      assemble_stiff();
    }

    // ************************************************************************
    // Newton-Raphson iteration
    // ************************************************************************

    int niter = 0;
    int miter = 50;
    Real cnorm = 1.0e-6;
    Real rnorm = 2 * cnorm;

    r.assemble();

    while ((rnorm > cnorm) && (niter < miter))
    {
      niter++;

      // Solve Ku = -R
      ir = kspsolver.solve(K,u,r,1.e-16);

      // Update geometry
      //udisp = 0.0;
      u.get_data(udisp.ref());
      update_geometry(udisp);

      // Compute new residual forces and stiffness matrix
      init_resid_stiff();
      boundary_normal_pressure();
      assemble_stiff();

      // Check for equilibrium convergence
      newton_converged(rnorm);

      // Output information
      cout << " Newton iteration " << setw(3) << niter;
      cout << " ( ksp its " << setw(3) << ir.first;
      cout << " rnorm " << setw(3) << fixed << scientific << ir.second;
      cout << " ) , ( resid " << rnorm << " )" << endl;
    }

    if (niter >= miter)
    {
      kspsolver.converged_reason();
      error("updated_lagrangian.cpp","solve","newton solver did not converge");
    }
    else
    {
      // Save deformed mesh at each step
      std::stringstream ssi;
      ssi << cont;
      ssi << "_";
      ssi << increment;
      std::string output_file = basename + "_" + ssi.str();
      this->output_vtk(output_file);
    }
  }

  cout << "-----------------------------------------";
  cout << "-----------------------------------------" << endl;
  cout << "End of load increment loop" << endl;
  cout << "-----------------------------------------";
  cout << "-----------------------------------------" << endl;

  dump("Deformed configuration and nodal forces\n");

  cont++;
  first = false;

  // Save nodal loads
  for(int i=0; i<fespace.get_ndofs(); i++)
    eload0(i) += eload(i);
}

/*
void UpdatedLagrangian::check_convergence(Real * rnorm)
{
  Real enorm = 0.0, fnorm = 0.0;
  Vector raux(fespace.get_ndofs());
  Vector eaux(fespace.get_ndofs());

  r.get_data(raux.ref());
  react.get_data(eaux.ref());

  *rnorm = 0.0;

  for(int i=0; i<raux.size(); i++)
  {
    if(ldgof[i])
    {
      *rnorm += raux[i] * raux[i];
      fnorm += eload(i) * eload(i);
    }

    if (!ldgof[i])
      enorm += eaux[i] * eaux[i];
  }

  fnorm = fnorm * xlamb * xlamb;

  *rnorm = sqrt(*rnorm/(enorm+fnorm));

  if (isnan(*rnorm))
  {
    std::cerr << "\n Error: residual norm is NaN.\n";
    std::cerr << "   rnorm=" << *rnorm << std::endl;
    std::cerr << "   enorm=" << enorm << std::endl;
    std::cerr << "   fnorm=" << fnorm << std::endl;

    error("nonlinear_elasticity.cpp",
          "checking convergence of Newton-Raphson method",
          "NaN");
  }
}

// MOVED
void UpdatedLagrangian::dump(const std::string & msg)
{
  int aux = msh.get_n_points();
  Vector vreact(fespace.get_ndofs());
  react.get_data(vreact.ref());

  log->write_line();
  log->write(msg.c_str());
  log->write_line();

  for(int i=0; i<msh.get_n_points(); i++)
  {
    log->write("%f \t %f", x[i][0], x[i][1]);

    if (ldgof[i])
      log->write("\t %f ", xlamb*eload(i));
    else
      log->write("\t %f ", vreact(i));

    if (ldgof[i+aux])
      log->write("\t %f ", xlamb*eload(i+aux));
    else
      log->write("\t %f ", vreact(i+aux));

    log->write("\n");
  }
}

// MOVED
void UpdatedLagrangian::evaluate_forces()
{
  for(int i=0; i<fespace.get_ndofs(); i++)
    r.add(i, dlamb * eload[i]);
}

void UpdatedLagrangian::get_elvec(const std::vector<int> & vdnums, const Vector & u, Vector & elu)
{
  int nvd = vdnums.size();
  elu.resize(nvd);
  for(int i=0;i<nvd;i++)
    elu(i) = u(vdnums[i]);
}

// MOVED
void UpdatedLagrangian::get_x(int e, std::vector<Point3d> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_element_pt_nums(e,ptnums);

#ifdef USE_SNES
  uint j;
  uint nd = msh.get_ndim();
  uint np = msh.get_n_points();

  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    for(j=0; j<nd; j++)
    {
      xe[i](j) = xaux[k+j*np];
    }
  }
#else
  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x[k];
  }
#endif

}

void UpdatedLagrangian::get_x0(int e, std::vector<Point3d> & xe)
{
  uint i, k;
  std::vector<int> ptnums;
  msh.get_element_pt_nums(e,ptnums);

  for(i=0; i<ptnums.size(); i++)
  {
    k = ptnums[i];
    xe[i] = x0[k];
  }
}

// MOVED
int UpdatedLagrangian::get_num_nz_prescribed()
{
  int cont=0;
  Real val;
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

void UpdatedLagrangian::output_vtk(const string & name)
{
  Vector vmstress(msh.get_n_elements());
  string vtuname = name + ".vtu";
  von_mises(udisp,vmstress);
  vtkout.write_def_mesh(x, vmstress, vtuname);
}

void UpdatedLagrangian::output_vtk_only(const string & name)
{
  filename = name;
  string meshname = name + ".vtu";
  init_mesh();
  vtkout.write_mesh(meshname);
}

void UpdatedLagrangian::output_vtk_solution(const string & name, const Vector & u)
{
  string vtuname = name + ".vtu";
  vtkout.write_nodal_data_def_mesh(x, u, vtuname);
}

// MOVED
void UpdatedLagrangian::prescribe_displacements()
{
  int nod, dir;
  Real val;
  NodalData inf;
  std::map<int,Real> boundary_values;
  std::map<int,NodalData>::iterator it;

  // Loop over nodes with prescribed displacement
  cout << " Prescribing displacements\n";

#ifndef USE_SNES
  ubar = 0.0;
  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    val = inf.second;

    x[nod][dir] = x0[nod][dir] + xlamb * val;

    //idx = (dir * msh.get_n_points()) + nod;
    //ubar.set(idx, xlamb * val);
  }
#else
  int idx;
  ubar = 0.0;
  for(it=fixed_nodes_map.begin(); it!=fixed_nodes_map.end(); ++it)
  {
    nod = it->first;
    inf = it->second;
    dir = inf.first;
    val = inf.second;
    idx = (dir * msh.get_n_points()) + nod;

    x[nod][dir] = x0[nod][dir] + xlamb * val;
    X.add(idx, dlamb * val);
    ubar.set(idx, xlamb * val);
  }
#endif
}
*/


#ifdef USE_SNES
void UpdatedLagrangian::solve_snes()
{

  Mat J;
  Vec res, x;
  SNES snes;
  KSP ksp;
  PC pc;
  PetscErrorCode ierr;
  PetscInt its;
  SNESConvergedReason reason;

  int num_nz_prescribed = get_num_nz_prescribed();

  // Function and solution vectors
  ierr = VecCreate(PETSC_COMM_WORLD,&x);
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,fespace.get_ndofs());
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = VecSetFromOptions(x);
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = VecDuplicate(x,&res);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  // Matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&J);
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,fespace.get_ndofs(),fespace.get_ndofs());
  CHKERRABORT(MPI_COMM_WORLD, ierr);
  ierr = MatSetFromOptions(J);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  // SNES setup
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = SNESSetFunction(snes, res, FormFunction, this);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = SNESSetJacobian(snes,J,J,FormJacobian, this);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = SNESGetKSP(snes,&ksp);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = KSPGetPC(ksp,&pc);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = PCSetType(pc,PCNONE);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = KSPSetTolerances(ksp,1.e-8,PETSC_DEFAULT,PETSC_DEFAULT,500);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  ierr = SNESSetFromOptions(snes);
  CHKERRABORT(MPI_COMM_WORLD, ierr);

  cout << "------------------------------------------------------------------------------" << endl;
  cout << "Solving nonlinear problem\n";

  while ( (xlamb < xlmax) && (increment < num_increments))
  {
    xlamb += dlamb;
    increment++;

    cout << "------------------------------------------------------------------------------" << endl;
    cout << "Load increment " << increment << endl;
    cout << "------------------------------------------------------------------------------" << endl;

    //if (fixed_nodes_map.size() > 0)
    //  prescribe_displacements();

    if (num_nz_prescribed > 0)
    {
      prescribe_displacements();
    }


    // ************************************************************************
    // SNES SOLVE
    // ************************************************************************

    // Initial guess
    ierr = VecCopy(X.vec(), x);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

    ierr = SNESSolve(snes, PETSC_NULL, x);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

    ierr = SNESGetIterationNumber(snes, &its);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

    ierr = SNESGetConvergedReason(snes, &reason);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of Newton iterations = %D, %s\n",
        its,SNESConvergedReasons[reason]);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

    // Save deformed mesh at each step
    X.get_data(xaux.ref());
    uint nd = msh.get_ndim();
    uint np = msh.get_n_points();

    for(int i=0; i<np; i++)
    {
      for(int j=0; j<nd; j++)
      {
        this->x[i](j) = xaux[i+j*np];
      }
    }

    std::stringstream ss;
    ss << increment;
    std::string output_file = basename + "_" + ss.str();
    this->output_vtk(output_file);

    //for(int i=0; i<msh.get_n_points(); i++)
    //  cout << scientific << xaux[i] << "\t" << xaux[i+msh.get_n_points()] << "\n";

  }

  X.get_data(xaux.ref());

  cout << "------------------------------------------------------------------------------" << endl;
  cout << "End of load increment loop" << endl;
  cout << "------------------------------------------------------------------------------" << endl;

  log->write_line();
  log->write("Deformed configuration");
  log->write_line();
  for(int i=0; i<msh.get_n_points(); i++)
    log->write("%f \t %f", xaux[i], xaux[i+msh.get_n_points()]);

  //cout << "Deformed configuration\n";
  //cout << "------------------------------------------------------------------------------" << endl;
  //for(int i=0; i<msh.get_n_points(); i++)
  //  cout << scientific << xaux[i] << "\t" << xaux[i+msh.get_n_points()] << "\n";

}
#endif

/*
// Calculate Von Mises Stress of Mixed FE Variable
void UpdatedLagrangian::von_mises(const Vector & u, Vector & sigma)
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
    Vector elvecu, elsigma;

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
*/

#ifdef USE_SNES

#ifdef USE_SNES
void UpdatedLagrangian::form_function()
{
  int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();

  Matrix elmat_const(n,n);
  Matrix elmat_sigma(n,n);
  Vector elvec(n);

  std::vector<int> dnums;

  r = 0.0;
  react = 0.0;

#ifdef USE_SNES
  X.get_data(xaux.ref());
#endif

  MixedFiniteElement * fe = fespace.createFE();

  // Loop over all the elements
  for(int i=0; i<nelem; i++)
  {
    calc_elmatvec (i, fe, elmat_const, elmat_sigma, elvec, dnums);

    fespace.get_element_dofs_u(i,dnums);

    // Residual vector and reaction vector
    for(int k=0;k<n;k++){
      if (ldgof[dnums[k]])
        r.add(dnums[k], -elvec(k)); // -R
      else
        react.add(dnums[k], elvec(k));
    }
  }
  r.assemble();
  delete fe;
}

void UpdatedLagrangian::form_jacobian()
{
  int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();

  Matrix elmat_const(n,n);
  Matrix elmat_sigma(n,n);
  Vector elvec(n);
  std::vector<int> dnums;

  K = 0.0;

#ifdef USE_SNES
  X.get_data(xaux.ref());
#endif

  MixedFiniteElement * fe = fespace.createFE();

  // Loop over all the elements
  for(int i=0; i<nelem; i++)
  {
    calc_elmatvec (i, fe, elmat_const, elmat_sigma, elvec, dnums);

    fespace.get_element_dofs_u(i,dnums);

    // K_const + K_sigma
    for(int k=0;k<n;k++)
      for(int l=0;l<n;l++)
        K.add(dnums[k],dnums[l], elmat_const(k,l)+elmat_sigma(k,l) );
  }

  K.assemble();
  delete fe;
}
#endif

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Monitor"
PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal norm, void* mctx)
{
  Mat J;
  Vec x, f;

  cout << "JACOBIAN MATRIX\n";
  SNESGetJacobian(snes, &J, &J, PETSC_NULL, PETSC_NULL);
  MatView(J,PETSC_VIEWER_STDOUT_WORLD);

  cout << "FUNCTION VECTOR\n";
  SNESGetFunction(snes, &f, PETSC_NULL, PETSC_NULL);
  VecView(f,PETSC_VIEWER_STDOUT_WORLD);

  cout << "SOLUTION\n";
  SNESGetSolution(snes,&x);
  VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  cout << "SOLUTION UPDATE\n";
  SNESGetSolutionUpdate(snes,&x);
  VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *nlobj)
{
  PetscScalar *xx, *ff;
  PetscErrorCode ierr;
  UpdatedLagrangian * nlp = (UpdatedLagrangian *) nlobj;

  // Copy x to big X within NonlinearElasticity problem
  // in order to compute the function correctly
  ierr = VecCopy(x,nlp->X.vec());CHKERRQ(ierr);

  ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

  // Compute function
  nlp->form_function();

  PETSc::Vector & r = nlp->get_residual_vector();

  // ****************************************************************
  // Allocating memory everytime here -> IMPROVE!!!
  // ****************************************************************
  Vector uu(nlp->fespace.get_ndofs());
  r.get_data(uu.ref());
  // ****************************************************************

  for(int i=0; i<nlp->fespace.get_ndofs(); i++)
    ff[i] = uu[i];

  // Restore vectors
  ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat *jac,Mat *B,MatStructure *flag, void *nlobj)
{
  PetscErrorCode ierr;
  UpdatedLagrangian * nlp = (UpdatedLagrangian *) nlobj;
  PETSc::Matrix & A = nlp->get_jacobian_matrix();

  ierr = VecCopy(x,nlp->X.vec());CHKERRQ(ierr);

  nlp->form_jacobian();
  nlp->apply_boundary();

  ierr = MatCopy(A.mat(), *jac, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  //*jac = nlp->get_jacobian_matrix().mat();

  *flag = SAME_NONZERO_PATTERN;

  ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*jac != *B){
    ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  //ierr = MatView(*jac,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  return 0;
}

#endif

