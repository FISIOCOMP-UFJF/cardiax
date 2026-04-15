#include "total_lagrangian_snes.hpp"

extern PetscErrorCode FormJacobian1(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction1(SNES,Vec,Vec,void*);
extern PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal norm, void* mctx);
extern PetscErrorCode MyConvergenceTest(SNES snes, PetscInt it, PetscReal xnorm,
					PetscReal snorm, PetscReal fnorm, 
					SNESConvergedReason *reason, void *ctx);

void TotalLagrangianSNES::update_coordinates(PetscScalar *xx)
{
  uint np = x.size();
  uint nd = msh.get_ndim();

  for(uint i=0; i<np; i++)
    for(uint j=0; j<nd; j++)
    {
      x[i][j] = xx[i+j*np];
    }

  /*
  // DEBUG ***
  cout << "COORDINATES updateFunc" << endl;
  for(int i=0; i<np; i++)
  {
    for(int d=0; d<nd; d++)	
    {   
      cout << x[i][d] << " ";
    }
    cout << endl;
  }
  */
}

void TotalLagrangianSNES::init_resid()
{
  r = 0;
  react = 0.0;

  // TL using SNES uses:  -lc.load() * fext
  for(int i=0; i<fespace.get_ndofs(); i++)
    r.set(i, fext0(i) - lc.load() * fext(i) );
}

void TotalLagrangianSNES::init_stiff()
{
  K = 0.0;
  u = 0.0;
}

void TotalLagrangianSNES::func()
{
  int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();
  arma::vec Re(n);
  arma::mat Ke(n,n);
  std::vector<int> dnums;

  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type());

  for(int i=0; i<nelem; i++)
  {
    calc_elmatvec (i, fe, qd, Ke, Re);
    fespace.get_element_dofs_u (i,dnums);

    for(int j=0; j<n; j++)
    {
      int k = dnums[j];
      if (ldgof[k])
        r.add(k, Re(j));
      else
	react.add(k, Re(j));
    }
  }

  r.assemble();

  delete qd;
  delete fe;
}

void TotalLagrangianSNES::jacob()
{
 int n = msh.get_ndim() * msh.get_nen();
  int nelem = msh.get_n_elements();
  arma::vec Re(n);
  arma::mat Ke(n,n);
  std::vector<int> dnums;

  MxFE * fe = fespace.createFE();
  Quadrature * qd = Quadrature::create(2*fe->get_order_u()-2, fe->get_type());

  for(int i=0; i<nelem; i++)
  {
    calc_elmatvec (i, fe, qd, Ke, Re);
    fespace.get_element_dofs_u(i,dnums);

    for(int k=0; k<n; k++)
      for(int l=0; l<n; l++)
        K.add(dnums[k], dnums[l], Ke(k,l));
  }

  K.assemble();

  delete qd;
  delete fe;  
}

void TotalLagrangianSNES::solve()
{
  static int cont = 1;
  int num_prescribed = get_num_nz_prescribed();
  int N = fespace.get_ndofs();

  // --------------------------------------------------------------------------
  // Setup SNES
  // --------------------------------------------------------------------------

  SNES snes;
  KSP ksp;
  PC pc;
  Vec x;
  //Vec residual;
  PetscInt its;
  SNESConvergedReason reason;

  // My stuff
  uint np = msh.get_n_points();
  uint nd = msh.get_ndim();
  arma::vec vvaux(nd*np);
  petsc::Vector vaux;
  vaux.create(nd*np);

  SNESCreate(PETSC_COMM_WORLD,&snes); 
  
  VecCreate(PETSC_COMM_WORLD,&x);     
  VecSetSizes(x,PETSC_DECIDE, N);     
  VecSetFromOptions(x);            
  //VecDuplicate(x,&residual);       
  
  x_aux.create(N);

  // Set function evaluation routine and vector.
  SNESSetFunction(snes, this->r.vec(), FormFunction1, this);

  // Set Jacobian matrix data structure and Jacobian evaluation routine
  SNESSetJacobian(snes, this->K.mat(), this->K.mat(), FormJacobian1, this);

  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCNONE);

  SNESMonitorSet(snes,Monitor,this,0);
  SNESSetFromOptions(snes);

  // No line search
  //SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL);

  // Print tolerances
  PetscScalar atol, rtol, stol;    
  PetscInt maxit, maxf;
  SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf);
  cout << "SNES Tolerances" << endl;
  cout << " Relative tolerance : " << rtol << endl;
  cout << " Absolute toleracen : " << atol << endl;
  cout << " Step tolerance     : " << stol << endl;

  cout << "Solving nonlinear problem using SNES" << endl;

  if(output_step) 
    output_vtk(cont,0);
   
  pre_computation(); 
  body_forces();
  assemble_traction();

  // --------------------------------------------------------------------------
  // Load increment loop
  // --------------------------------------------------------------------------

  while ( lc.has_load() )
  {
    lc.update();

    evaluate_forces(r);

    cout << endl << "Load increment " << lc.increment() << endl;

    if (num_prescribed > 0 || neumann_map.size() > 0)
    {
      prescribe_displacements();
      init_resid_stiff();
      //assemble_pressure();

      func();
      jacob();
      apply_boundary(K);      
    }

    // --------------------------------------------------------------------------
    // SNES
    // --------------------------------------------------------------------------  

    cout << "STARTING SNES" << endl;

    K.assemble();
    r.assemble(); 

    // Initial guess    
    for(uint i=0; i<np; i++)
    {
      for(uint j=0; j<nd; j++)
      {
	vvaux(i+j*np) = this->x[i][j];
      }
    }
    vaux.set_data(vvaux.memptr());
    VecCopy(vaux.vec(), x_aux.vec());

    // Change TOL
    //SNESSetTolerances(snes,atol,rtol,stol,PETSC_DEFAULT,PETSC_DEFAULT);

    // Solve
    SNESSolve(snes, PETSC_NULL, x_aux.vec());
    SNESGetIterationNumber(snes, &its);
    SNESGetConvergedReason(snes, &reason);
    
    PetscBool err=PETSC_TRUE;
    SNESGetErrorIfNotConverged(snes,&err);   
    cout << " [SNESConvergedReason] " << SNESConvergedReasons[reason] << endl;

    if(output_step)
      output_vtk(cont,lc.increment());

  } // end of load increment loop

  cont  += 1;
  fext0 += fext;
  
  // Clean 
  SNESDestroy(&snes);
}

// ----------------------------------------------------------------------------
// Callbacks for SNES 
// ----------------------------------------------------------------------------

PetscErrorCode FormFunction1(SNES snes, Vec x, Vec f, void *ctx)
{
  PetscScalar *ff, *xarray;
  PetscErrorCode ierr;

  // Get application context
  TotalLagrangianSNES *user = (TotalLagrangianSNES*) ctx;

  ierr = VecGetArray(x,&xarray); CHKERRQ(ierr);
  ierr = VecGetArray(f,&ff); CHKERRQ(ierr);

  // Copy x from SNES to x within TL problem
  user->update_coordinates(xarray);
  user->init_resid();
  user->func();
  
  petsc::Vector & rr = user->residual();
  arma::vec uu(user->fespace.get_ndofs());
  rr.get_data(uu.memptr());

  for(int i=0; i<user->fespace.get_ndofs(); i++)
    ff[i] = uu[i];

  // Restore vectors
  ierr = VecRestoreArray(f,&ff); CHKERRQ(ierr);
  ierr = VecRestoreArray(x,&xarray); CHKERRQ(ierr); 
  
  return 0;
}


PetscErrorCode FormJacobian1(SNES snes, Vec x, Mat * jac, Mat * B, 
			     MatStructure *flag, void * ctx)
{
  PetscScalar *xarray;
  PetscErrorCode ierr;
  
  // Get application context
  TotalLagrangianSNES  *user = (TotalLagrangianSNES*)ctx;

  VecGetArray(x,&xarray);

  user->update_coordinates(xarray);
  user->init_stiff();
  user->jacob();
  user->apply_boundary(user->get_K());

  // Preconditioning matrix has the same nonzero structure during
  // successive linear solves
  *flag = SAME_NONZERO_PATTERN;

  // Assemble matrix
  ierr = MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  if (*jac != *B){
    ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  }

  VecRestoreArray(x,&xarray);

  return 0;
}

PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal norm, void* mctx)
{
  //Mat J;
  //Vec x, f;

  //cout << " [SNES Monitor]";
  //cout << " [step " << setw(2) << its;
  //cout << ", fnorm = " << scientific << norm << "]" << endl;

  //cout << "JACOBIAN" << endl;
  //SNESGetJacobian(snes, &J, &J, PETSC_NULL, PETSC_NULL);
  //MatView(J,PETSC_VIEWER_STDOUT_WORLD);
  
  //cout << "RESIDUAL" << endl;
  //SNESGetFunction(snes, &f, PETSC_NULL, PETSC_NULL);
  //VecView(f,PETSC_VIEWER_STDOUT_WORLD);

  //cout << "SOLUTION" << endl;
  //SNESGetSolution(snes,&x);
  //VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  //cout << "SOLUTION UPDATE" << endl;
  //SNESGetSolutionUpdate(snes,&x);
  //VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  //PetscScalar *xarray;
  //VecGetArray(x,&xarray);
  //for(int i=36; i<60; i++)
  //for(int i=0; i<395; i++)
  //  cout << i << " " << scientific << xarray[i] << endl;
  //VecRestoreArray(x,&xarray);

  return 0;
}

PetscErrorCode MyConvergenceTest(SNES snes, PetscInt it, PetscReal xnorm,
				 PetscReal snorm, PetscReal fnorm, 
				 SNESConvergedReason *reason, void *ctx)
{
  //AppCtx        *user = (AppCtx *) ctx;
  //Parameter     *param = user->param;
  //KSP            ksp;
  
  /*
    if (param->interrupted) {
     param->interrupted = PETSC_FALSE;
     PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: exiting SNES solve. \n");
     *reason = SNES_CONVERGED_FNORM_ABS;
     return(0);
   } else if (param->toggle_kspmon) {
     param->toggle_kspmon = PETSC_FALSE;
     SNESGetKSP(snes, &ksp);
     if (param->kspmon) {
       KSPMonitorCancel(ksp);
       param->kspmon = PETSC_FALSE;
       PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: deactivating ksp singular value monitor. \n");
     } else {
       KSPMonitorSet(ksp,KSPMonitorSingularValue,PETSC_NULL,PETSC_NULL);
       param->kspmon = PETSC_TRUE;
       PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: activating ksp singular value monitor. \n");
     }
   }
*/
  cout << "RNORM = " << fnorm << endl;
  cout << "SNORM = " << snorm << endl;
   
  //*reason = SNES_DIVERGED_MAX_IT;
  //PetscFunctionReturn(SNESDefaultConverged(snes,it,xnorm,snorm,fnorm,reason,ctx));

  return 0;
}
