#include "bfgs.hpp"
   
BFGSSolver::BFGSSolver(NonlinearProblem * p) : NonlinearSolver(p)
{
  maxups = 20;
  maxref = 15;
  cmax   = 1e5;
  
  // LS parameters
  ls_min = 0.01;
  ls_tol = 0.9;
  ls_mxiter = 5;
  
  // pointer to linear solver
  parsolver = 0;
}

BFGSSolver::~BFGSSolver()
{
  delete K0;
}

void BFGSSolver::init()
{
  // nothing
  // init == pre_process
}
  
void BFGSSolver::init(int n, PardisoSolver * ps)
{  
  neq  = n;
  nups = 0;
  parsolver = ps; 
  
  // allocate memory
  V.resize(neq, maxups); //column-major
  W.resize(neq, maxups);
  
  D.resize(neq);
  G.resize(neq);
  H.resize(neq);

  ui.resize(neq); ui.zeros();
  r0.resize(neq); r0.zeros();
  r1.resize(neq); r1.zeros();
  
  aux_r0.create(neq); aux_r0 = 0.0;
  aux_r1.create(neq); aux_r1 = 0.0;
  aux_ui.create(neq); aux_ui = 0.0;
      
  K0 = new petsc::Matrix();
  K0->create(neq,neq,120); 
}

double BFGSSolver::line_search(double s)
{
	double smin = s;
	double a, A, B, D;
	double e0, e1, er;

	// max nr of line search iterations
	int nmax = ls_mxiter;
	int n = 0;

	// initial energy
	e0 = fabs( arma::dot(ui,r0) );

	double rmin = fabs(e0);  

	// ul = ls*ui
  int size = ui.size();
  petsc::Vector ul;
  ul.create(size);
    
	do
	{
		// update geometry    
    ul.copy_values(size, aux_ui);
    ul.scale(s);   
    nlp->update(ul, s);
		
		// calculate residual at this point
		nlp->evaluate(aux_r1);

		// make sure we are still in a valid range
		if (s < ls_min) 
		{			
			s = 0.5;
			// reupdate  
      ul.copy_values(size, aux_ui);
      ul.scale(s);
      nlp->update(ul, s);
			
			// recalculate residual at this point
			nlp->evaluate(aux_r1);
			break;
		}

		// calculate energies
		e1 = fabs( arma::dot(ui, r1) );

		if ((n==0) || (fabs(e1) < rmin))
		{
			smin = s;
			rmin = fabs(e1);
		}

		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
    er = (fabs(e1) < 1.e-20) ? 0 : fabs(e1/e0);
    
    cout << scientific << setprecision(4);
    cout << "   energy0 = " << e0 << "   energy1 = " << e1;
    cout << "   abs(e1/e0) = " << er << endl;

		if (er > ls_tol)
		{
			// calculate the line search step
			a = e0/e1;
			A = 1 + a*(s-1);
			B = a*(s*s);
			D = B*B - 4*A*B;

			// update the line step
			if (D >= 0) 
			{
				s = (B + sqrt(D))/(2*A);
				if (s < 0) s = (B - sqrt(D))/(2*A);
				if (s < 0) s = 0;
			}
			else 
			{
				s = 0.5*B/A;
			}

      cout << "  LS step = " << s << endl;
			n++;
		}    
	}
	while ((er > ls_tol) && (n < nmax));
  
  cout << "NUPS = " << nups << endl;

	if (n >= nmax)
	{
		// max nr of iterations reached
    // line step with the smallest energy
		s = smin;
    ul.copy_values(size, aux_ui);
    ul.scale(s);
    nlp->update(ul, s);
    nlp->evaluate(aux_r1);
	}

	return s;
}

void BFGSSolver::solve_eqs(double * x, double * b)
{
  int i,j;
  double *vi, *wi, vr, wr;
  static arma::vec tmp(neq);
  for(i = 0; i < neq; i++) tmp(i) = b[i];
  	
	for (i=nups-1; i>=0; i--)
	{
		vi = V.colptr(i);
		wi = W.colptr(i);

		wr = 0;
		for (j=0; j<neq; j++) wr     += wi[j]*tmp[j];
		for (j=0; j<neq; j++) tmp[j] += vi[j]*wr;
	}

	// backsubstitution
  parsolver->back_solve(x, tmp.memptr());

	// loop again over all update vectors
	for (i=0; i<nups; i++)
	{
		vi = V.colptr(i);
		wi = W.colptr(i);

		vr = 0;
		for (j=0; j<neq; ++j) vr   += vi[j]*x[j];
		for (j=0; j<neq; ++j) x[j] += wi[j]*vr;
	}

}

int BFGSSolver::solve()
{   
  static bool bsetup = true;
 
  nups = 0;
  
  aux_r0.place_array(r0.memptr());
  aux_r1.place_array(r1.memptr());
  aux_ui.place_array(ui.memptr());
 
  // compute and factor the jacobian matrix
	nlp->jacobian(*K0);	
  
  // create the initial jacobian matrix
  if(bsetup)
  {
    parsolver->set_matrix(K0);  
    parsolver->init();
    bsetup = false;
  }
  else parsolver->update_matrix(K0);
  
	parsolver->factor();

	// create the initial RHS vector
	nlp->evaluate(aux_r0);
  
	do
	{    
		solve_eqs(ui.memptr(), r0.memptr());
    
		// update solution
		if (ls_tol > 0) line_search(1.0);
		else nlp->update(aux_ui, 1.0); // ATENCAO aqui...checar esse s=1   
    
		// check for convergence
		bool bconv = true;
		for (int i=0; i<neq; i++) 
			if (fabs(ui(i)) > 1e-3) { bconv = false; break;	}
    
    //bconv = nlp->converged(aux_ui);
       
		if (bconv == false)
		{
			nlp->evaluate(aux_r1); 

			update(1.0, ui, r0, r1);       
			r0 = r1;
		}
		else break;
	}
	while (true);
  
  aux_r0.reset_array();
  aux_r1.reset_array();
  aux_ui.reset_array();      
  cout << " BFGS solve done" << endl;  
  return 0;
}

bool BFGSSolver::update(double s, arma::vec & UI, arma::vec & R0, arma::vec & R1)
{  
  double dg, dh, dgi, c, r;
  //double *vn, *wn;
  
  neq = UI.size();
  
  D = s * UI;  //xk = xk-1 + s ui ---> D = xk - xk-1 = s ui
  G = R0 - R1;
  H = R0 * s;
  
  dg  = arma::dot(D,G);
  dh  = arma::dot(D,H);
  dgi = 1./dg;
  r   = dg/dh;
  
  // condition number and check if is less than max
  c = sqrt(fabs(r));
  if(c > cmax) 
  { 
    cout << "BFGS: condition number of update matrix too large" << endl;
    return false;
  }
    
  if(nups > maxups) throw std::runtime_error("maxups reached");
  
  // compute new 
  V.col(nups) = -H*c - G;
  W.col(nups) = D * dgi;
   
//  vn = V.colptr(nups);
//  wn = W.colptr(nups);
//  for(i = 0; i < neq; i++)
//  {
//    vn[i] = -H(i)*c - G(i);
//    wn[i] = D(i)*dgi;
//  }
    
  // update performed, so increment counter
  nups = nups + 1;
  if(nups >= maxups)
  {
    cout << "Recomputing K matrix" << endl;
    nlp->jacobian(*K0);	
  }
  nups = nups % maxups;
  
  return true;
}