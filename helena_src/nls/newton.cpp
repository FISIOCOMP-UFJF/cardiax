#include "newton.hpp"
//#include "../pdes/nonlinear_elasticity.hpp"

Newton::Newton(NonlinearProblem * p) 
  : NonlinearSolver(p), nits(0), nmax(50)
{
  // nothing
}

Newton::~Newton()
{
  delete K;
  delete u;
  delete r;  
}

void Newton::init()
{
  assert(nlp != NULL);
   
  int size = nlp->size();
  cout << "Size of the problem " << size << endl;
  
  K = new petsc::Matrix();
  u = new petsc::Vector();
  r = new petsc::Vector();
  K->create(size,size,120);
  u->create(size);
  r->create(size);
  cout << "Matrices and vectors creation done" << endl;
}

void Newton::ViewMatrix()
{
  ls->MatlabView(*K);
}

int Newton::solve()
{
  int step = 0;
  bool done = false;
  std::pair<int,double> itres; 

  nits = 0;

  nlp->jacobian(*K);
  nlp->evaluate(*r); 
  
  while (!done && (nits < nmax))
  {    
    step = step + 1;
    nits = nits + 1;
    cout << " " << nits << " Newton step ";

    // solve linear system Ku = -r
    itres = ls->solve(*K, *u, *r);

       
    // compute initial residual (u != 0)
    if(nits==1)
      done = nlp->converged(*u);
    
    // update solution x = x + u
    nlp->update(*u, 1);    
    
    // recompute tangent matrix and residual
    nlp->evaluate(*r);
    nlp->jacobian(*K);
       
    // check for convergence
    done = nlp->converged(*u);
  }  

  if (nits >= nmax)
  {
    ls->converged_reason();
    error("newton solver did not converge");
  }    
  
  return nits;
}
