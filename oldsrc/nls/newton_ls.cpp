#include "newton_ls.hpp"

#include "pdes/nonlinear_elasticity.hpp"

NewtonLineSearch::NewtonLineSearch(NonlinearProblem * p)
  : Newton(p), num_search(0), max_search(5), search(0.5)
{
  // nothing
}

NewtonLineSearch::~NewtonLineSearch()
{
  // nothing
}

void NewtonLineSearch::line_search(double & eta0, double & eta,
                                   double & rtu0, double & rtu)
{
  double rtu1, alfa, q;

  rtu = abs( u->dot(*r) );
  if (rtu == 0.0) return;

  rtu1 = (rtu - rtu0 * (1 - eta))/(eta*eta);
  alfa = rtu0/rtu1;
  eta0 = eta;

  // take different actions, depending on the value of alfa
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

int NewtonLineSearch::solve()
{
  int step = 0;
  bool done = false;
  std::pair<int,double> itres;

  nits = 0;
  num_search = 0;

  nlp->jacobian(*K);
  nlp->evaluate(*r);

  while (!done && (nits < nmax))
  {
    step = step + 1;
    nits = nits + 1;
    cout << " " << nits << " Newton-LS step";

    // solve linear system Ku = -r
    timer.enter("Linear Solver");
    itres = ls->solve(*K, *u, *r);
    timer.leave();
    //itres = ls->solve_3_FieldSplit(*K, *u, *r);

    cout << " KSP its " << itres.first;

    // compute initial residual (u != 0)
    if(nits==1) done = nlp->converged(*u);

    // line search loop
    double eta0, eta, rtu0, rtu;
    eta0 = 0.0;
    eta  = 1.0;
    rtu0 = abs( u->dot(*r) ); // energy
    rtu  = rtu0 * search * 2;
    num_search = 0;

    while( (abs(rtu) > abs(rtu0 * search)) && (num_search < max_search) )
    {
      num_search++;

      // update solution x = x + u
      nlp->update(*u, eta-eta0);

      // re-evaluate
      nlp->evaluate(*r);

      // perform line search
      line_search(eta0, eta, rtu0, rtu);
    }

    // check for convergence
    done = nlp->converged(*u);
    if(done) break;

    // recompute tangent matrix and residual
    nlp->jacobian(*K);
    //nlp->evaluate(*r);

    // TODO: pressure
    //NonlinearElasticity * nelas = dynamic_cast<NonlinearElasticity*>(nlp);
    ///nelas->assemble_pressure();
    //Derived<int> *d = dynamic_cast<Derived<int> *>(b);
  }

  if (nits >= nmax)
  {
    ls->converged_reason();
    error("newton solver did not converge");
  }

  return nits;
}

