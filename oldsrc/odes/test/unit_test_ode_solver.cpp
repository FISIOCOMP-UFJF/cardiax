#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <typeinfo>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "../../util/datatype.hpp"
#include "cellmodel.h"
#include "cellmodel_factory.h"
#include "fitz_hugh_nagumo.h"
#include "ten_tusscher.h"
#include "luo_rudy.h"

using namespace std;

Real compute_error(uint N, Real *u, Real *v)
{
  Real tmp = 0.0;
  Real err = 0.0;
  Real den = 0.0;
  ifstream apdata("ap");
  ifstream exdata("exact");

  // read and close files
  for(uint i=0; i<N; i++)
  {
    apdata >> tmp >> u[i];    
    exdata >> tmp >> v[i];
  }
  
  apdata.close();
  exdata.close();

  // compute discrete L-2 error
  for(uint i=0; i<N; i++)
  {
    err += (u[i]-v[i])*(u[i]-v[i]);
    den += (v[i]*v[i]);
  }
  
  cout << "\n L2 Error: " << sqrt(err)/sqrt(den) << endl;

  return sqrt(err)/sqrt(den);
}

BOOST_AUTO_TEST_SUITE(test_suite_ode_solver)

BOOST_AUTO_TEST_CASE(create_test)
{
  cout << "\n*** Running ODESolver unit test " << endl;

  // tolerance
  const Real tol = 0.01;

  // used variables
  Real T = 5.0;
  Real step = 0.001;
  Real erro1, erro2;
  uint N = (uint)T/step;
  Real *exact = new Real[N];
  Real *aprox = new Real[N];
  
  //
  // Test 1 with Forward Euler
  // Compute the relative error between exact and approximated
  // solution and check to see if it is smaller than tolerance
  //  
  SimpleODE *sode1 = new SimpleODE();
  sode1->solver("ForwardEuler");
  sode1->setup(step,T,step);
  sode1->solve();
  sode1->dump_exact();
  erro1 = compute_error(N,aprox,exact);
  BOOST_CHECK_LE(erro1,tol);  

  //
  // Test 2 with RK4
  //
  SimpleODE *sode2 = new SimpleODE();
  sode2->solver("RungeKutta4");
  sode2->setup(step,T,step);
  sode2->solve();
  sode2->dump_exact();
  erro2 = compute_error(N,aprox,exact);
  BOOST_CHECK_LE(erro2,tol); 

  delete [] exact;
  delete [] aprox;

}

BOOST_AUTO_TEST_SUITE_END()
