/*
#include <iostream>
#include <time.h>
#include "vec.hpp"
#include "datatype.h"

void add(uint n, Real *z, Real *x, Real *y)
{
  for(uint i=0; i<n; i++)
    z[i] = x[i] + y[i];
}

void saxpy(uint n, Real *z, Real a, Real *x, Real *y)
{
  for(uint i=0; i<n; i++)
    z[i] = a*x[i] + y[i];
}

int main()
{
  int k = 50000;
  int num = 10000;

  Real *p = new Real[k];
  Real *z = new Real[k];
  Real *x = new Real[k];
  Real *y = new Real[k];
  
  for(uint i=0; i<k; i++)
  {
    p[i] = i;
    z[i] = 0.0;
    x[i] = i;
    y[i] = i+1;
  }

  Vector v(k);
  v=0.1;
  Vector w(k);
  w = 0.6;

  Vector u(k);
  u = 2.0*w;

  time_t tm0 = time(0);
  clock_t ck0 = clock();

  for(uint i=0; i<num; i++)
    v = 3.14 * v + w;

  time_t tm1 = time(0);
  clock_t ck1 = clock();
  
  cout << "Using expression templates\n";
  cout << " wall time = " << difftime(tm1,tm0) << " seconds.\n";
  cout << " CPU time  = " << double(ck1-ck0)/CLOCKS_PER_SEC 
       << " seconds.\n";

  for(uint i=0; i<num; i++)
    //add(k,z,x,y);
    saxpy(k,z,3.14,x,y);
  //z[i] = 3.14 * x[i];

  time_t tm2 = time(0);
  clock_t ck2 = clock();
  cout << "Using plain C loop\n";
  cout << " wall time = " << difftime(tm2,tm1) << " seconds.\n";
  cout << " CPU time  = " << double(ck2-ck1)/CLOCKS_PER_SEC 
       << " seconds.\n";

  //  Benchmark for defereed-evaluation, old version!
  //  SAXPY
  //  Output of code compiled with -O3 flag.

  //  deferred-evaluation
  //  wall time = 4 seconds.
  //  CPU time  = 4.01 seconds.
  //  traditional saxpy operation
  //  wall time = 4 seconds.
  //  CPU time  = 4.03 seconds.

  //  Output of code compiled without optimizations.

  //  deferred-evaluation
  //  wall time = 17 seconds.
  //  CPU time  = 16.84 seconds.
  //  traditional saxpy operation
  //  wall time = 6 seconds.
  //  CPU time  = 5.83 seconds.
  //

  delete [] p;
  delete [] x;
  delete [] y;
  delete [] z;

  return 0;
}

*/
