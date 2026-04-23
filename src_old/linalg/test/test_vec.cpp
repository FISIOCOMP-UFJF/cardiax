#include <iostream>
#include <time.h>
#include "vec.hpp"
#include "vecops.hpp"

using namespace std;

template<class T, int N>
inline T dot3(Vec<T,N> & a, Vec<T,N> & b)
{
  T s=0.0;
  for(int i=0; i<3; i++) 
    s = s + a[i]*b[i];
  return s;
}

int main()
{
  Vec<float,3> v1;
  Vec<float,3> v2;  
  Vec<float,3> v3;

  v3 = 19.0;
  v2 = 1.0;
  
  v1 = v2 + v3;

  cout << "V1 ";
  cout << v1;

  cout << "V2 ";
  cout << v2;

  cout << "V3 ";
  cout << v3;

  // benchmark traditional dot and template metaprogram version
  float r1, r2; 
  unsigned int k = 90000000;

  time_t tm0 = time(0);
  clock_t ck0 = clock();

  for(unsigned int i=0; i<k; i++)
    r1 = dot3(v1,v2);
  
  time_t tm1 = time(0);
  clock_t ck1 = clock();

  for(unsigned int i=0; i<k; i++)
    r2 = dot(v1,v2);

  time_t tm2 = time(0);
  clock_t ck2 = clock();

  cout << "\ntraditional dot_product with loop\n";
  cout << " wall time = " << scientific << difftime(tm1,tm0) << " seconds.\n";
  cout << " CPU time  = " << scientific <<  double(ck1-ck0)/CLOCKS_PER_SEC 
       << " seconds.\n";

  cout << "\ntemplate metaprogram version of dot product\n";
  cout << " wall time = " << scientific << difftime(tm2,tm1) << " seconds.\n";
  cout << " CPU time  = " << scientific << double(ck2-ck1)/CLOCKS_PER_SEC 
       << " seconds.\n\n";


  cout << r1 << endl;
  cout << r2 << endl;


  return 0;
}
