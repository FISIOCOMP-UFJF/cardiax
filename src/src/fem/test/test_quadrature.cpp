#include <iostream>
#include "../quadrature.hpp"
#include "../../linalg/linalg.h"

inline double func1(double x){
  return 2.0*x + 4.0;
}

inline double func2(double x){
  return 2.0*x*x + 4*x + 8;
}

int main()
{
  double xi, wi;

  // Test 1 -> R:5
  Quadrature1d qd1 = Quadrature1d(1);
  double r1 = 0.0;
  for(int i=0; i<qd1.get_num_ipoints(); i++){
    xi = qd1.get_point(i)(0);
    wi = qd1.get_weight(i);
    r1 = r1 + wi * func1(xi);
  }
  cout << "Integral of func1 is: " << r1 << endl;

  // Test 2 -> R:10+2/3
  Quadrature1d qd2 = Quadrature1d(3);
  double r2 = 0.0;
  for(int i=0; i<qd2.get_num_ipoints(); i++){
    xi = qd2.get_point(i)(0);
    wi = qd2.get_weight(i);
    r2 = r2 + wi * func2(xi);
  }
  cout << "Integral of func2 is: " << r2 << endl;
  
  // TODO: create tests for triangles and quadrilaterals

  return 0;
}
