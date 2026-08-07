/*
#include <iostream>
#include "../vector.hpp"
#include "../matrix.hpp"
#include "../../util/datatype.hpp"

int main()
{

  int k = 300;
  
  Real **mt = new double*[k];
  for(int i=0; i<k; i++) mt[i] = new double[k];

  for(int i=0; i<k; i++)
    for(int j=0; j<k; j++) 
      mt[i][j] = 2*i*j + i + 10;

  
  Matrix m1(k,k,mt);
  Matrix m2(k,k,5);
  
  
  Matrix m3(k,k);
  for(int i=0; i<k; i++)
    for(int j=0; j<k; j++) 
      m3[i][j] = 1/(2*i + j + 5.7);

  m3 += m1 + m2;
  m1 -= m3;

  Vector v(k);
  for(int i=0; i<k; i++)
    v[i] = 5*i + 3;

  v = m3*v;

  cout << v;

}
*/
