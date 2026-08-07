// cgTest.cpp 
// simple test program to demonstrate and test cg functionality

#include <iostream>

#include "spmat.hpp"
#include "vec.hpp"
#include "matops.hpp"
#include "cg.hpp" 

int main()
{
  using std::cout;
  using std::cin;
  using std::endl;

  // Solve -u''(x) = 1 on (0,1) with u(0)=0 du/dn(1) = 0  
  
  int nel = 10;
  double h = 1./nel; 
  int ndofs = nel; 
  
  SparseMatrix A(ndofs,ndofs); 
  Matrix elmat(2,2); 
  elmat(0,0) = 1; elmat(1,1) = 1; elmat(0,1) = -1; elmat(1,0) = -1; 
  
  cout << "elmat " << endl << elmat << endl; 
  cout << " A " << endl << A << endl; 
  
  for (int i=0;i<nel-1; i++) sys
  {
     for(int l = 0; l<2; l++) 
        for(int k = 0; k<2; k++)
     {
           A(i+l,i+k) += 1/h * elmat(l,k);
     }      
  }
  A(ndofs-1,ndofs-1) = 1/h; 
  A(0,0) += 1/h; // u(0)=0; 
 
  A.compress(); 
  
   Vector f(ndofs);
   f = h; 
   f(ndofs-1) = h/2; 

   cout << " A " << endl << A << endl; 
   cout << " f " << endl << f << endl; 
   
   

   Vector x(ndofs); 
   cg_solve(A,f,x,100,100,100,1e-8); 
   cout << " x " << x << endl;  
   
   Vector y(ndofs); 
   SparseMatrix invC(ndofs,ndofs,1);
   for(int i=0;i<ndofs;i++) invC(i,i) = 1/A(i,i);
   
   invC.compress(); 
   
   pcg_solve(A,invC,f,y,100,100,100,1e-8); 
   cout << " y " << y << endl;  
   
   return 0;
}


