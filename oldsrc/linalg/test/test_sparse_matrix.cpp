#include <iostream>
#include "spmat.hpp"

int main(int argc, char*argv[])
{
  SparseMatrix A;
  A.resize(10,10,10);
  
  for(int i=0; i<10; i++)
    for(int j=0; j<3; j++)
      A(i,j) = i+j + 100;

  std::cout << "Before compressing...\n";
  std::cout << A;
  A.Info();

  A.compress();
  
  std::cout << endl;
  std::cout << "After compressing...\n";
  std::cout << A;
  A.Info();
  

}
