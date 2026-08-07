#ifndef ARRAY3_H
#define ARRAY3_H

#include <cassert>

/** Auxiliary class for storing data in a 3 dimensional array.
*/
class Array3
{

private:

  // Dimensions
  int nx;  
  int ny;
  int nz;

  // Array to hold the data
  double * data;

  // Compute the index of the (i,j,k) entry in the array
  int get_array_index(int i, int j, int k)
  {
    assert(i<nx && j<ny && k<nz);
    return k*nx*ny + i*ny + j;
  }

public:

  // Default constructor
  Array3() : data(NULL)
  {}

  // Constructor
  Array3(int a, int b, int c) : nx(a), ny(b), nz(c), data(new double[a*b*c])
  {}

  // Destructor
  ~Array3()
  {
    delete [] data;
  }

  // Allocate memory for the data array
  void allocate(int a, int b, int c)
  {
    if (data==NULL) 
    {
      nx = a;
      ny = b;
      nz = c;
      data = new double[a*b*c];      
    }
    for(int i=0; i<nx*ny*nz; i++) data[i] = 0.0;
  }

  // Get the corresponding component of the array
  double& operator()(int i, int j, int k)
  {
    int index = get_array_index(i,j,k);
    return data[index];
  }

  // Print all data
  void view()
  {
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)	
	for(int k=0; k<3; k++)
	  {
	    std::cout << " (" << i << "," << j << "," << k << ") = ";
	    std::cout << data[get_array_index(i,j,k)];
	    std::cout << std::endl;
	  }
  }


  // Set all components of the array to zero.
  void zero()
  {
    for (int i=0; i<27; i++) data[i] = 0.0;
  }

};


#endif
