#include "tensor.hpp"

void Tensor4::get_matrix(arma::mat & D)
{
  Tensor4 & T = (*this);

  if(D.n_rows == 3 && D.n_cols == 3) // 2D case
  {
    D(0,0) = 2*T(0,0,0,0);
    D(0,1) = 2*T(0,0,1,1);
    D(0,2) = T(0,0,0,1) + T(0,0,1,0);

    D(1,0) = 2*T(1,1,0,0);
    D(1,1) = 2*T(1,1,1,1);
    D(1,2) = T(1,1,0,1) + T(1,1,1,0);

    D(2,0) = D(2,0);
    D(2,1) = D(2,1);
    D(2,2) = T(0,1,0,1) + T(0,1,1,0);

    D = 0.5*D;
  }
  else if(D.n_rows == 6 && D.n_cols == 6) // 3D case  
  {
    D(0,0) = T(0,0,0,0);
    D(0,1) = T(0,0,1,1);
    D(0,2) = T(0,0,2,2);
    D(0,3) = 0.5 * ( T(0,0,0,1) + T(0,0,1,0) );
    D(0,4) = 0.5 * ( T(0,0,1,2) + T(0,0,2,1) );
    D(0,5) = 0.5 * ( T(0,0,0,2) + T(0,0,2,0) );

    D(1,0) = T(1,1,0,0);
    D(1,1) = T(1,1,1,1);
    D(1,2) = T(1,1,2,2);
    D(1,3) = 0.5 * ( T(1,1,0,1) + T(1,1,1,0) );
    D(1,4) = 0.5 * ( T(1,1,1,2) + T(1,1,2,1) );
    D(1,5) = 0.5 * ( T(1,1,0,2) + T(1,1,2,0) );

    D(2,0) = T(2,2,0,0);
    D(2,1) = T(2,2,1,1);
    D(2,2) = T(2,2,2,2);
    D(2,3) = 0.5 * ( T(2,2,0,1) + T(2,2,1,0) );
    D(2,4) = 0.5 * ( T(2,2,1,2) + T(2,2,2,1) );
    D(2,5) = 0.5 * ( T(2,2,0,2) + T(2,2,2,0) );

    D(3,0) = D(0,3);
    D(3,1) = D(1,3);
    D(3,2) = D(2,3);
    D(3,3) = 0.5 * ( T(0,1,0,1) + T(0,1,1,0) );
    D(3,4) = 0.5 * ( T(0,1,1,2) + T(0,1,2,1) );
    D(3,5) = 0.5 * ( T(0,1,0,2) + T(0,1,2,0) );

    D(4,0) = D(0,4);
    D(4,1) = D(1,4);
    D(4,2) = D(2,4);
    D(4,3) = 0.5 * ( T(1,2,0,1) + T(1,2,1,0) );
    D(4,4) = 0.5 * ( T(1,2,1,2) + T(1,2,2,1) );
    D(4,5) = 0.5 * ( T(1,2,0,2) + T(1,2,2,0) );

    D(5,0) = D(0,5);
    D(5,1) = D(1,5);
    D(5,2) = D(2,5);
    D(5,3) = 0.5 * ( T(0,2,0,1) + T(0,2,1,0) );
    D(5,4) = 0.5 * ( T(0,2,1,2) + T(0,2,2,1) );
    D(5,5) = 0.5 * ( T(0,2,0,2) + T(0,2,2,0) );

    //D = 0.5*D;
  }
}

void Tensor4::view()
{
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)	
      for(int k=0; k<3; k++)
	for(int l=0; l<3; l++)
	{
	  std::cout << " i=" << i << " j=" << j;
	  std::cout << " k=" << k << " l=" << l;
	  std::cout << " ---> ";
	  std::cout << data[get_array_index(i,j,k,l)];
	  std::cout << std::endl;
	}
}

void Tensor4::zero()
{
  for (int i=0; i<81; i++) data[i] = 0.0;
}


// Operations ------------------------------------------------------------------

Tensor4 tensor_product(const arma::mat33 & A, const arma::mat33 & B)
{
  Tensor4 T;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)	
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          T(i,j,k,l) = A(i,j) * B(k,l);
  return T;
}

Tensor4 unit_tensor()
{
  arma::mat33 delta = arma::eye(3,3);
  Tensor4 I;
  I.zero();
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)	
      for(int k=0; k<3; k++)
	    for(int l=0; l<3; l++)
	      I(i,j,k,l) = delta(i,k) * delta(j,l);

  return I;
}

// Tensor4s --------------------------------------------------------------------

void Tensor4s::view()
{
  using namespace std;
  cout << data[0]  << " " << data[1]  << " " << data[3]  << " "
       << data[6]  << " " << data[10] << " " << data[15] << endl ;
  cout << data[1]  << " " << data[2]  << " " << data[4]  << " "
       << data[7]  << " " << data[11] << " " << data[16] << endl ;
  cout << data[3]  << " " << data[4]  << " " << data[5]  << " "
       << data[8]  << " " << data[12] << " " << data[17] << endl ;
  cout << data[6]  << " " << data[7]  << " " << data[8]  << " "
       << data[9]  << " " << data[13] << " " << data[18] << endl ;
  cout << data[10] << " " << data[11] << " " << data[12] << " "
       << data[13] << " " << data[14] << " " << data[19] << endl ;
  cout << data[15] << " " << data[16] << " " << data[17] << " "
       << data[18] << " " << data[19] << " " << data[20] << endl ;
  cout << endl;
}
