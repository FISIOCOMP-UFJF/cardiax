#include "incompressible_material.hpp"

IncompressibleMaterial::IncompressibleMaterial(double K)
  : kappa(K)
{
  cout << "Creating Incompressible Material" << endl;
  cout << "Bulk modulus: " << K << endl;
  uncoupled = true;  
}

double IncompressibleMaterial::strain_energy(MaterialData * md, 
                                             const arma::mat &) const
{
  return 0;
}

double IncompressibleMaterial::active_strain_energy(MaterialData * md,
                                             const arma::mat &) const
{
  return 0;
}

void IncompressibleMaterial::add_pressure(double press, 
                                          arma::mat & sigma)
{
  for(int i=0; i<ndim; i++)
    sigma(i,i) += press;
}

void IncompressibleMaterial::cauchy_stress(MaterialData * md, 
                                           arma::mat & sigma) const
{
  // TODO
}

void IncompressibleMaterial::piola2_stress(MaterialData * md, 
                                           arma::mat & S) const
{
  // TODO
}

void IncompressibleMaterial::sp_elastensor(MaterialData * md, 
                                           arma::mat & D) const
{
  // TODO
}

void IncompressibleMaterial::mt_elastensor(MaterialData * md, 
                                           arma::mat & D) const
{
  // TODO
}

void IncompressibleMaterial::sp_volumetric_elastensor(const double pressure,
                                                      Tensor4 & A) const
{
  static const arma::mat delta = arma::eye(ndim,ndim);
  static const Tensor4 II = unit_tensor();

  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
	      for(int l=0; l<ndim; l++)
            A(i,j,k,l) += pressure * delta(i,j) * delta(k,l) - (2.0*pressure*II(i,j,k,l));

	        //A(i,j,k,l) += pressure * delta(i,j) * delta(k,l)
  	        //           -  pressure * delta(i,k) * delta(j,l)
	        //           -  pressure * delta(i,l) * delta(j,k);
}


void IncompressibleMaterial::copy_elastensor(Tensor4 B, Tensor4 & A) const
{
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
        for(int l=0; l<ndim; l++)
          A(i,j,k,l) += B(i,j,k,l);
}

void IncompressibleMaterial::copy_stress(arma::mat A, arma::mat & S) const
{
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
          S(i,j) += A(i,j);
}


