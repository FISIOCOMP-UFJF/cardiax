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
                                             const arma::mat & E) const
{
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();
  arma::mat I = arma::eye<arma::mat>(3,3);

  // right Cauchy-Green deformation tensor
  arma::mat C = 2*E + I;
  double J = sqrt(arma::det(C));

  arma::mat Cbar = pow(J,-(2.0/3.0)) * C;
  double I4f  = dot(f0, Cbar*f0);

  double psi = 0.5 * md->get_active_stress() * (I4f - 1.0);

  // Stabilisation of the velocity-dependent active tension (Regazzoni &
  // Quarteroni, CMAME 373 (2021) 113506). The scheme replaces Ta by
  //     Ta + Ka * (lambda - lambda_prev)
  // in the stress. Adding it here as the elastic energy of the attached
  // crossbridges,
  //     psi_stab = 0.5 * Ka * (lambda - lambda_prev)^2,
  // is equivalent at the level of the stress -- d(psi_stab)/dE gives exactly
  // Ka*(lambda - lambda_prev)*dlambda/dE -- but has one decisive advantage:
  // the stress and the elasticity tensor are obtained here by finite
  // differences of this function, so the stabilisation enters the TANGENT as
  // well, for free.
  //
  // That matters more than it might seem. Ka = Lfac*(Tref/dr)*A*(XS+XW) with
  // A = 10, so Ka runs roughly 20x larger than Ta itself and, in systole,
  // two to three orders of magnitude above the passive fibre stiffness --
  // which is precisely Regazzoni's Ka > Kp condition. Feeding a stiffness
  // that large to the residual while hiding it from the Jacobian turns
  // Newton into a modified-Newton iteration with a very poor operator: it
  // converges linearly at best, and a single overshoot inverts an element.
  //
  // Written as a quadratic in lambda the term is also a genuine potential,
  // so the tangent stays symmetric and positive semi-definite, which is what
  // the CG/AMG solver expects.
  const double Ka = md->get_active_stiffness();
  if (Ka > 0.0)
  {
    // Full I4f, not the isochoric one: the same definition used to build
    // lambda_prev in NonlinearElasticity::fiber_stretch_elements. Mixing the
    // two would leave a spurious offset of order J^(-1/3) inside the
    // increment.
    const double I4f_full = dot(f0, C*f0);
    if (I4f_full > 0.0)
    {
      const double lam = sqrt(I4f_full);
      const double dlam = lam - md->get_lambda_prev();
      psi += 0.5 * Ka * dlam * dlam;
    }
  }

  return psi;
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


