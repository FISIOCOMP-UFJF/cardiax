#include "hyperelastic_material.hpp"
#include "material_data.hpp"
#include "compressible_neo_hookean.hpp"
#include "saint_venant_kirchoff.hpp"
#include "exponential.hpp"
#include "neo_hookean.hpp"
#include "holzapfel_ogden.hpp"
#include "guccione.hpp"

HyperelasticMaterial::HyperelasticMaterial(const std::vector<double> & prm)
  : uncoupled(false)
{
  for(uint i=0; i<prm.size(); i++)
    parameters.push_back(prm[i]);

  // initialize active stress
  active_stress = 0.0;
}

double HyperelasticMaterial::get_param(int i) const
{
  assert((uint)i < (uint)parameters.size());
  return parameters[i];
}

void HyperelasticMaterial::calc_fd_stress(int iel, MaterialData * md, arma::mat & S)					  
{
  const double eps = 1.0e-8;
  double psi_m, psi_p;
  arma::mat33 E, dij;

  E = md->lagrangian_strain();

  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
    {
      delta_tensor(i, j, eps, dij);
      psi_p = strain_energy(md, E + dij);
      psi_m = strain_energy(md, E - dij);
      S(i,j) = (psi_p - psi_m)/(2.*eps);
    }

  // copy symmetric entries
  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
      S(j,i) = S(i,j);

  set_active_stress(md,S,1);
}


void HyperelasticMaterial::set_active_stress(arma::mat & S, double loadFactor)
{
  S(0,0) += active_stress*loadFactor;
}

void HyperelasticMaterial::set_active_stress(MaterialData * md, arma::mat & S, double loadFactor)
{
  arma::vec3 f = md->fiber();
  //arma::vec3 s = md->sheet();
  //f = f / arma::norm(f,2);
  //S += active_stress*loadFactor*f*f.t();
  S += md->get_J()*active_stress*loadFactor*(f*f.t()); //+ 0.4*s*s.t());
  //S += active_stress*loadFactor*f*f.t();
}

void HyperelasticMaterial::active_stress_elastensor(int nincs, arma::vec3 fib,
                                                      Tensor4 & A) const
{
  arma::mat delta = fib*fib.t();
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
        for(int l=0; l<ndim; l++)
        {
          A(i,j,k,l) += -(delta_active_stress/nincs) * delta(i,j) * delta(k,l); 

        }

}

void HyperelasticMaterial::active_stress_elastensor(int nincs, MaterialData * md,
                                                    Tensor4 & A)
{
  arma::vec3 f = md->fiber();
  arma::mat delta = f*f.t();
  Tensor4 II = unit_tensor();

  //cout << "f: " << f << endl;
  //cout << "delta: " << delta <<endl;

  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
        for(int l=0; l<ndim; l++){
          A(i,j,k,l) += md->get_J()*(active_stress/nincs) * delta(i,j) * delta(k,l);
          //A(i,j,k,l) += (active_stress/nincs) * delta(i,j) * delta(k,l);
         //cout << 1000. * delta(i,j) * delta(k,l) << endl;
	//A(i,j,k,l) += 1000.*md->get_J()*(1./nincs) * delta(i,k) * delta(j,l);
	}
}

void HyperelasticMaterial::map_to_global(MaterialData * md, arma::mat & S)
{
  arma::mat33 R, aux;
  R.col(0) = md->fiber();
  R.col(1) = md->sheet();
  R.col(2) = md->normal();
  S = R * S * R.t();
}

void HyperelasticMaterial::map_to_local(MaterialData * md, arma::mat & S)
{
  arma::mat33 R, aux;
  R.col(0) = md->fiber();
  R.col(1) = md->sheet();
  R.col(2) = md->normal();
  S = R.t() * S * R;
}

void HyperelasticMaterial::calc_fd_active_elastensor(MaterialData * md, Tensor4 & A)
{
  const double eps  = 1.0e-4;
  const double eps2 = eps*eps;
  double psi0, psip, psim, psipp, psimp, psipm, psimm;
  arma::mat33 E, Epdij, Emdij, dij, dkl;

  E = md->lagrangian_strain();

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<=i; j++) // old j<3
    {
      delta_tensor(i, j, eps, dij);
      for(int k=0; k<3; k++)
      {
        for(int l=0; l<=k; l++) //old l<3
        {
          delta_tensor(k, l, eps, dkl);
          Epdij = E + dij;
          Emdij = E - dij;
          if(i==k && j==l)
          {
            psip = active_strain_energy(md, Epdij);
            psim = active_strain_energy(md, Emdij);
            psi0 = active_strain_energy(md, E);
            A(i,j,k,l) = (psip - 2.*psi0 + psim)/eps2;
          }
          else
          {
            psipp = active_strain_energy(md, Epdij + dkl);
            psimp = active_strain_energy(md, Emdij + dkl);
            psipm = active_strain_energy(md, Epdij - dkl);
            psimm = active_strain_energy(md, Emdij - dkl);
            A(i,j,k,l) = (psipp - (psimp + psipm) + psimm)/(4.*eps2);
          }
        }
      }
    }
  }

  // TODO: check if UNROLLING this will improve performance

  // copy minor symmetries (A_ijkl = A_jikl = A_ijlk)
  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          A(j,i,k,l) = A(i,j,k,l);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<=k; l++)
          A(i,j,l,k) = A(i,j,k,l);
}

void HyperelasticMaterial::calc_fd_elastensor(MaterialData * md, Tensor4 & A)
{
  const double eps  = 1.0e-4;
  const double eps2 = eps*eps;
  double psi0, psip, psim, psipp, psimp, psipm, psimm;
  arma::mat33 E, Epdij, Emdij, dij, dkl;

  E = md->lagrangian_strain();

  for(int i=0; i<3; i++)
  {    
    for(int j=0; j<=i; j++) // old j<3
    {
      delta_tensor(i, j, eps, dij);
      for(int k=0; k<3; k++)
      {	      
	      for(int l=0; l<=k; l++) //old l<3
        {
          delta_tensor(k, l, eps, dkl);	  
          Epdij = E + dij;
          Emdij = E - dij;
          if(i==k && j==l)
          {
            psip = strain_energy(md, Epdij);
            psim = strain_energy(md, Emdij);
            psi0 = strain_energy(md, E);
            A(i,j,k,l) = (psip - 2.*psi0 + psim)/eps2;
          }
          else
          {
            psipp = strain_energy(md, Epdij + dkl);
            psimp = strain_energy(md, Emdij + dkl);
            psipm = strain_energy(md, Epdij - dkl);
            psimm = strain_energy(md, Emdij - dkl);
            A(i,j,k,l) = (psipp - (psimp + psipm) + psimm)/(4.*eps2);
          }
	      }
      }
    }
  }

  // TODO: check if UNROLLING this will improve performance
  
  // copy minor symmetries (A_ijkl = A_jikl = A_ijlk)  
  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
      for(int k=0; k<3; k++)
	      for(int l=0; l<3; l++) 
          A(j,i,k,l) = A(i,j,k,l);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<=k; l++) 
          A(i,j,l,k) = A(i,j,k,l);

  active_stress_elastensor(1,md,A);

}

void HyperelasticMaterial::calc_fd_elastensor(MaterialData * md, arma::mat & D)
{
  // first compute using finite difference
  Tensor4 A;
  calc_fd_elastensor(md, A);
  //active_stress_elastensor(1,md,A);

  // then transform 4th order tensor to matrix (Voigt) form
  A.get_matrix(D);
}

void HyperelasticMaterial::push_forward(const arma::mat33 & F,
                                        const arma::mat33 & S,
                                        arma::mat33 & sigma)
{
  double J = arma::det(F);
  sigma = (1./J) * F * S * F.t(); 
}

void HyperelasticMaterial::push_forward(const arma::mat33 & F, Tensor4 & Am,  
                                        Tensor4 & As)
{
  const double J = arma::det(F);
  const double oneJ = 1.0/J;

  As.zero();

  for(int a=0; a<3; a++)
    for(int b=0; b<3; b++)
      for(int c=0; c<3; c++)
	      for(int d=0; d<3; d++)
          for(int A=0; A<3; A++)
            for(int B=0; B<3; B++)
              for(int C=0; C<3; C++)
                for(int D=0; D<3; D++)
		              As(a,b,c,d) += oneJ * F(a,A)*F(b,B)*F(c,C)*F(d,D)*Am(A,B,C,D);
}


void HyperelasticMaterial::map_elas_to_global(MaterialData * md, Tensor4 & Am,
                                        Tensor4 & As)
{
    arma::mat33 R;

    R.col(0) = md->fiber();
    R.col(1) = md->sheet();
    R.col(2) = md->normal();

    As.zero();

    for(int a=0; a<3; a++)
        for(int b=0; b<3; b++)
            for(int c=0; c<3; c++)
                for(int d=0; d<3; d++)
                    for(int A=0; A<3; A++)
                        for(int B=0; B<3; B++)
                            for(int C=0; C<3; C++)
                                for(int D=0; D<3; D++)
                                        As(A,B,C,D) += R(A,a)*R(B,b)*R(C,c)*R(D,d)*Am(a,b,c,d);

}

void HyperelasticMaterial::map_elas_to_local(MaterialData * md, Tensor4 & Am,
                                              Tensor4 & As)
{
  arma::mat33 R;

  R.col(0) = md->fiber();
  R.col(1) = md->sheet();
  R.col(2) = md->normal();

  As.zero();

  for(int a=0; a<3; a++)
    for(int b=0; b<3; b++)
      for(int c=0; c<3; c++)
        for(int d=0; d<3; d++)
          for(int A=0; A<3; A++)
            for(int B=0; B<3; B++)
              for(int C=0; C<3; C++)
                for(int D=0; D<3; D++)
                  As(a,b,c,d) += R(a,A)*R(b,B)*R(c,C)*R(d,D)*Am(A,B,C,D);

}

HyperelasticMaterial *
HyperelasticMaterial::create(std::string matname,
                             ElasticityType & elastype,
                             const std::vector<double> & matprop)
{
  if(matprop.size() == 0)
    throw runtime_error("HyperelasticMaterial without parameters!");

  HyperelasticMaterial *mat = NULL;

  if (matname == "CompressibleNeoHookean")
    mat = new CompressibleNeoHookean(matprop);
  else if (matname == "SaintVenantKirchoff")
    mat = new SaintVenantKirchoff(matprop);
  else if (matname == "NeoHookean")
    mat = new NeoHookean(matprop);
  else if (matname == "Exponential")
    mat = new Exponential(matprop);
  else if (matname == "Guccione")
    mat = new Guccione(matprop);
  else if (matname == "HolzapfelOgden")
    mat = new HolzapfelOgden(matprop);
  else
    throw runtime_error("unknown HyperelasticMaterial");

  // setup number of dimensions
  if (elastype == PLANE_STRAIN)
    mat->set_ndim(2);
  else if (elastype == THREE_DIM)
    mat->set_ndim(3);

  return mat;
};

HyperelasticMaterial *
HyperelasticMaterial::create(std::string matname,
                             ElasticityType & elastype,
                             const std::vector<std::vector<double>> & matprop, int num_materials, std::vector<int> & map_AHAmat)
{
  if(matprop.size() == 0)
    throw runtime_error("HyperelasticMaterial without parameters!");

  HyperelasticMaterial *mat = NULL;

  if (matname == "Guccione")
    mat = new Guccione(matprop, num_materials, map_AHAmat);
  else
    if(matname == "HolzapfelOgden" )
      mat = new HolzapfelOgden(matprop, num_materials, map_AHAmat);
  else
    throw runtime_error("unknown HyperelasticMaterial: Mutiple materials not suported");

  // setup number of dimensions
  if (elastype == PLANE_STRAIN)
    mat->set_ndim(2);
  else if (elastype == THREE_DIM)
    mat->set_ndim(3);

  return mat;
};

double HyperelasticMaterial::penalty_term(MaterialData * md,
                                          const arma::mat & E) const
{
  arma::mat33 I = arma::eye<arma::mat>(3,3);
  arma::mat33 C = 2*E + I;
  double J = sqrt(arma::det(C));  
  double K = parameters[2];
  double t = (K/2.0)*((J-1.0)*(J-1.0));
  return t;
}

void HyperelasticMaterial::calc_fd_stress_penalty(MaterialData * md, arma::mat & S)					  
{
  const double eps = 1.0e-8;
  double psi_m, psi_p;
  arma::mat33 E, dij;

  E = md->lagrangian_strain();

  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
    {
      delta_tensor(i, j, eps, dij);
      psi_p = penalty_term(md, E + dij);
      psi_m = penalty_term(md, E - dij);
      S(i,j) = (psi_p - psi_m)/(2.*eps);
    }

  // copy symmetric entries
  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
      S(j,i) = S(i,j);
}

void HyperelasticMaterial::calc_fd_elastensor_penalty(MaterialData * md, arma::mat & D)
{
  Tensor4 A;
  calc_fd_elastensor_penalty(md, A);
  A.get_matrix(D);
}
  
void HyperelasticMaterial::calc_fd_elastensor_penalty(MaterialData * md, Tensor4 & A)
{
  const double eps  = 1.0e-4;
  const double eps2 = eps*eps;
  double psi0, psip, psim, psipp, psimp, psipm, psimm;
  arma::mat33 E, Epdij, Emdij, dij, dkl;

  E = md->lagrangian_strain();

  for(int i=0; i<3; i++)
  {    
    for(int j=0; j<=i; j++) // old j<3
    {
      delta_tensor(i, j, eps, dij);
      for(int k=0; k<3; k++)
      {	      
	      for(int l=0; l<=k; l++) //old l<3
        {
          delta_tensor(k, l, eps, dkl);	  
          Epdij = E + dij;
          Emdij = E - dij;
          if(i==k && j==l)
          {
            psip = penalty_term(md, Epdij);
            psim = penalty_term(md, Emdij);
            psi0 = penalty_term(md, E);
            A(i,j,k,l) = (psip - 2.*psi0 + psim)/eps2;
          }
          else
          {
            psipp = penalty_term(md, Epdij + dkl);
            psimp = penalty_term(md, Emdij + dkl);
            psipm = penalty_term(md, Epdij - dkl);
            psimm = penalty_term(md, Emdij - dkl);
            A(i,j,k,l) = (psipp - (psimp + psipm) + psimm)/(4.*eps2);
          }
	      }
      }
    }
  }
   
  // copy minor symmetries (A_ijkl = A_jikl = A_ijlk)  
  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
      for(int k=0; k<3; k++)
	      for(int l=0; l<3; l++) 
          A(j,i,k,l) = A(i,j,k,l);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<=k; l++) 
          A(i,j,l,k) = A(i,j,k,l);
}


// ----------------------------------------------------------------------------

/*
void HyperelasticMaterial::calc_fd_stress(arma::mat & E, arma::mat & S,
					  MaterialData * md)
{
  double psi_m, psi_p, eold;
  const double eps = 1.0e-8;
  //arma::mat33 dij;

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      //delta_tensor(i, j, eps, dij);
      eold = E(i,j);
      
      // E + dij
      if(i==j) E(i,j) += eps;
      else
      {
	E(i,j) += eps/2.;
	E(j,i) += eps/2.;
      }
      psi_p = strain_energy(md, E);

      E(i,j) = eold;
      if(i!=j) E(j,i) = eold;
      //psi_p = strain_energy(md, E + dij);

      // E - dij
     if(i==j) E(i,j) -= eps;
      else
      {
	E(i,j) -= eps/2.;
	E(j,i) -= eps/2.;
      }
      psi_m = strain_energy(md, E);

      E(i,j) = eold;
      if(i!=j) E(j,i) = eold;
      //psi_m = strain_energy(md, E - dij);

      S(i,j) = (psi_p - psi_m)/(2.0*eps);
    }
}

void HyperelasticMaterial::calc_fd_elastensor(arma::mat & E, arma::mat & D,
						MaterialData * md)
{
  double psi0, psip, psim, psipp, psimp, psipm, psimm, eold, eold2;
  const double eps  = 1.0e-4;
  const double eps2 = eps*eps;
  arma::mat33 dij, dkl, Epdij, Emdij;

  Tensor4 elastensor;
  
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      //delta_tensor(i, j, eps, dij);
      for(int k=0; k<3; k++)
      {
	for(int l=0; l<3; l++)
	{
	  //delta_tensor(k, l, eps, dkl);
	  
	  //Epdij = E + dij;
	  //Emdij = E - dij;

	  if(i==k && j==l)
	  {
	    // E + Dij
	    eold = E(i,j);
	    if(i==j) E(i,j) += eps;
	    else
	    {
	      E(i,j) += eps/2.;
	      E(j,i) += eps/2.;
	    }
	    psip = strain_energy(md, E);
	    E(i,j) = eold;
	    if(i!=j) E(j,i) = eold;
	    //psip = strain_energy(md, Epdij);

	    // E - Dij
	    eold = E(i,j);
	    if(i==j) E(i,j) -= eps;
	    else
	    {
	      E(i,j) -= eps/2.;
	      E(j,i) -= eps/2.;
	    }
	    psim = strain_energy(md, E);
	    E(i,j) = eold;
	    if(i!=j) E(j,i) = eold;	    
	    //psim = strain_energy(md, Emdij);

	    psi0 = strain_energy(md, E);

	    elastensor(i,j,k,l) = (psip-2.0*psi0+psim)/eps2;
	  }
	  else
	  {
	    // E + Dij + Dkl -----------------------------------
	    eold  = E(i,j);
	    eold2 = E(k,l);
	    if(i==j) E(i,j) += eps;
	    else
	    {
	      E(i,j) += eps/2.;
	      E(j,i) += eps/2.;
	    }
	    if(k==l) E(k,l) += eps;
	    else
	    {
	      E(k,l) += eps/2.;
	      E(l,k) += eps/2.;
	    }	    
	    psipp = strain_energy(md, E);
	    E(i,j) = eold;
	    if(i!=j) E(j,i) = eold;
	    E(k,l) = eold2;
	    if(k!=l) E(l,k) = eold2;
	    //psipp = strain_energy(md, E+dij + dkl);

	    // E - Dij + Dkl -----------------------------------
	    eold = E(i,j);
	    eold2 = E(k,l);
	    if(i==j) E(i,j) -= eps;
	    else
	    {
	      E(i,j) -= eps/2.;
	      E(j,i) -= eps/2.;
	    }
	    if(k==l) E(k,l) += eps;
	    else
	    {
	      E(k,l) += eps/2.;
	      E(l,k) += eps/2.;
	    }	    
	    psimp = strain_energy(md, E);
	    E(i,j) = eold;
	    if(i!=j) E(j,i) = eold;
	    E(k,l) = eold2;
	    if(k!=l) E(l,k) = eold2;
	    //psimp = strain_energy(md, E-dij + dkl);

	    // E + Dij - Dkl ----------------------------------------
	    eold = E(i,j);
	    eold2 = E(k,l);
	    if(i==j) E(i,j) += eps;
	    else
	    {
	      E(i,j) += eps/2.;
	      E(j,i) += eps/2.;
	    }
	    if(k==l) E(k,l) -= eps;
	    else
	    {
	      E(k,l) -= eps/2.;
	      E(l,k) -= eps/2.;
	    }	    
	    psipm = strain_energy(md, E);
	    E(i,j) = eold;
	    if(i!=j) E(j,i) = eold;
	    E(k,l) = eold2;
	    if(k!=l) E(l,k) = eold2;
	    //psipm = strain_energy(md, E + dij - dkl);
	    
	    // E - Dij - Dkl
	    eold = E(i,j);
	    eold2 = E(k,l);
	    if(i==j) E(i,j) -= eps;
	    else
	    {
	      E(i,j) -= eps/2.;
	      E(j,i) -= eps/2.;
	    }
	    if(k==l) E(k,l) -= eps;
	    else
	    {
	      E(k,l) -= eps/2.;
	      E(l,k) -= eps/2.;
	    }	    
	    psimm = strain_energy(md, E);
	    E(i,j) = eold;
	    if(i!=j) E(j,i) = eold;
	    E(k,l) = eold2;
	    if(k!=l) E(l,k) = eold2;
	    //psimm = strain_energy(md, E-dij - dkl);
	    elastensor(i,j,k,l) = (psipp-(psimp+psipm)+psimm)/(4.0*eps2);
	  }
	}
      }
    }
  }

  tensor4_to_matrix(elastensor, D);   
}
*/
