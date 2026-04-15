#include "holzapfel_ogden.hpp"

//
//#define DEBUG
//

double HolzapfelOgden::strain_energy(MaterialData * md, const arma::mat & E) const
{
#ifndef DEBUG
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();
#endif

#ifdef DEBUG
  arma::vec3 f0;
  arma::vec3 s0;

  f0 = md->fiber();
  s0(0) = -f0(1); s0(1) = f0(0); s0(2) = 0;
  //s0 = arma::norm(s0,2);

#endif

  double a = v_a[map_mat[md->get_marker()]] ;
  double af = v_af[map_mat[md->get_marker()]] ;
  double as = v_as[map_mat[md->get_marker()]] ;
  double afs = v_afs[map_mat[md->get_marker()]] ;
  double b = v_b[map_mat[md->get_marker()]] ;
  double bf = v_bf[map_mat[md->get_marker()]] ;
  double bs = v_bs[map_mat[md->get_marker()]] ;
  double bfs = v_bfs[map_mat[md->get_marker()]] ;

  // right Cauchy-Green deformation tensor
  arma::mat C = 2*E + I;
  double J = sqrt(arma::det(C));


  // isochoric right Cauchy-Green deformation tensor and its trace
  arma::mat Cbar = pow(J,-(2.0/3.0)) * C;
  double Icbar = arma::trace(Cbar);

  // invariants
  const arma::vec3 cs = Cbar*s0;
  double I4f  = dot(f0, Cbar*f0);
  double I4s  = dot(s0, cs);
  double I8fs = dot(f0, cs);

  const double K = parameters[9];

  // strain energy function calculation
  double tiso = (a/(2.0*b)) * (exp(b*(Icbar-3.0)) - 1.0);
  double tfib = (I4f > 1.0) ? (af/(2*bf))*(exp(bf*(I4f-1.)*(I4f-1.))-1.) : 0.0;
  //double tfib = (af/(2*bf))*(exp(bf*(I4f-1.)*(I4f-1.))-1.);
  double tshe = (I4s > 1.0) ? (as/(2*bs))*(exp(bs*(I4s-1.)*(I4s-1.))-1.) : 0.0;
  double tfs  = (afs/(2*bfs)) * (exp(bfs*(I8fs*I8fs))-1.0);
  double tvol = (K/2.0)*((J-1.0)*(J-1.0));

  return tiso + tfib + tvol;// + active_stress*I4f;
  //return tiso + tfib + tshe + tfs + tvol;
  //return tiso + tvol;
}

double HolzapfelOgden::active_strain_energy(int iel, MaterialData * md, const arma::mat & E) const
{
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();

  // right Cauchy-Green deformation tensor
  arma::mat C = 2*E + I;
  double J = sqrt(arma::det(C));

  arma::mat Cbar = pow(J,-(2.0/3.0)) * C;
  double I4f  = dot(f0, Cbar*f0);

  return I4f*active_stress;
}

void HolzapfelOgden::deviatoric_stress(MaterialData * md, arma::mat & sigma) const
{
#ifndef DEBUG
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();
#endif

#ifdef DEBUG
  arma::vec3 f0;
  arma::vec3 s0;
  f0.zeros();
  s0.zeros();
  f0(0) = 1.;
  s0(1) = 1.;
#endif

  double a = v_a[map_mat[md->get_marker()]] ;
  double af = v_af[map_mat[md->get_marker()]] ;
  double as = v_as[map_mat[md->get_marker()]] ;
  double afs = v_afs[map_mat[md->get_marker()]] ;
  double b = v_b[map_mat[md->get_marker()]] ;
  double bf = v_bf[map_mat[md->get_marker()]] ;
  double bs = v_bs[map_mat[md->get_marker()]] ;
  double bfs = v_bfs[map_mat[md->get_marker()]] ;

  // Kinematics
  double J = md->get_J();
  arma::mat33 F = md->get_F();
  arma::mat33 Fbar = md->isochoric_def_grad();
  arma::mat33 Bbar = md->isochoric_lcg();
  arma::mat33 Cbar = md->isochoric_rcg();

  // Modified Invariants
  arma::vec3 Cbar_s0 = Cbar * s0;
  double I1   = arma::trace(Cbar);
  double I4f  = dot(f0, Cbar * f0);
  double I4s  = dot(s0, Cbar_s0);
  double I8fs = dot(f0, Cbar_s0);

  // Derivatives of the W_iso with respect to the modified invariants
  double psi1 = (I1 != 3.0) ? (a / 2.) * exp(b * (I1 - 3.)) : 0;
  double psi4f = af * (I4f - 1.) * exp(bf * (I4f - 1.) * (I4f - 1.));
  double psi4s = as * (I4s - 1.) * exp(bs * (I4s - 1.) * (I4s - 1.));
  double psi8fs = afs * I8fs * exp(bfs * I8fs * I8fs);

  // Tensor products
  arma::vec3 f_bar = Fbar * f0;
  arma::vec3 s_bar = Fbar * s0;
  arma::mat33 fs = f_bar * s_bar.t();
  arma::mat33 sf = s_bar * f_bar.t();
  arma::mat33 fssf = fs + sf;

  // Compute stress
  //sigma = (2.0/J) * (psi1*dev(Bbar) + psi4f*dev(ff) + psi4s*dev(ss) + 0.5*psi8fs*dev(fssf));

  sigma = (2. / J) * (psi1 * sp_dev_operator(Bbar));

  if (I4f > 1)
  {
    arma::mat33 ff = f_bar * f_bar.t();
    sigma += (2. / J) * (psi4f * sp_dev_operator(ff));
  }

  if (I4s > 1)
  {
    arma::mat33 ss = s_bar * s_bar.t();
    sigma += (2. / J) * (psi4s * sp_dev_operator(ss));
  }

  sigma += (2.0 / J) * (1. / 2.) * psi8fs * sp_dev_operator(fssf);
}

void HolzapfelOgden::deviatoric_elastensor(MaterialData * md, Tensor4 & A) const
{
  static const arma::mat33 I = arma::eye(3,3);
  static const Tensor4 II = unit_tensor();

#ifndef DEBUG
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();
#endif

#ifdef DEBUG
  arma::vec3 f0;
  arma::vec3 s0;
  f0.zeros();
  s0.zeros();
  f0(0)=1.;
  s0(1)=1.;
#endif

  double a = v_a[map_mat[md->get_marker()]] ;
  double af = v_af[map_mat[md->get_marker()]] ;
  double as = v_as[map_mat[md->get_marker()]] ;
  double afs = v_afs[map_mat[md->get_marker()]] ;
  double b = v_b[map_mat[md->get_marker()]] ;
  double bf = v_bf[map_mat[md->get_marker()]] ;
  double bs = v_bs[map_mat[md->get_marker()]] ;
  double bfs = v_bfs[map_mat[md->get_marker()]] ;

  // kinematics
  //arma::mat33 Ctens = md->right_cauchy_green();

  arma::mat33 Fbar  = md->isochoric_def_grad();
  arma::mat33 Bbar  = md->isochoric_lcg() ;
  arma::mat33 Cbar  = md->isochoric_rcg();

  // modified Invariants
  arma::vec3 Cbar_s0 = Cbar * s0;

  double I1   = arma::trace(Cbar);
  double I4f  = dot(f0, Cbar*f0);
  double I4s  = dot(s0, Cbar_s0);
  double I8fs = dot(f0, Cbar_s0);

  // derivatives of the W_iso with respect to the modified invariants
  double I8fs_2    = I8fs * I8fs;
  double I4fm1_2   = (I4f-1.0) * (I4f-1.0);
  double I4sm1_2   = (I4s-1.0) * (I4s-1.0);

  double psi1      = (I1 != 3.0) ? (a/2.) * exp(b*(I1-3.)) : 0;
  double psi4f     = af * (I4f-1.) * exp(bf * I4fm1_2);
  double psi4s     = as * (I4s-1.) * exp(bs * I4sm1_2);
  double psi8fs    = afs * I8fs * exp(bfs*I8fs_2);

  double psi11     = (I1 != 3.0) ? a * b * exp(b * (I1 - 3.)) : 0;
  double psi4f4f   = af * (1. + 2.* bf * I4fm1_2) * exp(bf * I4fm1_2);
  double psi4s4s   = as * (1. + 2.* bs * I4sm1_2) * exp(bs * I4sm1_2);
  double psi8fs8fs = afs * (1. + 2.* bfs * I8fs_2) * exp(bfs * I8fs_2);

  // tensor products
  arma::vec3 f_bar = Fbar * f0;
  arma::vec3 s_bar = Fbar * s0;
  arma::mat33 ff   = f_bar * f_bar.t();
  arma::mat33 ss   = s_bar * s_bar.t();
  arma::mat33 fssf = f_bar * s_bar.t() + s_bar * f_bar.t();

  // prepare some tensors
  arma::mat33 devb    = sp_dev_operator(Bbar);
  arma::mat33 devff   = sp_dev_operator(ff);
  arma::mat33 devss   = sp_dev_operator(ss);
  arma::mat33 devfssf = sp_dev_operator(fssf);

  // compute c_iso_I1 + c_iso_I4f + c_iso_I4s + c_iso_I8fs
  A.zero();

  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
        for(int l=0; l<ndim; l++)
        {
          double P = II(i,j,k,l) - (1./3.)*(I(i,j)*I(k,l));

          // c_iso_I1
          A(i,j,k,l) += 4.*psi11*(devb(i,j)*devb(k,l)) - (4./3.)*psi1*( (I(i,j)*devb(k,l)) + (devb(i,j)*I(k,l)) - I1*P );

          // c_iso_I4f
          if(I4f > 1)
            A(i,j,k,l) += 4.*psi4f4f*(devff(i,j)*devff(k,l)) - (4./3.)*psi4f*( (I(i,j)*devff(k,l)) + (devff(i,j)*I(k,l)) - I4f*P );

          // c_iso_I4s
          if(I4s > 1)
            A(i,j,k,l) += 4.*psi4s4s*(devss(i,j)*devss(k,l)) - (4./3.)*psi4s*( (I(i,j)*devss(k,l)) + (devss(i,j)*I(k,l)) - I4s*P );

          // c_iso_I8fs
          A(i,j,k,l) += psi8fs8fs*(devfssf(i,j)*devfssf(k,l)) - (2./3.)*psi8fs*( (I(i,j)*devfssf(k,l)) + (devfssf(i,j)*I(k,l)) - 2.*I8fs*P );
        }

/*
  // Compute c_iso_I1
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
	    for(int l=0; l<ndim; l++){
          double P = II(i,j,k,l) - (1./3.)*(I(i,j)*I(k,l));
          A(i, j, k, l) += 4. * psi11 * (devb(i, j) * devb(k, l)) - (4. / 3.) * psi1 * ((I(i, j) * devb(k, l)) + (devb(i, j) * I(k, l)) - I1 * P);
        }

  // Compute c_iso_I4f
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
	    for(int l=0; l<ndim; l++)
        {
          double P = II(i,j,k,l) - (1./3.)*(I(i,j)*I(k,l));
          A(i, j, k, l) += 4. * psi4f4f * (devff(i, j) * devff(k, l)) - (4. / 3.) * psi4f * ((I(i, j) * devff(k, l)) + (devff(i, j) * I(k, l)) - I4f * P);
        }

  // Compute c_iso_I4s
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
	    for(int l=0; l<ndim; l++)
        {
          double P = II(i,j,k,l) - (1./3.)*(I(i,j)*I(k,l));
          A(i, j, k, l) += 4. * psi4s4s * (devss(i, j) * devss(k, l)) - (4. / 3.) * psi4s * ((I(i, j) * devss(k, l)) + (devss(i, j) * I(k, l)) - I4s * P);
        }

  // Compute c_iso_I8fs
  for(int i=0; i<ndim; i++)
    for(int j=0; j<ndim; j++)
      for(int k=0; k<ndim; k++)
	    for(int l=0; l<ndim; l++)
        {
          double P = II(i,j,k,l) - (1./3.)*(I(i,j)*I(k,l));
          A(i, j, k, l) += psi8fs8fs * (devfssf(i, j) * devfssf(k, l)) - (2. / 3.) * psi8fs * ((I(i, j) * devfssf(k, l)) + (devfssf(i, j) * I(k, l)) - 2. * I8fs * P);
        }
*/
}


void HolzapfelOgden::piola2_stress(MaterialData * md, arma::mat & pk2) const {
#ifndef DEBUG
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();
#endif

#ifdef DEBUG
  arma::vec3 f0;
  arma::vec3 s0;
  f0.zeros();
  s0.zeros();
  f0(0) = 1.;
  s0(1) = 1.;
#endif

  double a = v_a[map_mat[md->get_marker()]] ;
  double af = v_af[map_mat[md->get_marker()]] ;
  double as = v_as[map_mat[md->get_marker()]] ;
  double afs = v_afs[map_mat[md->get_marker()]] ;
  double b = v_b[map_mat[md->get_marker()]] ;
  double bf = v_bf[map_mat[md->get_marker()]] ;
  double bs = v_bs[map_mat[md->get_marker()]] ;
  double bfs = v_bfs[map_mat[md->get_marker()]] ;

  static const arma::mat33 I = arma::eye(3,3);
  // Kinematics
  double J = md->get_J();
  //arma::mat33 Ctens = md->right_cauchy_green();


  arma::mat33 Fbar  = md->isochoric_def_grad();
  arma::mat33 Bbar  = md->isochoric_lcg();
  arma::mat33 Cbar  = md->isochoric_rcg();
  arma::mat33 F = md->get_F();
  arma::mat33 C    = md->right_cauchy_green();

  // Modified Invariants
  arma::vec3 Cbar_s0 = C * s0;
  double I1 = arma::trace(C);
  double I4f = dot(f0, C * f0);
  double I4s = dot(s0, Cbar_s0);
  double I8fs = dot(f0, Cbar_s0);

  // Derivatives of the W_iso with respect to the modified invariants
  double psi1   = (I1 != 3.0) ? (a / 2.) * exp(b * (I1 - 3.) ) : 0;
  double psi4f  = (I4f > 1.0) ? af * (I4f - 1.) * exp(bf * (I4f - 1.) * (I4f - 1.)) : 0;
  double psi4s  = (I4s > 1.0) ? as * (I4s - 1.) * exp(bs * (I4s - 1.) * (I4s - 1.)) : 0;
  double psi8fs = afs * I8fs * exp(bfs * I8fs * I8fs);


  // Tensor products
  arma::mat33 f0f0     = f0 * f0.t();
  arma::mat33 s0s0     = s0 * s0.t();
  arma::mat33 f0s0s0f0 = ( f0 * s0.t() ) + (s0 * f0.t());

  arma::mat33 DevI    = mt_dev_operator(I,C);
  arma::mat33 Devf0f0 = mt_dev_operator(f0f0,C);
  arma::mat33 Devs0s0 = mt_dev_operator(s0s0,C);
  arma::mat33 Devf0s0 = mt_dev_operator(f0s0s0f0,C);

  double jm23 = pow(J,-2.0/3.0);

  double lamb = 1.05;
  double lamb2 = lamb*lamb;
  double PSI1 = (a/2.0)*exp(b*(lamb2 + (2.0/lamb) - 3.0));
  double PSI4f = af * (lamb2 - 1.0) * exp(bf*(lamb2 - 1)*(lamb2 - 1));

  double sigma11 = 2*PSI1*(lamb2 - (1.0/lamb)) + 2.0*PSI4f*lamb2;
  double ph = -2.0*PSI1*(1.0/lamb);
  double sigma22 = ph + 2.0*PSI1*(1.0/lamb);
  double sigma33 = ph + 2.0*PSI1*(1.0/lamb);

  pk2.zeros();
  // volumetric
  // isochoric
  pk2 += 2.0 * jm23 * ( psi1*DevI + psi4f*Devf0f0 + psi4s*Devs0s0 + psi8fs*0.5*Devf0s0 );
  pk2 += ph * J * C.i();


}
/*
void HolzapfelOgden::deviatoric_elastensor(MaterialData * md, Tensor4 & A) const
{
  static const arma::mat33 I = arma::eye(3,3);
  static const Tensor4 II = unit_tensor();

#ifndef DEBUG
  const arma::vec3 & f0 = md->fiber();
  const arma::vec3 & s0 = md->sheet();
#endif

#ifdef DEBUG
  arma::vec3 f0;
  arma::vec3 s0;
  f0.zeros();
  s0.zeros();
  f0(0)=1.;
  s0(1)=1.;
#endif

  A.zero();

}

*/
