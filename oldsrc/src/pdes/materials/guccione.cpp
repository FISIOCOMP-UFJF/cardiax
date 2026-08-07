#include "guccione.hpp"



double Guccione::strain_energy(MaterialData * md, const arma::mat & E) const
{
  arma::mat33 El, R;


  double Cg = v_Cg[map_AHA[md->get_marker()]] ;
  double bf = v_bf[map_AHA[md->get_marker()]] ;
  double bt = v_bt[map_AHA[md->get_marker()]] ;
  double bfs = v_bfs[map_AHA[md->get_marker()]] ;


  R.col(0) = md->fiber();
  R.col(1) = md->sheet();
  R.col(2) = md->normal();

  El = R.t() * E * R;

  // Right Cauchy-Green deformation tensor and its trace
  arma::mat C = 2*El + I;

  // Jacobian
  double J = sqrt(arma::det(C));

  arma::mat Cbar = pow(J,-(2.0/3.0)) * C;
  arma::mat Ebar = 0.5*(Cbar - I);

  double Q = bf*El(0,0)*El(0,0) \
           + bt*( El(1,1)*El(1,1) + El(2,2)*El(2,2) + El(1,2)*El(1,2) + El(2,1)*El(2,1) ) \
           + bfs*( El(0,1)*El(0,1) + El(1,0)*El(1,0) + El(0,2)*El(0,2) + El(2,0)*El(2,0) );

  const double K = parameters[5];
  double term1 = (Cg/2.)*( exp(Q) - 1. );
  double term2 = (K/2.0)*((J-1)*(J-1));
  return term1 + term2;
}


double Guccione::active_strain_energy(MaterialData * md,
                                                    const arma::mat &) const
{
  return 0;
}

void Guccione::cauchy_stress(MaterialData * md, arma::mat & sigma)  const
{
  throw runtime_error("Guccione: not implemented");
}

void Guccione::piola2_stress(MaterialData * md, arma::mat & S) const
{

  double Cg = v_Cg[map_AHA[md->get_marker()]] ;
  double bf = v_bf[map_AHA[md->get_marker()]] ;
  double bt = v_bt[map_AHA[md->get_marker()]] ;
  double bfs = v_bfs[map_AHA[md->get_marker()]] ;

  arma::mat33 E = md->lagrangian_strain();
  arma::mat33 C = 2 * E + I;
  double J = sqrt(arma::det(C));
  arma::mat Cbar = pow(J, -(2.0 / 3.0)) * C;
  E = 0.5 * (Cbar - I);

  arma::mat33 El, R;

  R.col(0) = md->fiber();
  R.col(1) = md->sheet();
  R.col(2) = md->normal();



  El = R.t() * E * R;


  const double K = parameters[5];
  double p = K * (J - 1);

  // some constants to improve
  const double ce00p2 = pow(El(0, 0), 2);
  const double ce01p2 = pow(El(0, 1), 2);
  const double ce02p2 = pow(El(0, 2), 2);
  const double ce10p2 = pow(El(1, 0), 2);
  const double ce11p2 = pow(El(1, 1), 2);
  const double ce12p2 = pow(El(1, 2), 2);
  const double ce20p2 = pow(El(2, 0), 2);
  const double ce21p2 = pow(El(2, 1), 2);
  const double ce22p2 = pow(El(2, 2), 2);

  const double expterm = exp(bf * ce00p2
                           + bfs* (ce01p2 + ce02p2 + ce10p2 + ce20p2)
                           + bt * (ce11p2 + ce12p2 + ce21p2 + ce22p2));

  S(0, 0) = Cg * El(0, 0) * bf * expterm;
  S(0, 1) = Cg * El(0, 1) * bfs * expterm;
  S(0, 2) = Cg * El(0, 2) * bfs * expterm;
  S(1, 0) = Cg * El(1, 0) * bfs * expterm;
  S(1, 1) = Cg * El(1, 1) * bt * expterm;
  S(1, 2) = Cg * El(1, 2) * bt * expterm;
  S(2, 0) = Cg * El(2, 0) * bfs * expterm;
  S(2, 1) = Cg * El(2, 1) * bt * expterm;
  S(2, 2) = Cg * El(2, 2) * bt * expterm;

  // Old, working but not efficient
  //  S(0,0) = Cg*El(0,0)*bf*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(0,0);
  //  S(0,1) = Cg*El(0,1)*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(0,1);
  //  S(0,2) = Cg*El(0,2)*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(0,2);
  //  S(1,0) = Cg*El(1,0)*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(1,0);
  //  S(1,1) = Cg*El(1,1)*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(1,1);
  //  S(1,2) = Cg*El(1,2)*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(1,2);
  //  S(2,0) = Cg*El(2,0)*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(2,0);
  //  S(2,1) = Cg*El(2,1)*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(2,1);
  //  S(2,2) = Cg*El(2,2)*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));// + p*J*invC(2,2);
}

void Guccione::sp_elastensor(MaterialData * md, arma::mat & D) const
{
  throw runtime_error("Guccione: not implemented");
}

void Guccione::mt_elastensor(MaterialData * md, arma::mat & D) const
{
  throw runtime_error("Guccione: not implemented");
}

void Guccione::deviatoric_stress(MaterialData * md, arma::mat & stress) const
{
  throw runtime_error("Guccione: not implemented");
}

void Guccione::deviatoric_elastensor(MaterialData * md, Tensor4 & A) const
{

  double Cg = v_Cg[map_AHA[md->get_marker()]] ;
  double bf = v_bf[map_AHA[md->get_marker()]] ;
  double bt = v_bt[map_AHA[md->get_marker()]] ;
  double bfs = v_bfs[map_AHA[md->get_marker()]] ;



  //throw runtime_error("Guccione: not implemented");

  arma::mat33 E = md->lagrangian_strain();
  arma::mat33 C = 2*E + I;
  double J = sqrt(arma::det(C));
  arma::mat Cbar = pow(J,-(2.0/3.0)) * C;
  E = 0.5*(Cbar - I);
  arma::mat33 El, R;

  R.col(0) = md->fiber();
  R.col(1) = md->sheet();
  R.col(2) = md->normal();

  El = R.t() * E * R;
  arma::mat33 invC = arma::inv(C);
  double K = parameters[5], p = K*(J-1), Ck, Cp;

  const double expterm = exp(bf*pow(El(0,0), 2)
                           + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2))
                           + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));

  A(0,0,0,0) = 2*Cg*pow(El(0,0), 2)*pow(bf, 2)*expterm + Cg*bf*expterm;
  A(0,0,0,1) = 2*Cg*El(0,0)*El(0,1)*bf*bfs*expterm;
  A(0,0,0,2) = 2*Cg*El(0,0)*El(0,2)*bf*bfs*expterm;
  A(0,0,1,0) = 2*Cg*El(0,0)*El(1,0)*bf*bfs*expterm;
  A(0,0,1,1) = 2*Cg*El(0,0)*El(1,1)*bf*bt*expterm;
  A(0,0,1,2) = 2*Cg*El(0,0)*El(1,2)*bf*bt*expterm;
  A(0,0,2,0) = 2*Cg*El(0,0)*El(2,0)*bf*bfs*expterm;
  A(0,0,2,1) = 2*Cg*El(0,0)*El(2,1)*bf*bt*expterm;
  A(0,0,2,2) = 2*Cg*El(0,0)*El(2,2)*bf*bt*expterm;
  A(0,1,0,0) = 2*Cg*El(0,0)*El(0,1)*bf*bfs*expterm;
  A(0,1,0,1) = 2*Cg*pow(El(0,1), 2)*pow(bfs, 2)*expterm + Cg*bfs*expterm;
  A(0,1,0,2) = 2*Cg*El(0,1)*El(0,2)*pow(bfs, 2)*expterm;
  A(0,1,1,0) = 2*Cg*El(0,1)*El(1,0)*pow(bfs, 2)*expterm;
  A(0,1,1,1) = 2*Cg*El(0,1)*El(1,1)*bfs*bt*expterm;
  A(0,1,1,2) = 2*Cg*El(0,1)*El(1,2)*bfs*bt*expterm;
  A(0,1,2,0) = 2*Cg*El(0,1)*El(2,0)*pow(bfs, 2)*expterm;
  A(0,1,2,1) = 2*Cg*El(0,1)*El(2,1)*bfs*bt*expterm;
  A(0,1,2,2) = 2*Cg*El(0,1)*El(2,2)*bfs*bt*expterm;
  A(0,2,0,0) = 2*Cg*El(0,0)*El(0,2)*bf*bfs*expterm;
  A(0,2,0,1) = 2*Cg*El(0,1)*El(0,2)*pow(bfs, 2)*expterm;
  A(0,2,0,2) = 2*Cg*pow(El(0,2), 2)*pow(bfs, 2)*expterm + Cg*bfs*expterm;
  A(0,2,1,0) = 2*Cg*El(0,2)*El(1,0)*pow(bfs, 2)*expterm;
  A(0,2,1,1) = 2*Cg*El(0,2)*El(1,1)*bfs*bt*expterm;
  A(0,2,1,2) = 2*Cg*El(0,2)*El(1,2)*bfs*bt*expterm;
  A(0,2,2,0) = 2*Cg*El(0,2)*El(2,0)*pow(bfs, 2)*expterm;
  A(0,2,2,1) = 2*Cg*El(0,2)*El(2,1)*bfs*bt*expterm;
  A(0,2,2,2) = 2*Cg*El(0,2)*El(2,2)*bfs*bt*expterm;

  A(1,0,0,0) = 2*Cg*El(0,0)*El(1,0)*bf*bfs*expterm;
  A(1,0,0,1) = 2*Cg*El(0,1)*El(1,0)*pow(bfs, 2)*expterm;
  A(1,0,0,2) = 2*Cg*El(0,2)*El(1,0)*pow(bfs, 2)*expterm;
  A(1,0,1,0) = 2*Cg*pow(El(1,0), 2)*pow(bfs, 2)*expterm + Cg*bfs*expterm;
  A(1,0,1,1) = 2*Cg*El(1,0)*El(1,1)*bfs*bt*expterm;
  A(1,0,1,2) = 2*Cg*El(1,0)*El(1,2)*bfs*bt*expterm;
  A(1,0,2,0) = 2*Cg*El(1,0)*El(2,0)*pow(bfs, 2)*expterm;
  A(1,0,2,1) = 2*Cg*El(1,0)*El(2,1)*bfs*bt*expterm;
  A(1,0,2,2) = 2*Cg*El(1,0)*El(2,2)*bfs*bt*expterm;
  A(1,1,0,0) = 2*Cg*El(0,0)*El(1,1)*bf*bt*expterm;
  A(1,1,0,1) = 2*Cg*El(0,1)*El(1,1)*bfs*bt*expterm;
  A(1,1,0,2) = 2*Cg*El(0,2)*El(1,1)*bfs*bt*expterm;
  A(1,1,1,0) = 2*Cg*El(1,0)*El(1,1)*bfs*bt*expterm;
  A(1,1,1,1) = 2*Cg*pow(El(1,1), 2)*pow(bt, 2)*expterm + Cg*bt*expterm;
  A(1,1,1,2) = 2*Cg*El(1,1)*El(1,2)*pow(bt, 2)*expterm;
  A(1,1,2,0) = 2*Cg*El(1,1)*El(2,0)*bfs*bt*expterm;
  A(1,1,2,1) = 2*Cg*El(1,1)*El(2,1)*pow(bt, 2)*expterm;
  A(1,1,2,2) = 2*Cg*El(1,1)*El(2,2)*pow(bt, 2)*expterm;
  A(1,2,0,0) = 2*Cg*El(0,0)*El(1,2)*bf*bt*expterm;
  A(1,2,0,1) = 2*Cg*El(0,1)*El(1,2)*bfs*bt*expterm;
  A(1,2,0,2) = 2*Cg*El(0,2)*El(1,2)*bfs*bt*expterm;
  A(1,2,1,0) = 2*Cg*El(1,0)*El(1,2)*bfs*bt*expterm;
  A(1,2,1,1) = 2*Cg*El(1,1)*El(1,2)*pow(bt, 2)*expterm;
  A(1,2,1,2) = 2*Cg*pow(El(1,2), 2)*pow(bt, 2)*expterm + Cg*bt*expterm;
  A(1,2,2,0) = 2*Cg*El(1,2)*El(2,0)*bfs*bt*expterm;
  A(1,2,2,1) = 2*Cg*El(1,2)*El(2,1)*pow(bt, 2)*expterm;
  A(1,2,2,2) = 2*Cg*El(1,2)*El(2,2)*pow(bt, 2)*expterm;

  A(2,0,0,0) = 2*Cg*El(0,0)*El(2,0)*bf*bfs*expterm;
  A(2,0,0,1) = 2*Cg*El(0,1)*El(2,0)*pow(bfs, 2)*expterm;
  A(2,0,0,2) = 2*Cg*El(0,2)*El(2,0)*pow(bfs, 2)*expterm;
  A(2,0,1,0) = 2*Cg*El(1,0)*El(2,0)*pow(bfs, 2)*expterm;
  A(2,0,1,1) = 2*Cg*El(1,1)*El(2,0)*bfs*bt*expterm;
  A(2,0,1,2) = 2*Cg*El(1,2)*El(2,0)*bfs*bt*expterm;
  A(2,0,2,0) = 2*Cg*pow(El(2,0), 2)*pow(bfs, 2)*expterm + Cg*bfs*expterm;
  A(2,0,2,1) = 2*Cg*El(2,0)*El(2,1)*bfs*bt*expterm;
  A(2,0,2,2) = 2*Cg*El(2,0)*El(2,2)*bfs*bt*expterm;
  A(2,1,0,0) = 2*Cg*El(0,0)*El(2,1)*bf*bt*expterm;
  A(2,1,0,1) = 2*Cg*El(0,1)*El(2,1)*bfs*bt*expterm;
  A(2,1,0,2) = 2*Cg*El(0,2)*El(2,1)*bfs*bt*expterm;
  A(2,1,1,0) = 2*Cg*El(1,0)*El(2,1)*bfs*bt*expterm;
  A(2,1,1,1) = 2*Cg*El(1,1)*El(2,1)*pow(bt, 2)*expterm;
  A(2,1,1,2) = 2*Cg*El(1,2)*El(2,1)*pow(bt, 2)*expterm;
  A(2,1,2,0) = 2*Cg*El(2,0)*El(2,1)*bfs*bt*expterm;
  A(2,1,2,1) = 2*Cg*pow(El(2,1), 2)*pow(bt, 2)*expterm + Cg*bt*expterm;
  A(2,1,2,2) = 2*Cg*El(2,1)*El(2,2)*pow(bt, 2)*expterm;
  A(2,2,0,0) = 2*Cg*El(0,0)*El(2,2)*bf*bt*expterm;
  A(2,2,0,1) = 2*Cg*El(0,1)*El(2,2)*bfs*bt*expterm;
  A(2,2,0,2) = 2*Cg*El(0,2)*El(2,2)*bfs*bt*expterm;
  A(2,2,1,0) = 2*Cg*El(1,0)*El(2,2)*bfs*bt*expterm;
  A(2,2,1,1) = 2*Cg*El(1,1)*El(2,2)*pow(bt, 2)*expterm;
  A(2,2,1,2) = 2*Cg*El(1,2)*El(2,2)*pow(bt, 2)*expterm;
  A(2,2,2,0) = 2*Cg*El(2,0)*El(2,2)*bfs*bt*expterm;
  A(2,2,2,1) = 2*Cg*El(2,1)*El(2,2)*pow(bt, 2)*expterm;
  A(2,2,2,2) = 2*Cg*pow(El(2,2), 2)*pow(bt, 2)*expterm + Cg*bt*expterm;


  /*
  A(0,0,0,0) = 2*Cg*pow(El(0,0), 2)*pow(bf, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bf*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,0,1) = 2*Cg*El(0,0)*El(0,1)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,0,2) = 2*Cg*El(0,0)*El(0,2)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,1,0) = 2*Cg*El(0,0)*El(1,0)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,1,1) = 2*Cg*El(0,0)*El(1,1)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,1,2) = 2*Cg*El(0,0)*El(1,2)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,2,0) = 2*Cg*El(0,0)*El(2,0)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,2,1) = 2*Cg*El(0,0)*El(2,1)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,0,2,2) = 2*Cg*El(0,0)*El(2,2)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,0,0) = 2*Cg*El(0,0)*El(0,1)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,0,1) = 2*Cg*pow(El(0,1), 2)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,0,2) = 2*Cg*El(0,1)*El(0,2)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,1,0) = 2*Cg*El(0,1)*El(1,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,1,1) = 2*Cg*El(0,1)*El(1,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,1,2) = 2*Cg*El(0,1)*El(1,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,2,0) = 2*Cg*El(0,1)*El(2,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,2,1) = 2*Cg*El(0,1)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,1,2,2) = 2*Cg*El(0,1)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,0,0) = 2*Cg*El(0,0)*El(0,2)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,0,1) = 2*Cg*El(0,1)*El(0,2)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,0,2) = 2*Cg*pow(El(0,2), 2)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,1,0) = 2*Cg*El(0,2)*El(1,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,1,1) = 2*Cg*El(0,2)*El(1,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,1,2) = 2*Cg*El(0,2)*El(1,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,2,0) = 2*Cg*El(0,2)*El(2,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,2,1) = 2*Cg*El(0,2)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(0,2,2,2) = 2*Cg*El(0,2)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));

  A(1,0,0,0) = 2*Cg*El(0,0)*El(1,0)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,0,1) = 2*Cg*El(0,1)*El(1,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,0,2) = 2*Cg*El(0,2)*El(1,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,1,0) = 2*Cg*pow(El(1,0), 2)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,1,1) = 2*Cg*El(1,0)*El(1,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,1,2) = 2*Cg*El(1,0)*El(1,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,2,0) = 2*Cg*El(1,0)*El(2,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,2,1) = 2*Cg*El(1,0)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,0,2,2) = 2*Cg*El(1,0)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,0,0) = 2*Cg*El(0,0)*El(1,1)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,0,1) = 2*Cg*El(0,1)*El(1,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,0,2) = 2*Cg*El(0,2)*El(1,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,1,0) = 2*Cg*El(1,0)*El(1,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,1,1) = 2*Cg*pow(El(1,1), 2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,1,2) = 2*Cg*El(1,1)*El(1,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,2,0) = 2*Cg*El(1,1)*El(2,0)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,2,1) = 2*Cg*El(1,1)*El(2,1)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,1,2,2) = 2*Cg*El(1,1)*El(2,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,0,0) = 2*Cg*El(0,0)*El(1,2)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,0,1) = 2*Cg*El(0,1)*El(1,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,0,2) = 2*Cg*El(0,2)*El(1,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,1,0) = 2*Cg*El(1,0)*El(1,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,1,1) = 2*Cg*El(1,1)*El(1,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,1,2) = 2*Cg*pow(El(1,2), 2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,2,0) = 2*Cg*El(1,2)*El(2,0)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,2,1) = 2*Cg*El(1,2)*El(2,1)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(1,2,2,2) = 2*Cg*El(1,2)*El(2,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));

  A(2,0,0,0) = 2*Cg*El(0,0)*El(2,0)*bf*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,0,1) = 2*Cg*El(0,1)*El(2,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,0,2) = 2*Cg*El(0,2)*El(2,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,1,0) = 2*Cg*El(1,0)*El(2,0)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,1,1) = 2*Cg*El(1,1)*El(2,0)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,1,2) = 2*Cg*El(1,2)*El(2,0)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,2,0) = 2*Cg*pow(El(2,0), 2)*pow(bfs, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bfs*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,2,1) = 2*Cg*El(2,0)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,0,2,2) = 2*Cg*El(2,0)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,0,0) = 2*Cg*El(0,0)*El(2,1)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,0,1) = 2*Cg*El(0,1)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,0,2) = 2*Cg*El(0,2)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,1,0) = 2*Cg*El(1,0)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,1,1) = 2*Cg*El(1,1)*El(2,1)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,1,2) = 2*Cg*El(1,2)*El(2,1)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,2,0) = 2*Cg*El(2,0)*El(2,1)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,2,1) = 2*Cg*pow(El(2,1), 2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,1,2,2) = 2*Cg*El(2,1)*El(2,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,0,0) = 2*Cg*El(0,0)*El(2,2)*bf*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,0,1) = 2*Cg*El(0,1)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,0,2) = 2*Cg*El(0,2)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,1,0) = 2*Cg*El(1,0)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,1,1) = 2*Cg*El(1,1)*El(2,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,1,2) = 2*Cg*El(1,2)*El(2,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,2,0) = 2*Cg*El(2,0)*El(2,2)*bfs*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,2,1) = 2*Cg*El(2,1)*El(2,2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  A(2,2,2,2) = 2*Cg*pow(El(2,2), 2)*pow(bt, 2)*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2))) + Cg*bt*exp(pow(El(0,0), 2)*bf + bfs*(pow(El(0,1), 2) + pow(El(0,2), 2) + pow(El(1,0), 2) + pow(El(2,0), 2)) + bt*(pow(El(1,1), 2) + pow(El(1,2), 2) + pow(El(2,1), 2) + pow(El(2,2), 2)));
  */

  /*
  Cprod = tensor_product(invC, invC);

  for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
          for(int k=0; k<ndim; k++)
              for(int l=0; l<ndim; l++) {
                  _I(i,j,k,l) = 0.5*(invC(i, k) * invC(j, l) + invC(i, l) * invC(j, k));

              }
  for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
          for(int k=0; k<ndim; k++)
              for(int l=0; l<ndim; l++) {
                  Cp = p * J * (Cprod(i,j,k,l) - 2*_I(i,j,k,l));
                  Ck = K*J*J*Cprod(i,j,k,l);
                  A(i,j,k,l) += Cp + Ck;

              }

  */
}
