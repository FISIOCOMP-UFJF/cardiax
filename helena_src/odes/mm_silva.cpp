#include "mm_silva.hpp"

MMSilva::MMSilva() : CellModel(6) {
  var_names.insert( std::pair<int, std::string>(0, "u") );
  var_names.insert( std::pair<int, std::string>(1, "v") );
  var_names.insert( std::pair<int, std::string>(2, "w") );
  var_names.insert( std::pair<int, std::string>(3, "s") );
  var_names.insert( std::pair<int, std::string>(4, "ca") );
  var_names.insert( std::pair<int, std::string>(5, "ta") );
}

void MMSilva::init(double * values) const
{
  assert(values != nullptr);
  values[0] = 0.0; //u
  values[1] = 1.0; //v
  values[2] = 1.0; //w
  values[3] = 0.0; //s
  values[4] = 0.0; //ca
  values[5] = 0.0; //ta
}

void MMSilva::equation(const double time, const double *rY, double *rDY)
{
  // State variables
  const double u  = rY[0];
  const double v  = rY[1];
  const double w = rY[2];
  const double s = rY[3];
  const double ca = rY[4];
  const double ta = rY[5];
 
  double tauv=0, tauw=0, taus0=0, taus=0, tau0=0, vinf=0, winf=0, jfi=0, jso=0, jsi=0, c1, fat_ca, fina_ca, fat_ta, fina_ta;
  double du_dt=0, dv_dt=0, ds_dt=0, dw_dt=0, dca_dt=0, dta_dt = 0;

  double u0 =0.615234 ;
  double uu =1.100000;
  double thetav =0.012207 ;
  double thetaw =0.001172 ;
  double thetavlinha =0.169629 ;
  double theta0 =0.041602 ;
  double tauv1 =74.023438 ;
  double tauv2 =969.140625 ;
  double tauvmais =1.333984 ;
  double tauw1 = 92.382812 ;
  double tauw2 = 35.117188 ;
  double kw = 60.253906 ;
  double uw = 0.606445 ;
  double tauwmais = 160.839844 ;
  double taufi = 0.093750 ;
  double tau01 = 406.640625 ;
  double tau02 = 7.958984 ;
  double taus01 = 17.480469 ;
  double taus02 = 1.097656 ;
  double kso = 3.803711 ;
  double uso = 0.826172 ;
  double taus1 = 4.819336 ;
  double taus2 = 17.656250 ;
  double ks = 3.959961 ;
  double us = 0.193359 ;
  double tausi = 4.902344 ;
  double tauwinf = 0.793945 ;
  double winfstar = 0.518555 ;
  double gamma0 = 0.124023 ;
  double gamma2 = 9.503906 ;
  double p1 = 3.941406 ;
  double p2 = 0.003809 ;
  double gamma3 = 7.820312 ;
  double a = 0.375000 ;
  double CA_Distolico = 0 ;
  double gamma7 = 0.014746 ;
  double gamma4 = 0.830078 ;
  double gamma5 = 0.281250 ;
  double p11 = 0.037207 ;
  double p22 = 0.004883 ;
  double gamma33 = 3.066406 ;
  double aa = 0.092773;

  tauv = (1 - H(u, thetavlinha)) * tauv1 + H(u, thetavlinha) * tauv2;
  tauw = tauw1 + (tauw2 - tauw1) * (1 + tanh(kw * (u - uw))) / 2.0;
  taus0 = taus01 + (taus02 - taus01) * (1 + tanh(kso * (u - uso))) / 2.0;
  taus = (1 - H(u, thetaw)) * taus1 + H(u, thetaw) * taus2;
  tau0 = (1 - H(u, theta0)) * tau01 + H(u, theta0) * tau02;

  vinf = vnoinf(u, thetavlinha);
  winf = (1 - H(u, theta0)) * (1 - u / tauwinf) + H(u, theta0) * winfstar;

  jfi = -v * H(u, thetav) * (u - thetav) * (uu - u) / (taufi);
  //jjfi[i] = jfi
  jso = (u - u0) * (1 - H(u, thetaw)) / tau0 + H(u, thetaw) / taus0;
  //jjso[i] = jso
  jsi = -H(u, thetaw) * w * s / tausi;
  //jjsi[i] = jsi;

  fat_ca = (pow(-jfi, gamma0))/gamma2 ;
  //cout << time << " " << -jfi << " " << fat_ca << " " << pow(-jfi, gamma0) <<endl;
  //fat_ca = ((-jjfi[i])**gamma0)/gamma2;
  fina_ca = 1.*(((p1+p2)/2) -((p1-p2)/2)*tanh(-gamma3*(jsi - a)));

  fat_ta = gamma7*(pow(ca, gamma4));
  //fat_ta = gamma7*(ca**gamma4);

  //if((-jjfi[i])**(gamma5) > 0.001)
  if( pow( -jfi, gamma5) > 0.001 )
    c1 = 0;
  else
    c1 = 1;

  fina_ta = ((p11+p22)/2) + ((p11-p22)/2)*tanh(-gamma33*(ca - aa));

  du_dt = -(jfi + jso + jsi) + i_stim;
  dv_dt = (1 - H(u, thetav)) * (vinf - v) / tauv - H(u, thetav) * v / tauvmais;
  dw_dt = (1 - H(u, thetaw)) * (winf - w) / tauw - H(u, thetaw) * w / tauwmais;
  ds_dt = ((1 + tanh(ks * (u - us))) / 2.0 - s) / taus;

  dca_dt = fat_ca - fina_ca * (ca - CA_Distolico);
  dta_dt = c1*fat_ta - fina_ta*(ta);

  
  rDY[0] = du_dt;
  rDY[1] = dv_dt;
  rDY[2] = dw_dt;
  rDY[3] = ds_dt;
  rDY[4] = dca_dt;
  rDY[5] = dta_dt;
}

double MMSilva::H(double x, double y){
    if (x > y)
        return 1;
    else if(x < y)
        return 0;
    else
        return 0.5;
}
double MMSilva::vnoinf(double x, double y){
    if(x < y)
        return 1;
    else
        return 0;
}

