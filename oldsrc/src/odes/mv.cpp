#include "mv.hpp"

// Heaviside function
inline double H(double x) { return (x < 0) ? 0.0 : 1.0; }

MinimalVentricular::MinimalVentricular() : CellModel(4), paramset(false)
{
  var_names.insert( std::pair<int, std::string>(0, "u") );
  var_names.insert( std::pair<int, std::string>(1, "V") );
  var_names.insert( std::pair<int, std::string>(2, "w") );
  var_names.insert( std::pair<int, std::string>(3, "s") );
}
    
void MinimalVentricular::equation(const double time, const double * statevars, 
                                  double * values)
{
  // state variables
  const double u  = statevars[0];
  const double v  = statevars[1];
  const double w  = statevars[2];
  const double s  = statevars[3];
  //const double Ta = statevars[4];
  
  if(!paramset) init_params();

  // time constants - functions of u
  double tau_vb = (1-H(u-thetavb))*tauv1b + H(u-thetavb)*tauv2b;
  double tau_wb = tauw1b + (tauw2b - tauw1b)*(1+tanh(kwb*(u-uwb)))/2.0;
  double tau_so = tauso1 + (tauso2 - tauso1)*(1+tanh(kso*(u-uso)))/2.0;
  double tau_s  = (1-H(u-thetaw))*taus1 + H(u-thetaw)*taus2;
  double tau_o  = (1-H(u-thetao))*tauo1 + H(u-thetao)*tauo2;  
  double v_inf  = (u < thetavb) ? 1.0 : 0.0;
  double w_inf  = (1-H(u-thetao))*(1-(u/tau_winf)) + H(u-thetao)*w_infs; 
  
  // calculations
  double j_fi   = -v*H(u-thetav)*(u-thetav)*(uu-u) / tau_fi;
  double j_so   = (u-uo)*(1-H(u-thetaw))/tau_o + H(u-thetaw)/tau_so;
  double j_si   = -H(u-thetaw)*w*s/tau_si;
  double j_stim = i_stim;
    
  values[0] = -(j_fi + j_so + j_si + j_stim); // V(mv) = 85.7*u - 84
  values[1] = (1.0-H(u-thetav))*(v_inf-v)/tau_vb - H(u-thetav)*v/tauvp;
  values[2] = (1.0-H(u-thetaw))*(w_inf-w)/tau_wb - H(u-thetaw)*w/tauwp;
  values[3] = ((1.0+tanh(ks*(u-us)))/2.0 - s)/tau_s;
  
  //values[4] = 0.075 * (0.55*s - 0.3*Ta);
}

void MinimalVentricular::init(double * values) const
{
  values[0] = 0.0; // u
  values[1] = 1.0; // v
  values[2] = 1.0; // w
  values[3] = 0.0; // s
  //values[4] = 0.0; // Ta
}
   
void MinimalVentricular::init_params()
{
  paramset = true;
  
  if(type == MCELL)
  {
    uo       = 0.0;
    uu       = 1.61;
    thetav   = 0.3;
    thetaw   = 0.13;
    thetavb  = 0.1;
    thetao   = 0.005;
    tauv1b   = 80;
    tauv2b   = 1.4506;  
    tauvp    = 1.4506;
    tauw1b   = 70;
    tauw2b   = 8;
    kwb      = 200;
    uwb      = 0.016;
    tauwp    = 280;
    tau_fi   = 0.078;
    tauo1    = 410;
    tauo2    = 7;
    tauso1   = 91;
    tauso2   = 0.8;
    kso      = 2.1;
    uso      = 0.6;
    taus1    = 2.7342;
    taus2    = 4;
    ks       = 2.0994;  
    us       = 0.9087;
    tau_si   = 3.3849;
    tau_winf = 0.01;
    w_infs   = 0.5; 
  }
  else if(type == EPI)
  {
    uo       = 0.0;
    uu       = 1.55;
    thetav   = 0.3;
    thetaw   = 0.13;
    thetavb  = 0.006;
    thetao   = 0.006;
    tauv1b   = 60;
    tauv2b   = 1150;  
    tauvp    = 1.4506;
    tauw1b   = 60;
    tauw2b   = 15;
    kwb      = 65;
    uwb      = 0.03;
    tauwp    = 200;
    tau_fi   = 0.11;
    tauo1    = 400;
    tauo2    = 6;
    tauso1   = 30.0181;
    tauso2   = 0.9957;
    kso      = 2.0458;
    uso      = 0.65;
    taus1    = 2.7342;
    taus2    = 16;
    ks       = 2.0994;  
    us       = 0.9087;
    tau_si   = 1.8875;
    tau_winf = 0.07;
    w_infs   = 0.94;  
  }
  else if(type == ENDO)
  {  
    uo       = 0.0;
    uu       = 1.56;
    thetav   = 0.3;
    thetaw   = 0.13;
    thetavb  = 0.2;
    thetao   = 0.006;
    tauv1b   = 75;
    tauv2b   = 10;  
    tauvp    = 1.4506;
    tauw1b   = 6;
    tauw2b   = 140;
    kwb      = 200;
    uwb      = 0.016;
    tauwp    = 280;
    tau_fi   = 0.1;
    tauo1    = 470;
    tauo2    = 6;
    tauso1   = 40;
    tauso2   = 1.2;
    kso      = 2;
    uso      = 0.65;
    taus1    = 2.7342;
    taus2    = 2;
    ks       = 2.0994;  
    us       = 0.9087;
    tau_si   = 2.9013;
    tau_winf = 0.0273;
    w_infs   = 0.78; 
  } 
}
