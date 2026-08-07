#include "rice_ten_tusscher.hpp"

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) < (b) ? (b) : (a))
#define sign(a) ((a) < (0.) ? (-1.) : (1.))
#define heav(a) ((a) < (0.) ? (0.) : (1.))
#define abs(a) ((a) < (0.) ? (-a) : (a))



#ifdef HF
RiceTenTusscher::RiceTenTusscher() : CellModel(30)
#endif
#ifndef HF
RiceTenTusscher::RiceTenTusscher() : CellModel(28)
#endif
{
  // state variable names
  var_names.insert( std::pair<int, std::string>(0, "V") );
  var_names.insert( std::pair<int, std::string>(1, "Xr1") );
  var_names.insert( std::pair<int, std::string>(2, "Xr2") );
  var_names.insert( std::pair<int, std::string>(3, "Xs") );
  var_names.insert( std::pair<int, std::string>(4, "m") );
  var_names.insert( std::pair<int, std::string>(5, "h") );
  var_names.insert( std::pair<int, std::string>(6, "j") );
  var_names.insert( std::pair<int, std::string>(7, "d") );
  var_names.insert( std::pair<int, std::string>(8, "f") );
  var_names.insert( std::pair<int, std::string>(9, "fCa") );
  var_names.insert( std::pair<int, std::string>(10, "s") );
  var_names.insert( std::pair<int, std::string>(11, "r") );
  var_names.insert( std::pair<int, std::string>(12, "g") );
  var_names.insert( std::pair<int, std::string>(13, "Ca_SR") );
  var_names.insert( std::pair<int, std::string>(14, "Ca_i_total") );
  var_names.insert( std::pair<int, std::string>(15, "Na_i") );
  var_names.insert( std::pair<int, std::string>(16, "K_i") );
  var_names.insert( std::pair<int, std::string>(17, "TRPNCaL") );
  var_names.insert( std::pair<int, std::string>(18, "TRPNCaH") );
  var_names.insert( std::pair<int, std::string>(19, "N_NoXB") );
  var_names.insert( std::pair<int, std::string>(20, "P_NoXB") );
  var_names.insert( std::pair<int, std::string>(21, "N") );
  var_names.insert( std::pair<int, std::string>(22, "XBprer") );
  var_names.insert( std::pair<int, std::string>(23, "XBpostr") );
  var_names.insert( std::pair<int, std::string>(24, "xXBprer") );
  var_names.insert( std::pair<int, std::string>(25, "xXBpostr") );
  var_names.insert( std::pair<int, std::string>(26, "SL") );
  var_names.insert( std::pair<int, std::string>(27, "intf") );

  // setup
  type = EPI;
  isometric = true;

  monitored.push_back( &active );
  //monitored.push_back( &SL );
  //monitored.push_back( &INaCa );

  active = 0.00059167;
  
  // variables for RL integration in Explicit Euler
  for(int i=1; i<12; i++) if(i != 9) rlvars.insert(i);


#ifdef HF  
  rlvars.insert(28);
  rlvars.insert(29);
#endif

}

void RiceTenTusscher::init(double * values) const
{
  /*
  // OLD - OK for isometric
  values[0] =-85.23;                    // V
  values[1] = 0.00621;                  // Xr1
  values[2] = 0.4712;                   // Xr2
  values[3] = 0.0095;                   // Xs
  values[4] = 0.00172;                  // m
  values[5] = 0.7444;                   // h
  values[6] = 0.7045;                   // j
  values[7] = 3.373e-5;                 // d
  values[8] = 0.7888;                   // f
  values[9] =  1.0;                     // fCa
  values[10] = 1.0;                     // s
  values[11] = 2.42e-8;                 // r
  values[12] = 1.0;                     // g
  values[13] = 0.2;                     // Ca_SR
  values[14] = 9.28194487e-03;          // Ca_i_total
  values[15] = 8.604;                   // Na_i
  values[16] = 136.89;                  // K_i
  values[17] = 1.47730085e-02;          // TRPNCaL
  values[18] = 1.30660966e-01;          // TRPNCaH
  values[19] = 9.99999959e-01;          // N_NoXB
  values[20] = 4.07437174e-08;          // P_NoXB
  values[21] = 9.99997835e-01;          // N
  values[22] = 3.04949649e-07;          // XBprer
  values[23] = 1.81017564e-06;          // XBpostr
  values[24] = 3.41212829e-08;          // xXBprer
  values[25] = 7.00005395e-03;          // xXBpostr
  values[26] = (isometric) ? 2.4 : 1.9; // SL
  values[27] = -4.51134525e-06;         // intf
  */  

//#ifndef HF
  values[0] = -8.629169987876476e+01;  // V - ten Tusscher model variables
  values[1] = 1.817573054244821e-04;   // Xr1
  values[2] = 4.822125371672458e-01;   // Xr2
  values[3] = 2.998779317686131e-03;   // Xs
  values[4] = 1.368501861611859e-03;   // m
  values[5] = 7.728393292625330e-01;   // h
  values[6] = 3.694904454167577e-01;   // j
  values[7] = 1.962096961554300e-05;   // d
  values[8] = 9.998430021390249e-01;   // f
  values[9] = 1.006126715016802e+00;   // fCa
  values[10] = 6.699733924504092e-01;  // s
  values[11] = 2.024681561943644e-08;  // r
  values[12] = 9.999492474604431e-01;  // g
  values[13] = 5.205023026469432e-01;  // Ca_SR
  values[14] = 9.281944866155925e-03;  // Ca_i_total
  values[15] = 1.114903806343864e+01;  // Na_i
  values[16] = 1.383870498592104e+02;  // K_i
  values[17] = 1.311635148577327e-02;  // TRPNCaL - Rice myofilament variables
  values[18] = 1.175428809165677e-01;  // TRPNCaH
  values[19] = 9.999995709556512e-01;  // N_NoXB
  values[20] = 4.290431881685032e-07;  // P_NoXB
  values[21] = 9.996687602297853e-01;  // N
  values[22] = 2.211689363408625e-04;  // XBprer
  values[23] = 1.079302937362542e-04;  // XBpostr
  values[24] = 1.007048444140799e-18;  // xXBprer
  values[25] = 7.000000000000143e-03;  // xXBpostr
  values[26] = 1.9;                    // SL
  values[27] = -3.637088203841350e+07; // intf    
//#endif
  
#ifdef HF
  if (type==ENDO)
  {
    cout << "IC HF ENDO" << endl;
    values[0] = -8.528928e+01;
    values[1] = 2.155566e-04;
    values[2] = 4.717913e-01;
    values[3] = 3.220768e-03;
    values[4] = 1.693961e-03;
    values[5] = 7.464877e-01;
    values[6] = 7.450248e-01;
    values[7] = 2.242698e-05;
    values[8] = 9.995748e-01;
    values[9] = 1.005842e+00;
    values[10] = 6.053618e-01;
    values[11] = 2.393213e-08;
    values[12] = 9.999134e-01;
    values[13] = 8.617414e-02;
    values[14] = 1.008992e-02;
    values[15] = 1.328551e+01;
    values[16] = 1.368598e+02;
    values[17] = 1.432423e-02;
    values[18] = 1.275423e-01;
    values[19] = 9.999999e-01;
    values[20] = 1.282613e-07;
    values[21] = 9.999881e-01;
    values[22] = 5.375208e-06;
    values[23] = 6.344522e-06;
    //values[24] = 1.005483e-18;
    //values[25] = 7.000000e-03;
    values[24] = 1.007048444140799e-18;  // xXBprer
    values[25] = 7.000000000000143e-03;  // xXBpostr
    values[26] = 1.900000e+00;
    values[27] = -3.637284e+07;
    values[28] = 1.472719e-03;
    values[29] = 2.571574e-01;
  }
  else if(type==MCELL)
  {
    cout << "IC HF MID" << endl;
    values[0] = -8.481683e+01;
    values[1] = 2.642653e-04;
    values[2] = 4.668872e-01;
    values[3] = 3.331036e-03;
    values[4] = 1.872532e-03;
    values[5] = 7.332868e-01;
    values[6] = 7.300959e-01;
    values[7] = 2.388534e-05;
    values[8] = 9.989570e-01;
    values[9] = 1.005031e+00;
    values[10] = 9.999977e-01;
    values[11] = 2.589578e-08;
    values[12] = 9.996975e-01;
    values[13] = 1.398963e-01;
    values[14] = 1.225003e-02;
    values[15] = 1.280699e+01;
    values[16] = 1.372002e+02;
    values[17] = 1.758578e-02;
    values[18] = 1.526399e-01;
    values[19] = 9.999995e-01;
    values[20] = 5.039252e-07;
    values[21] = 9.999722e-01;
    values[22] = 1.521383e-05;
    values[23] = 1.198677e-05;
    //values[24] = 1.005483e-18;
    //values[25] = 7.000000e-03;
    values[24] = 1.007048444140799e-18;  // xXBprer
    values[25] = 7.000000000000143e-03;  // xXBpostr
    values[26] = 1.900000e+00;
    values[27] = -3.637315e+07;
    values[28] = 1.593312e-03;
    values[29] = 2.342089e-01;
  }
  else if(type==EPI)
  {    
    cout << "IC HF EPI" << endl;
    values[0] = -8.455109e+01;
    values[1] = 2.487432e-04;
    values[2] = 4.641325e-01;
    values[3] = 3.394631e-03;
    values[4] = 1.980930e-03;
    values[5] = 7.256835e-01;
    values[6] = 7.234677e-01;
    values[7] = 2.474674e-05;
    values[8] = 9.995314e-01;
    values[9] = 1.005301e+00;
    values[10] = 9.999975e-01;
    values[11] = 2.706767e-08;
    values[12] = 9.997927e-01;
    values[13] = 1.477493e-01;
    values[14] = 1.155700e-02;
    values[15] = 1.322450e+01;
    values[16] = 1.367668e+02;
    values[17] = 1.652898e-02;
    values[18] = 1.444025e-01;
    values[19] = 9.999996e-01;
    values[20] = 3.305227e-07;
    values[21] = 9.999809e-01;
    values[22] = 1.016265e-05;
    values[23] = 8.524538e-06;
    //values[24] = 1.005483e-18;
    //values[25] = 7.000000e-03;
    values[24] = 1.007048444140799e-18;  // xXBprer
    values[25] = 7.000000000000143e-03;  // xXBpostr
    values[26] = 1.900000e+00;
    values[27] = -3.637282e+07;
    values[28] = 1.665350e-03;
    values[29] = 2.339779e-01;
  }
#endif  

}

void RiceTenTusscher::equation(const double time, const double * statevars, 
                               double * values)
{
  // copy values
  V = statevars[0];
  Xr1 = statevars[1];
  Xr2 = statevars[2];
  Xs = statevars[3];
  m = statevars[4];
  h = statevars[5];
  j = statevars[6];
  d = statevars[7];
  f = statevars[8];
  fCa = statevars[9];
  s = statevars[10];
  r = statevars[11];
  g = statevars[12];
  Ca_SR = statevars[13];
  Ca_i_total = statevars[14];
  Na_i = statevars[15];
  K_i = statevars[16];
  TRPNCaL = statevars[17];
  TRPNCaH = statevars[18];
  N_NoXB = statevars[19];
  P_NoXB = statevars[20];
  N = statevars[21];
  XBprer = statevars[22];
  XBpostr = statevars[23];
  xXBprer = statevars[24];
  xXBpostr = statevars[25];
  SL = statevars[26];
  intf = statevars[27];
    
  // setup values
  const double R = 8314.47200;
  const double T = 310.0;
  const double F = 96485.3415;
  const double Cm = 0.185;  

  // External concentrations
  const double K_o  = 5.4;
  const double Ca_o = 2.0;
  const double Na_o = 140.0;

  // Intracellular volumes
  const double Vc = 1.6404e-02;
  const double Vsr = 1.094e-03;

  // Parameters for currents - maximal conductances
  double g_K1  = 5.405;
  double g_Kr  = 0.096;  
  double g_Ks  = 0.245;
  double g_Na  = 14.838;
  double g_bna = 2.9e-04;
  double g_CaL = 1.75e-04;
  double g_bca = 5.92e-04;  
  double g_to  = 0.294;
  double g_pCa = 0.825; 
  double g_pK  = 1.46e-02;
  
#ifdef HF  
  //double g_NaL = 0.015; // Trenor PLOS
  double g_NaL = 0.055;   // Trenor com TP06
#endif

  // Other parameters
  double P_NaK = 1.362;
  double P_kna = 0.03;
  double K_mk = 1.0;
  double K_mNa = 40.0;
  double K_NaCa = 1000.0;
  
  double Km_Nai   = 87.5;
  double Km_Ca    = 1.38; 
  double K_sat    = 0.1;
  double K_pCa    = 5.0e-04;  
  double alpha    = 2.5;
  double a_rel    = 1.6464e-02;
  double b_rel    = 0.25;
  double c_rel    = 8.232e-03;
  double gamma_tt = 0.35;  
  double Vmax_up  = 0.000425;  
  double K_up     = 2.5e-04;    
  double V_leak   = 8.0e-05; 
  double tau_g    = 2.0;    
  double Buf_c    = 1.31477e-01; // adjusted    
  double K_buf_c  = 1.0e-03;
  double Buf_sr   = 10.0;
  double K_buf_sr = 0.3;

  // ---------------------------------------------------------------------------
  // Rice parameters
  // ---------------------------------------------------------------------------
  
  const double TmpC = 37.0;
  
  // Sarcomere Geometry
  const double SLmax = 2.4;
  const double SLmin = 1.4;
  const double len_thick = 1.65;
  const double len_thin  = 1.2;
  const double len_hbare = 0.1;
  
  // Temperature Dependence
  const double Qkon = 1.5;
  const double Qkoff = 1.3;
  const double Qkn_p = 1.6;
  const double Qkp_n = 1.6;
  const double Qfapp = 6.25;
  const double Qgapp = 2.5;
  const double Qhf = 6.25;
  const double Qhb = 6.25;
  const double Qgxb = 6.25;
  
  const double kon = 4.96439500e-02;    
  const double koffL = 0.25;
  const double koffH = 0.025;
  const double perm50 = 0.5;
  const double nperm = 15.0;
  const double kn_p = 0.60778;  
  const double kp_n = 0.01602985;
  const double koffmod = 1.00634;  
  
  const double fapp = 4.80783;
  const double gapp = 0.093443;
  const double gslmod = 6.0;
  const double hf = 0.01027154;  
  const double hfmdc = 5.0;
  const double hb = 0.03403732;  
  const double hbmdc = 0.0;
  const double gxb = 0.03026968;
  const double sigmap = 8.0;
  const double sigman = 1.0;  
  const double xbmodsp = 0.2; 
  const double x_0 = 7.0e-03;
  const double xPsi = 2.0;  
  const double SLrest = 1.85;
  const double PCon_t = 0.002;
  const double PExp_t = 10.0;
  const double PCon_c = 0.02;
  const double PExp_c = 70.0;
  const double SL_c = 2.25;  
  const double SLset = 1.9;
  const double SEon = 1.0;
  
  //const double massf = 50.0;
  //const double visc = 3.0;
  const double KSE = 1.0;
  const double kxb = 120.0;
  const double Trop_conc = 70.0;
  
  // setup type
  g_Ks = 0.245;
  if (type == EPI)
    g_Ks = 0.245;
  else if (type == ENDO)
    g_Ks = 0.245;
  else if (type == MCELL)
    g_Ks = 0.062;

  if (type != ENDO)
    g_to = 0.294;
  else
    g_to = 0.073;

  // --------------------------------------------------------------------------
  // HYPERTROPHIC-HF MODEL
  // --------------------------------------------------------------------------      
#ifdef HF
  
  /*
  // Vivek Iyer 2011
  if(type==ENDO) 
  {    
    g_K1    *= 0.9;
    Vmax_up *= 0.5; 
    K_NaCa  *= 1.60;
    g_to    *= 1.0;
    g_CaL   *= 0.45;
    g_NaL   *= 1.1;
  }
  else if (type==MCELL)
  {
    g_K1    *= 0.9;
    Vmax_up *= 0.35;
    K_NaCa  *= 1.60;
    g_to    *= 0.7;
    g_CaL   *= 0.8;
    g_NaL   *= 1.3;
  }
  else if (type==EPI)
  {
    g_K1    *= 0.9;
    Vmax_up *= 0.3;
    K_NaCa  *= 2.0;
    g_to    *= 0.7;
    g_CaL   *= 1.0;
    g_NaL   *= 1.1;
  }   
  */  

  // MUITO BOM 
  //g_K1  *= 0.8; g_Ks *= 0.5;  V_leak *= 1.45; g_to *= 0.64;  Vmax_up *= 0.5;
  //g_NaL *= 5.0;
  //if(type==ENDO)        K_NaCa  *= 1.75;
  //else if (type==MCELL) K_NaCa  *= 1.75;
  //else if (type==EPI)   K_NaCa  *= 2.0;
  
  // Trevor - PLOS ------------------------------------------------------------
  g_K1    *= 0.67;
  P_NaK   *= 0.9;
  g_bna    = 0.0;
  g_bca   *= 1.53;
  V_leak  *= 5.00; 
  g_NaL   *= 2.0;
    
  if(type==ENDO)
  {    
    K_NaCa  *= 1.65;    
    Vmax_up *= 0.55;    
    g_CaL   *= 0.65;
  }
  else if (type==MCELL)
  {   
    K_NaCa  *= 1.65;
    Vmax_up *= 0.6;    
    g_to    *= 0.30;
    g_CaL   *= 0.8;
  }
  else if (type==EPI)
  {
    K_NaCa  *= 2.0;
    Vmax_up *= 0.7;
    g_to    *= 0.30;
    g_CaL   *= 1;
  }   
#endif    
 
  // Calculations
  const double RTOF = (R*T)/F;
  const double TMP3710 = (TmpC-37.0)/10.0;
  double VFORT = (V*F)/(R*T);
  double CaSRsquare = Ca_SR * Ca_SR;

  // Reversal potentials
  double E_Na = RTOF * log(Na_o/Na_i);
  double E_K  = RTOF * log(K_o/K_i);
  double E_Ks = RTOF * log(((K_o + (P_kna*Na_o)) / ( K_i + (P_kna * Na_i))));

  double xr1_inf   = 1.0/(1.0 + exp(((-26.0)-V)/7.0));
  double alpha_xr1 = 450.0/(1.0 + exp(((-45.0)-V)/10.0));
  double beta_xr1  = 6.0/(1.0 + exp(((V + 30.0)/11.50)));
  double tau_xr1   = alpha_xr1 * beta_xr1;
  
  double xr2_inf   = 1.0/(1.0 + exp(((V + 88.0)/24.0)));
  double alpha_xr2 = 3.0/(1.0 + exp((((-60.0)-V)/20.0)));
  double beta_xr2  = 1.1200/(1.0 + exp(((V-60.0)/20.0)));
  double tau_xr2   = alpha_xr2 * beta_xr2;
  
  double xs_inf   = 1.0/(1.0 + exp((((-5.0)-V)/14.0)));
  double alpha_xs = 1.1e+03/pow((1.0+exp((((-1.0e+01)-V)/6.0))),1.0/2.0);
  double beta_xs  = 1.0/(1.0 + exp(((V-6.0e+01)/2.0e+01)));
  double tau_xs   = 1.0 * alpha_xs * beta_xs;
  
  double m_inf   = 1.0/pow((1.0 + exp((((-56.86)-V)/9.03))),2.0);
  double alpha_m = 1.0/(1.0 + exp((((-60.0)-V)/5.0)));
  double beta_m  = (0.1/(1.0 + exp(((V+35.0)/5.0))))+(1.0e-01/(1.+exp(((V-50.0)/200.0))));
  double tau_m   = alpha_m * beta_m;
  
  double h_inf   = 1.0/pow((1.0 + exp(((V+71.55)/7.43))),2.0);
  double alpha_h = (V < -40.) ? 0.057 * exp( -(V+80.)/6.8 ) : 0.0;  
  double beta_h  = (V < -40.) ? ((2.7*exp((0.079*V))) + (3.1e+05 * exp(0.3485*V))) : ((0.77/(0.13*(1.+exp(((V + 10.66)/(-11.1)))))));
  double tau_h   = 1.0/(alpha_h + beta_h);
  
  double j_inf   = 1.0/pow((1.0 + exp(((V+71.55)/7.43))),2.0);
  double alpha_j = (V < -40.) ? ( ((-2.5428e+04)*exp(0.2444*V)-(6.9480e-06*exp(-0.04391*V)) ) * (V+37.78) )/(1. + exp(0.311*(V+79.23))) : 0.0;
  double beta_j  = (V < -40.) ? ( ((2.4240e-02*exp(((-1.0520e-02)*V)))/(1.+exp(((-0.13780)*(V+40.14)))))) : (0.6*exp(0.057*V))/(1.+exp((-0.1)*(V+32.)));
  double tau_j   = 1.0/(alpha_j + beta_j);       
  
  double d_inf   = 1.0/(1.0 + exp((((-5.0)-V)/7.5)));
  double alpha_d = (1.4/(1.0 + exp((((-35.)-V)/13.0)))) + 0.25;
  double beta_d  = 1.4/(1.0 + exp(((V+5.0)/5.0)));
  double gamma_d = 1.0/(1.0 + exp(((5.0e+01-V)/2.0e+01)));
  double tau_d   = (1.0 * alpha_d * beta_d) + gamma_d;
  
  double f_inf = 1.0/(1.+exp(((V+20.)/7.)));
  double tau_f = (1125.*exp(((-pow((V + 27.),2))/240.)))+80.+(165./(1.0+exp(((25.-V)/10.))));

  double s_inf, tau_s;  
  if (type != ENDO) // For epicardial and M cells
  {    
    s_inf = 1.0 / (1.0 + exp((V + 20.0) / 5.0));
    tau_s = (85.0*exp((-pow(V + 45.0, 2))/320.0)) + (5.0/(1.0 + exp((V - 20.0) / 5.0))) + 3.0;
  }
  else // For endocardial cells
  {
    s_inf = 1.0 / (1.0 + exp((V + 28.0) / 5.0));
    tau_s = (1000.0*exp((-pow(V + 67.0, 2))/1000.0)) + 8.0;
  }
  
  double r_inf = 1.0/(1.0 + exp((20.0-V)/6.0) );
  double tau_r = (9.5 * exp(((-pow((V + 40.0),2.0))/1800.0))) + 0.8;
  
  double alpha_K1 = 0.1/(1.0+exp(0.06*(V-E_K-200.0)));
  double beta_K1 = ((3.*exp(2.0e-04*(V-E_K+100.))) 
                 + (1.*exp(0.1*((V-E_K)-10.))))/(1.+exp(-0.5*(V-E_K)));
  double xK1_inf = alpha_K1/(alpha_K1 + beta_K1);
  
//  double i_stimAmplitude = -52.0;
//  double i_stimStart = 50.0;
//  double i_stimEnd = 500000.0;
//  double i_stimPeriod = 1000.0;
//  double i_stimPulseDuration = 1.0;
//  if(((time>=i_stimStart)&&(time<=i_stimEnd)&&
//     (((time-i_stimStart)-(floor(((time-i_stimStart)/i_stimPeriod))
//          *i_stimPeriod))<=i_stimPulseDuration)))
//    Istim = i_stimAmplitude;
//  else
//    Istim=0;  
   
  Istim = i_stim;
  
  IK1 = (g_K1*xK1_inf*pow((K_o/5.4000),1.0/2.0)*(V-E_K));
  INaK  = P_NaK*(K_o*Na_i)/((K_o+K_mk)*(Na_i+K_mNa)*(1.+0.1245*exp(-0.1*VFORT)+0.0353*exp(-VFORT)));
  IKr = g_Kr*pow(K_o/5.4,1./2.)*Xr1*Xr2*(V-E_K);
  IKs = g_Ks * pow(Xs,2.0) * (V - E_Ks);
  Ito = g_to*r*s*(V-E_K);
  IpK = (g_pK*(V-E_K))/(1.0+exp(((25.0-V)/5.98)));  
  INa = g_Na*pow(m,3.0)*h*j*(V-E_Na);
  IbNa = g_bna*(V-E_Na);
  Irel  = ((((a_rel*CaSRsquare)/((pow(b_rel,2.0)+CaSRsquare))+c_rel)*d*g));
  
#ifdef HF  
  double mL = statevars[28];
  double hL = statevars[29];    
  double hL_inf   = 1.0/(1.0+exp((V+91.0)/(6.1)));
  double tau_hL   = 233; 
  double alpha_mL = (0.32 * (V+47.13))/(1.0-exp(-0.1*(V+47.13)));
  double beta_mL  = 0.08 * exp(-V/11.0);
  double INaL = g_NaL * pow(mL,3) * hL * (V-E_Na);    
#endif
#ifndef HF
  double INaL = 0.0;
#endif    
   
  // Rice myofilament model calculations
  double P = (((1.0-N)-XBprer)-XBpostr);
  
  // Compute single overlap fractions
  double sovr_ze = min(len_thick/2,SL/2);
  double sovr_cle =  max(SL/2-(SL-len_thin),len_hbare/2);
  double len_sovr = sovr_ze - sovr_cle;
  double SOVFThick = (len_sovr*2.0)/(len_thick-len_hbare);
  double SOVFThin = len_sovr/len_thin;  

  // Compute combined Ca binding to high- (w/XB) and low- (no XB) sites  
  double Tropreg = (((1.0-SOVFThin)*TRPNCaL)+(SOVFThin*TRPNCaH));
  double permtot = pow(fabs((1.0/(1.0+pow((perm50/Tropreg),nperm)))),1.0/2.0);
  double inprmt = min(1/permtot, 100);  
  
  // Adjustments for Ca activation, temperature, SL, stress and strain
  double konT = kon * pow(Qkon,TMP3710);
  double koffLT = koffL * koffmod * pow(Qkoff,TMP3710);
  double koffHT = koffH * koffmod * pow(Qkoff,TMP3710);
  double kn_pT = kn_p * permtot * pow(Qkn_p,TMP3710);
  double kp_nT = kp_n * inprmt * pow(Qkp_n, TMP3710);
  double fappT = fapp * xbmodsp * pow(Qfapp,TMP3710);
  double gapslmd = (1.0+((1.0-SOVFThick)*gslmod));
  double gappT = gapp * gapslmd * xbmodsp * pow(Qgapp,TMP3710);
  double hfmd = exp((((-xXBprer)/fabs(xXBprer))*hfmdc*pow((xXBprer/x_0),2.)));
  double hbmd = exp((((xXBpostr-x_0)/fabs((xXBpostr-x_0)))*hbmdc*pow(((xXBpostr-x_0)/x_0),2.)));
  double hfT = hf * hfmd * xbmodsp * pow(Qhf,TMP3710);
  double hbT = hb * hbmd * xbmodsp * pow(Qhb,TMP3710);   
  double gxbmd = heav(x_0-xXBpostr)*exp(sigmap*((x_0-xXBpostr)/x_0)*((x_0-xXBpostr)/x_0))
               + (1.-heav(x_0-xXBpostr)*exp(sigman*(((xXBpostr-x_0)/x_0)*(xXBpostr-x_0)/x_0)));
  double gxbT = gxb * gxbmd * xbmodsp * pow(Qgxb,TMP3710);

  // regulation and crossbridge cycling state derivatives
  double dXBprer = (((fappT*P)+(hbT*XBpostr))-((gappT*XBprer)+(hfT*XBprer)));
  double dXBpostr = ((hfT*XBprer)-((hbT*XBpostr)+(gxbT*XBpostr)));
    
  // steady-state fractions in XBprer and XBpostr using King-Altman rule
  double SSXBprer = (((hb*fapp)+(gxb*fapp))/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));
  double SSXBpostr = ((fapp*hf)/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));
  
  // normalization for scaling active and passive force (maximal force)
  double Fnordv = kxb * x_0 * SSXBpostr;
  
  // Calculate Forces (active, passive, preload, afterload)
  double force = (kxb*SOVFThick*((xXBpostr*XBpostr)+(xXBprer*XBprer)));
  active = force/Fnordv; 
  
  double ppforce_t = sign(SL-SLrest)*PCon_t*(exp(PExp_t*abs(SL-SLrest))-1);
  double ppforce_c = heav(SL-SL_c)*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1);
  double ppforce = ppforce_t + ppforce_c;
  double preload = (fabs((SLset-SLrest))/(SLset-SLrest))*PCon_t*(exp((PExp_t*fabs((SLset-SLrest))))-1.);
  double afterload = (SEon == 1.0) ? (KSE * (SLset - SL)) : 0.0;  
  double dintf = (preload + afterload) - (ppforce + active);
  
  //dSL = ifnumber_14(); 
  const double stretch = 1.0;
  if(!isometric && ((SL <= SLmax) && (SL > SLmin)) )
    dSL = SL*(stretch-1.0)/0.05;
   //dSL = (((intf + ((SLset - SL)*visc))/massf));
  else
    dSL = 0.0;

  double dutyprer = ((hbT*fappT)+(gxbT*fappT))/((fappT*hfT)+(gxbT*hfT)+(gxbT*gappT)+(hbT*fappT)+(hbT*gappT)+(gxbT*fappT));
  double dutypostr = (fappT*hfT)/((fappT*hfT)+(gxbT*hfT)+(gxbT*gappT)+(hbT*fappT)+(hbT*gappT)+(gxbT*fappT));     
  double dxXBprer = ((dSL/2.0)+((xPsi/dutyprer)*((fappT*(-xXBprer))+(hbT*(xXBpostr-(x_0+xXBprer))))));
  double dxXBpostr = (dSL/2.0) + ((xPsi/dutypostr)*hfT*(xXBprer+x_0-xXBpostr)); 
  
  // Ca buffering by low-affinity troponin C (LTRPNCa)
  double FrSBXB = ((XBpostr+XBprer)/(SSXBpostr+SSXBprer));
  double dFrSBXB = (dXBpostr+dXBprer)/(SSXBpostr+SSXBprer);   
  
  double dsovr_ze = -(dSL/2.0)*heav(len_thick-SL);
  double dsovr_cle = -(dSL/2.0)*heav((2.0*len_thin-SL)-len_hbare);
  double dlen_sovr = dsovr_ze - dsovr_cle;   
  double dSOVFThin = (dlen_sovr/len_thin);
  double dSOVFThick = (2.*dlen_sovr)/(len_thick-len_hbare);      
  
  TropTot = (Trop_conc*(((1.-SOVFThin)*TRPNCaL)+(SOVFThin*((FrSBXB*TRPNCaH)+((1.-FrSBXB)*TRPNCaL)))));  
  double dTropTot = (Trop_conc*(((-dSOVFThin)*TRPNCaL)+((1.-SOVFThin)*dTRPNCaL)
                        +(dSOVFThin*((FrSBXB*TRPNCaH)+((1.-FrSBXB)*TRPNCaL)))
                        +(SOVFThin*(((dFrSBXB*TRPNCaH)+(FrSBXB*dTRPNCaH)+
                                    ((1.-FrSBXB)*dTRPNCaL))-(dFrSBXB*TRPNCaL)))));
  
  double dforce = ((kxb*dSOVFThick*((xXBpostr*XBpostr)+(xXBprer*XBprer)))+(kxb*SOVFThick*((dxXBpostr*XBpostr)+(xXBpostr*dXBpostr)+(dxXBprer*XBprer)+(xXBprer*dXBprer))));  
  
  // Calcium handling
  double cc = (((TropTot/1000.0)-Ca_i_total)*K_buf_c);
  double bc = (((K_buf_c+Buf_c)+(TropTot/1000.))-Ca_i_total);
  //Ca_i = (((((Ca_i_total-K_buf_c)-Buf_c)-(TropTot/1000.))+pow((pow(((K_buf_c+Buf_c+(TropTot/1000.))-Ca_i_total),2.)-(4.0*((TropTot/1000.)-Ca_i_total)*K_buf_c)),0.5))/2.);
  Ca_i = (((-bc)+pow(((bc*bc)-(4.0*cc)),0.5))/2.0);
  Ca_i_bufc = ((Buf_c*Ca_i)/(Ca_i+K_buf_c));
  Ca_sr_bufsr = (1.0/(1.0+((Buf_sr*K_buf_sr)/pow((Ca_SR+K_buf_sr),2.0))));   

  // convert units of Ca_i for coupling to Rice model
  Cai = 1000.0*Ca_i;   
  
  double E_Ca = 0.5 * RTOF * log(Ca_o/Ca_i);
  ICaL = g_CaL*d*f*fCa*4.*((V*F*F)/(R*T))*((Ca_i*exp(2.*VFORT)-0.341*Ca_o)/(exp(2.*VFORT)-1.));
  IbCa = g_bca * (V - E_Ca);
  
  double alpha_fCa = 1.0/(1.0 + pow((Ca_i/3.2500e-04),8.0));
  double beta_fCa = 0.1/(1.0 + exp((Ca_i - 5.0e-04)/1.0e-04));
  double gama_fCa = 0.2/(1.0 + exp((Ca_i - 7.5e-04)/8.0e-04));    
  double tau_fCa = 2.0;
  double fCa_inf = (alpha_fCa + beta_fCa + gama_fCa + 0.23)/1.460;
  double d_fCa = (fCa_inf - fCa)/tau_fCa;   
  
  double g_inf = (Ca_i < 3.5e-04) ? (1.0/(1.0+pow((Ca_i/3.5e-04),6.0))) : (1.0/(1.0+pow((Ca_i/3.5e-04),16.0)));
  double d_g = (g_inf - g)/tau_g;
  
  INaCa = (K_NaCa*((exp(gamma_tt*VFORT)*pow(Na_i,3.)*Ca_o)-(exp((gamma_tt-1.)*VFORT)*pow(Na_o,3.0)*Ca_i*alpha)))/((pow(Km_Nai,3.)+pow(Na_o,3.))*(Km_Ca+Ca_o)*(1.+(K_sat*exp((gamma_tt-1.0)*VFORT))));  
  IpCa = (g_pCa*Ca_i)/(Ca_i+K_pCa);  
  Iup = (Vmax_up/(1.0+(pow(K_up,2.0)/pow(Ca_i,2.0))));      
  Ileak = V_leak * (Ca_SR-Ca_i);    
  // End of calcium handling

  // Coupling TNNP + Rice
  dTRPNCaL = (konT * Cai * (1.-TRPNCaL)) - (koffLT*TRPNCaL);
  dTRPNCaH = (konT * Cai * (1.-TRPNCaH)) - (koffHT*TRPNCaH);  
  
  // derivatives
  d_dt_V = -(IK1+Ito+IKr+IKs+ICaL+INaK+INa+INaL+IbNa +INaCa+IbCa+IpK+IpCa+Istim);  
  d_dt_Ca_SR = ((Ca_sr_bufsr * Vc)/Vsr) * (Iup - (Irel + Ileak));
  d_dt_Ca_i_total = (((Ileak - Iup) + Irel)-(((1.*((ICaL + IbCa + IpCa)-(2.*INaCa)))/(2.*Vc*F))*Cm));
  d_dt_Na_i = (((-1.)*(INa + INaL + IbNa + (3.*INaK) + (3.*INaCa))*Cm)/(Vc*F)); 
  d_dt_K_i = (((-1.)*((IK1 + Ito + IKr + IKs + IpK + Istim) - (2.0*INaK))*Cm)/(1.0*Vc*F));  
  d_dt_N = (kp_nT * P) - (kn_pT * N);
  d_dt_N_NoXB = (kp_nT * P_NoXB) - (kn_pT * N_NoXB);
  d_dt_P_NoXB = (kn_pT * N_NoXB) - (kp_nT * P_NoXB);
  d_dt_Xr1 = (xr1_inf - Xr1)/tau_xr1;
  d_dt_Xr2 = (xr2_inf - Xr2)/tau_xr2;
  d_dt_Xs = (xs_inf - Xs)/tau_xs;
  d_dt_m = (m_inf - m)/tau_m;
  d_dt_h = (h_inf - h)/tau_h;
  d_dt_j = (j_inf - j)/tau_j;
  d_dt_d = (d_inf - d)/tau_d;
  d_dt_f = (f_inf - f)/tau_f;
  d_dt_fCa = (((0.01 * d_fCa) > 0.0) && (V > (-60.0))) ? 0 : d_fCa; 
  d_dt_s = (s_inf - s)/tau_s;
  d_dt_r = (r_inf - r)/tau_r;
  d_dt_g = (((0.01 * d_g) > 0.0) && (V > (-60.0))) ? 0.0 : d_g;
  
  // ---------------------------------------------------------------------------
 
  // copy new
  values[0]  = d_dt_V;  
  //values[1]  = d_dt_Xr1;
  //values[2]  = d_dt_Xr2;
  //values[3]  = d_dt_Xs;
  //values[4]  = d_dt_m;
  //values[5]  = d_dt_h;
  //values[6]  = d_dt_j;
  //values[7]  = d_dt_d;
  //values[8]  = d_dt_f;
  //values[9]  = d_dt_fCa;
  //values[10] = d_dt_s;
  //values[11] = d_dt_r;
  // values[12] = d_dt_g;
  // Rush-Larsen
  const double dt = 0.05;
  values[ 1] = xr1_inf + (Xr1 - xr1_inf) * exp(-dt / tau_xr1);
  values[ 2] = xr2_inf + (Xr2 - xr2_inf) * exp(-dt / tau_xr2);
  values[ 3] = xs_inf + (Xs - xs_inf) * exp(-dt / tau_xs);
  values[ 4] = m_inf + (m - m_inf) * exp(-dt / tau_m);
  values[ 5] = h_inf + (h - h_inf) * exp(-dt / tau_h);
  values[ 6] = j_inf + (j - j_inf) * exp(-dt / tau_j);
  values[ 7] = d_inf + (d - d_inf) * exp(-dt / tau_d);
  values[ 8] = f_inf + (f - f_inf) * exp(-dt / tau_f);
  values[ 9] = d_dt_fCa; 
  values[10] = s_inf + (s - s_inf) * exp(-dt / tau_s);
  values[11] = r_inf + (r - r_inf) * exp(-dt / tau_r);
  values[12] = d_dt_g; 
  
  values[13] = d_dt_Ca_SR;
  values[14] = d_dt_Ca_i_total;
  values[15] = d_dt_Na_i;
  values[16] = d_dt_K_i;
  values[17] = dTRPNCaL;
  values[18] = dTRPNCaH;
  values[19] = d_dt_N_NoXB;
  values[20] = d_dt_P_NoXB;
  values[21] = d_dt_N;
  values[22] = dXBprer;
  values[23] = dXBpostr;
  values[24] = dxXBprer;
  values[25] = dxXBpostr;
  values[26] = 0.0; //dSL;
  values[27] = dintf;
   
#ifdef HF
  //double d_dt_mL = alpha_mL * (1-mL) - beta_mL*mL;
  //double d_dt_hL = (hL_inf - hL)/tau_hL; 
  // uses the same m from INa (Beatriz Trevor, PLoS ONE, 2012)
  double d_dt_mL = m_inf  + (m  -  m_inf) * exp(-dt / tau_m);
  double d_dt_hL = hL_inf + (hL - hL_inf) * exp(-dt / tau_hL);
  values[28] = d_dt_mL;
  values[29] = d_dt_hL;
#endif
  
}

//inline double RiceTenTusscher::ifnumber_0()
//{
//  if(V < -40.)
//    return 0.057 * exp( -(V+80.)/6.8 );
//  else
//    return 0.0;
//}

//inline double RiceTenTusscher::ifnumber_1()
//{
//  if(V < -40.)
//    return ((2.7*exp((0.079*V))) + (3.1e+05 * exp(0.3485*V)));
//  else
//    return ((0.77/(0.13*(1.+exp(((V + 10.66)/(-11.1)))))));
//}

//inline double RiceTenTusscher::ifnumber_2()
//{
//  if(V < -40.)
//    return ( ((-2.5428e+04)*exp(0.2444*V)-
//	      (6.9480e-06*exp(-0.04391*V)) ) * (V+37.78) )/(1. + exp(0.311*(V+79.23)));
//  else
//    return 0.0;
//}

//inline double RiceTenTusscher::ifnumber_3()
//{
//  if(V < -40.)
//    return (((2.4240e-02*exp(((-1.0520e-02)*V)))/(1.+exp(((-0.13780)*(V+40.14))))));
//  else
//    return (0.6*exp(0.057*V))/(1.+exp((-0.1)*(V+32.)));
//}

//inline double RiceTenTusscher::ifnumber_4()
//{
//  if(((0.01 * d_fCa) > 0.0) && (V > (-60.0)))
//    return 0.0;
//  else
//    return d_fCa;
//}

//inline double RiceTenTusscher::ifnumber_5()
//{
//  if(Ca_i < 3.5e-04)
//    return (1.0/(1.0+pow((Ca_i/3.5e-04),6.0)));
//  else
//    return (1.0/(1.0+pow((Ca_i/3.5e-04),16.0)));
//}

//inline double RiceTenTusscher::ifnumber_6()
//{
//  if(((0.01 * d_g) > 0.0) && (V > (-60.0)))
//    return 0.0;
//  else
//    return d_g;
//}

//inline double RiceTenTusscher::ifnumber_8()
//{
//  if((len_thick/2.0) < (SL/2.))
//    return (len_thick/2.);
//  else
//    return (SL/2.);
//}
//
//inline double RiceTenTusscher::ifnumber_9()
//{
//  if(((SL/2.) - (SL-len_thin)) > (len_hbare/2.0))
//    return ((SL/2.) - (SL-len_thin));
//  else
//    return (len_hbare/2.);
//}

//inline double RiceTenTusscher::ifnumber_10()
//{
//  if((1./permtot) < 100.)
//    return (1./permtot);
//  else
//    return 100.;
//}

//inline double RiceTenTusscher::ifnumber_11(double sigmap, double sigman)
//{
//  if(xXBpostr < x_0)
//    return exp((sigmap*pow(((x_0-xXBpostr)/x_0),2.0)));
//  else 
//    return exp((sigman*pow(((xXBpostr-x_0)/x_0),2.0)));
//}

//inline double RiceTenTusscher::ifnumber_12()
//{
//  if(SL > SL_c)
//    return (PCon_c * (exp((PExp_c * fabs((SL-SL_c)))) - 1.0));
//  else
//    return 0.0;
//}

//inline double RiceTenTusscher::ifnumber_14()
//{
//  const double stretch=1.0;
//  if(!isometric && ((SL <= SLmax) && (SL > SLmin)) )
//    return SL*(stretch-1.0)/0.05;
//  //return (((intf + ((SLset - SL)*visc))/massf));
//  else
//    return 0.0;
//}

//inline double RiceTenTusscher::ifnumber_13()
//{
//  if(SEon == 1.0)
//    return (KSE * (SLset - SL));
//  else
//    return 0.0;
//}

//inline double RiceTenTusscher::ifnumber_15()
//{
//  if(SL < len_thick)
//    return -0.5*dSL;
//  else
//    return 0.0;
//}
//
//inline double RiceTenTusscher::ifnumber_16()
//{
//  if(((2.*len_thin)-SL) > len_hbare)
//    return -0.5*dSL;
//  else
//    return 0.0;
//}