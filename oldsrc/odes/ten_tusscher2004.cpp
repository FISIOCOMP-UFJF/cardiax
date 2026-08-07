#include "ten_tusscher2004.hpp"

#define SAC

TenTusscher2004::TenTusscher2004() : CellModel(19)
{
  const string vn[] = {"V", "Xr1", "Xr2", "Xs", "m", "h", "j", "d", "f",
                       "f2", "fCass", "s", "r", "Ca_i", "Ca_SR", "Ca_ss",
                       "R_prime", "Na_i", "K_i"};

  for(int i=0; i<19; i++)
    var_names.insert( std::pair<int, std::string>(i,  vn[i]) );

  // variables for RL integration in Explicit Euler
  for(int i=1; i<=12; i++) rlvars.insert(i);

}

void TenTusscher2004::init(double * values) const
{
  assert(values != 0);

  values[ 0] = -85.23;   // V;       millivolt
  values[ 1] = 0.00621;  // Xr1;     dimensionless
  values[ 2] = 0.4712;   // Xr2;     dimensionless
  values[ 3] = 0.0095;   // Xs;      dimensionless
  values[ 4] = 0.00172;  // m;       dimensionless
  values[ 5] = 0.7444;   // h;       dimensionless
  values[ 6] = 0.7045;   // j;       dimensionless
  values[ 7] = 3.373e-5; // d;       dimensionless
  values[ 8] = 0.7888;   // f;       dimensionless
  values[ 9] = 0.9755;   // f2;      dimensionless
  values[10] = 0.9953;   // fCass;   dimensionless
  values[11] = 0.999998; // s;       dimensionless
  values[12] = 2.42e-8;  // r;       dimensionless
  values[13] = 0.000126; // Ca_i;    millimolar
  values[14] = 3.64;     // Ca_SR;   millimolar
  values[15] = 0.00036;  // Ca_ss;   millimolar
  values[16] = 0.9073;   // R_prime; dimensionless
  values[17] = 8.604;    // Na_i;    millimolar
  values[18] = 136.89;   // K_i;     millimolar
}

void TenTusscher2004::equation(const double time, 
                               const double * statevars, 
                               double * values)
{
  // State variables
  const double V = statevars[0];        // Membrane variable
  const double Xr1 = statevars[1];      // Rapid time dependent potassium current Xr1
  const double Xr2 = statevars[2];      // Rapid time dependent potassium current Xr2
  const double Xs = statevars[3];       // Slow time dependent potassium current Xs
  const double m = statevars[4];        // Fast sodium current m
  const double h = statevars[5];        // Fast sodium current h gate
  const double j = statevars[6];        // Fast sodium current j gate
  const double d = statevars[7];        // L type Ca current d gate
  const double f = statevars[8];        // var_L_type_Ca_current_f_gate__f
  const double f2 = statevars[9];       // var_L_type_Ca_current_f2_gate__f2
  const double fCass = statevars[10];   // L_type_Ca_current__fCass
  const double s = statevars[11];       // var_transient_outward_current_s_gate__s
  const double r = statevars[12];       // _transient_outward_current__r
  const double Ca_i = statevars[13];    // var_calcium_dynamics__Ca_i
  const double Ca_SR = statevars[14];   // calcium sarcoplasmic reticulum
  const double Ca_ss = statevars[15];   // Ca_ss
  const double R_prime = statevars[16]; // R_prime
  const double Na_i = statevars[17];    // var_sodium_dynamics__Na_i
  const double K_i = statevars[18];     // var_potassium_dynamics__K_i

  // Some constants
  const double R    = 8314.472;
  const double T    = 310.0;
  const double F    = 96485.3415;
  const double Cm   = 0.185;
  const double V_c  = 0.016404;

  // Precalculations
  const double RToF  = (R*T)/F;
  const double VcF   = V_c * F;
  const double VFoRT = (V*F)/(R*T);

  // GAMBIARA RL
  double dt = 0.01;

  // Calculations            
  const double Ko  = 5.4;
  const double Nao = 140.0;
  const double Cao = 2.0;    

  const double P_kna  = 0.03;         
  const double Km_Nai = 87.5;
  const double K_mk = 1.0;
  const double P_NaK = 2.724;
  const double K_mNa = 40.0;
  const double V_sr = 0.001094;
  const double V_ss = 5.468e-05;
  const double V_rel = 0.102;
  const double k1_prime = 0.15;
  const double max_sr = 2.5;
  const double min_sr = 1.0;
  const double EC = 1.5;        
  const double k3 = 0.06;    
  const double k2_prime2 = 0.045;
  const double k4 = 0.005;
  const double Buf_c = 0.2;
  const double K_buf_c = 0.001;
  const double K_buf_sr = 0.3;
  const double Buf_sr = 10.0;
  const double Buf_ss = 0.4;
  const double K_buf_ss = 0.00025;    
  const double Vmax_up = 0.006375;
  const double K_up = 0.00025;    
  const double V_leak = 0.00036;   
  const double V_xfer = 0.0038;

  // reversal potentials
  double EK  = RToF * log(Ko/K_i);
  double EKs = RToF * log((Ko + (P_kna * Nao)) / (K_i + (P_kna * Na_i)));
  double ENa = RToF * log(Nao / Na_i);
  double ECa = 0.5 * RToF * log(Cao / Ca_i);

  // inward rectifier potassium current
  double beta_K1  = ((3.0*exp(0.0002*((V-EK)+100.0)))+exp(0.1*((V-EK)-10.0)))/(1.0+exp((-0.5)*(V-EK)));
  double alpha_K1 = 0.1 / (1.0 + exp(0.06 * ((V - EK) - 200.0)));
  double xK1_inf  = alpha_K1 / (alpha_K1 + beta_K1);
  double g_K1     = 5.405;
  double IK1 = g_K1 * xK1_inf * (V - EK);

  // transient outward current  
  double g_to = 0.294;
  if (type != ENDO) 
    g_to = 0.294;
  else
    g_to = 0.073;        
  double Ito = g_to * r * s * (V - EK);

  // rapid time dependent potassium current
  const double g_Kr = 0.153;
  double IKr = g_Kr * Xr1 * Xr2 * (V - EK) * sqrt(Ko / 5.4);

  // slow time depedent potassium current
  double g_Ks = 0.098;
  if (type == EPI)        g_Ks = 0.245;
  else if (type == ENDO)  g_Ks = 0.245;
  else if (type == MCELL) g_Ks = 0.062;   
  double IKs = g_Ks * pow(Xs, 2) * (V - EKs);   

  // L-type calcium current
  //double var_L_type_Ca_current__Ca_ss = Ca_ss;
  const double g_CaL = 3.98e-05;
  double ICaL = (V < 15.0 - 1.0e-5 || V > 15.0 + 1.0e-5) 
      ? ((((g_CaL*d*f*f2*fCass*4.0*(V-15.0)*pow(F,2))/(R*T))*((0.25*Ca_ss*exp((2.0*(V-15.0)*F)/(R*T)))-Cao))/(exp((2.0*(V-15.0)*F)/(R*T))-1.0)) 
      : g_CaL*d*f*f2*fCass*2.0*F*(0.25*Ca_ss-Cao);

  // sodium potassium pump current
  double INaK = ((((P_NaK*Ko)/(Ko+K_mk))*Na_i)/(Na_i+K_mNa))/(1.0+(0.1245*exp(((-0.1)*V*F)/(R*T)))+(0.0353*exp(((-V)*F)/(R*T))));

  // fast sodium current
  const double g_Na = 14.838;
  double INa = g_Na * pow(m, 3) * h * j * (V - ENa);

  // background sodium current
  const double g_bna = 0.00029;
  double IbNa = g_bna * (V - ENa);

  // sodium calcium exchanger
  const double alpha = 2.5;
  const double gamma = 0.35;
  const double K_sat = 0.1;
  const double Km_Ca = 1.38;    
  const double K_NaCa = 1000.0;
  //double INaCa = (K_NaCa*((exp((gamma*V*F)/(R*T))*pow(Na_i,3)*Cao)-(exp(((gamma-1.0)*V*F)/(R*T))*pow(Nao,3)*Ca_i*alpha)))/((pow(Km_Nai,3)+pow(Nao,3))*(Km_Ca+Cao)*(1.0+(K_sat*exp( ((gamma-1.0)*V*F)/(R*T) ))));
  double INaCa = (K_NaCa*((exp(gamma*VFoRT)*pow(Na_i,3)*Cao)-(exp((gamma-1.0)*VFoRT)*pow(Nao,3)*Ca_i*alpha)))/((pow(Km_Nai,3)+pow(Nao,3))*(Km_Ca+Cao)*(1.0+(K_sat*exp((gamma-1.0)*VFoRT))));    

  // calcium background current
  const double g_bca = 0.000592;
  double IbCa = g_bca * (V - ECa);

  // potassium pump current
  const double g_pK = 0.0146;
  double IpK = (g_pK*(V-EK))/(1.0+exp((25.0-V)/5.98));

  // calcium pump current
  const double K_pCa = 0.0005;
  const double g_pCa = 0.1238;
  double IpCa = (g_pCa*Ca_i)/(Ca_i+K_pCa);

  // stimulus
//    double Istim;
//    double amp    = -52.0;
//    double start  = 50.0;
//    double end    = 50000.0;
//    double period = 1000.0;
//    double dur    = 1.0;
//    
//    if( (time>=start) && (time<=end) && 
//        (((time-start)-(floor(((time-start)/period))*period))<=dur) )
//      Istim = amp;
//    else
//      Istim = 0.0;

  double Istim = i_stim;

  // Gating variables
  double xr1_inf   = 1.0 / (1.0 + exp(((-26.0) - V) / 7.0));
  double alpha_xr1 = 450.0 / (1.0 + exp(((-45.0) - V) / 10.0));
  double beta_xr1  = 6.0 / (1.0 + exp((V + 30.0) / 11.5));
  double tau_xr1   = 1.0 * alpha_xr1 * beta_xr1;

  double xr2_inf   = 1.0 / (1.0 + exp((V + 88.0) / 24.0));
  double alpha_xr2 = 3.0 / (1.0 + exp(((-60.0) - V) / 20.0));
  double beta_xr2  = 1.12 / (1.0 + exp((V - 60.0) / 20.0));
  double tau_xr2   = 1.0 * alpha_xr2 * beta_xr2;

  double xs_inf   = 1.0 / (1.0 + exp(((-5.0) - V) / 14.0));
  double alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0));
  double beta_xs  = 1.0 / (1.0 + exp((V - 35.0) / 15.0));
  double tau_xs   = (1.0 * alpha_xs * beta_xs) + 80.0;

  double m_inf   = 1.0 / pow(1.0 + exp(((-56.86) - V) / 9.03), 2);
  double alpha_m = 1.0 / (1.0 + exp(((-60.0) - V) / 5.0));
  double beta_m  = (0.1/(1.0+exp((V+35.0)/5.)))+(0.1/(1.0+exp((V-50.0)/200.)));
  double tau_m   = 1.0 * alpha_m * beta_m;

  double h_inf   = 1.0 / pow(1.0 + exp((V + 71.55) / 7.43), 2);
  double alpha_h = (V < (-40.0)) ? (0.057 * exp((-(V + 80.0)) / 6.8)) : 0.0;
  double beta_h  = (V < (-40.0)) ? ((2.7 * exp(0.079 * V)) + (310000.0 * exp(0.3485 * V))) : (0.77 / (0.13 * (1.0 + exp((V + 10.66) / (-11.1)))));
  double tau_h   = 1.0 / (alpha_h + beta_h);

  double j_inf   = 1.0 / pow(1.0 + exp((V + 71.55) / 7.43), 2);
  double alpha_j = (V < (-40.0)) ? ((((((-25428.0) * exp(0.2444 * V)) - (6.948e-06 * exp((-0.04391) * V))) * (V + 37.78)) / 1.0) / (1.0 + exp(0.311 * (V + 79.23)))) : 0.0;
  double beta_j  = (V < (-40.0)) ? ((0.02424 * exp((-0.01052) * V)) / (1.0 + exp((-0.1378) * (V + 40.14)))) : ((0.6 * exp(0.057 * V)) / (1.0 + exp((-0.1) * (V + 32.0))));
  double tau_j   = 1.0 / (alpha_j + beta_j);

  double d_inf   = 1.0 / (1.0 + exp(((-8.0) - V) / 7.5));
  double alpha_d = (1.4 / (1.0 + exp(((-35.0) - V) / 13.0))) + 0.25;
  double beta_d  = 1.4 / (1.0 + exp((V + 5.0) / 5.0));
  double gamma_d = 1.0 / (1.0 + exp((50.0 - V) / 20.0));
  double tau_d   = (1.0 * alpha_d * beta_d) + gamma_d;

  double f_inf = 1.0 / (1.0 + exp((V + 20.0) / 7.0));
  double tau_f = (1102.5 * exp((-pow(V + 27.0, 2)) / 225.0)) + (200.0 / (1.0 + exp((13.0 - V) / 10.0))) + (180.0 / (1.0 + exp((V + 30.0) / 10.0))) + 20.0;

  double f2_inf = (0.67/(1.0+exp((V+35.0)/7.0)))+0.33;
  double tau_f2 = (562.0*exp((-pow(V+27.0,2))/240.0))+(31.0/(1.0+exp((25.0-V)/10.0)))+(80.0/(1.0+exp((V+30.0)/10.)));

  double fCass_inf = (0.6/(1.0+pow(Ca_ss/0.05,2)))+0.4;
  double tau_fCass = (80.0/(1.0+pow(Ca_ss/0.05,2)))+2.0;

  double s_inf = 1.0 / (1.0 + exp((V+20.0)/5.0));
  double tau_s = (85.0*exp((-pow(V+45.0,2))/320.))+(5./(1.+exp((V-20.0)/5.0)))+3.0;

  double r_inf = 1.0 / (1.0+exp((20.0-V)/6.0));
  double tau_r = (9.5*exp((-pow(V+40.0,2))/1800.0))+0.8;

  double kcasr = max_sr - ((max_sr - min_sr) / (1.0 + pow(EC / Ca_SR, 2)));
  double k1 = k1_prime / kcasr;

  double O = (k1 * pow(Ca_ss, 2) * R_prime) / (k3 + (k1 * pow(Ca_ss, 2)));
  double Irel = V_rel * O * (Ca_SR - Ca_ss);    

  double Iup = Vmax_up / (1.0 + (pow(K_up, 2) / pow(Ca_i, 2)));          
  double Ileak = V_leak * (Ca_SR - Ca_i);   
  double Ixfer = V_xfer * (Ca_ss - Ca_i);    

  double k2 = k2_prime2 * kcasr;   
  double Ca_i_bufc   = 1.0 / (1.0 + ((Buf_c * K_buf_c) / pow(Ca_i + K_buf_c, 2)));    
  double Ca_sr_bufsr = 1.0 / (1.0 + ((Buf_sr * K_buf_sr) / pow(Ca_SR + K_buf_sr, 2)));    
  double Ca_ss_bufss = 1.0 / (1.0 + ((Buf_ss * K_buf_ss) / pow(Ca_ss + K_buf_ss, 2)));

  //double d_dt_rapid_time_dependent_potassium_current_Xr1_gate__Xr1 = (xr1_inf - Xr1) / tau_xr1;
  //double d_dt_rapid_time_dependent_potassium_current_Xr2_gate__Xr2 = (xr2_inf - Xr2) / tau_xr2;
  //double d_dt_slow_time_dependent_potassium_current_Xs_gate__Xs = (xs_inf - Xs) / tau_xs;
  //double d_dt_fast_sodium_current_m_gate__m = (m_inf - m) / tau_m;
  //double d_dt_fast_sodium_current_h_gate__h = (h_inf - h) / tau_h;
  //double d_dt_fast_sodium_current_j_gate__j = (j_inf - j) / tau_j;
  //double d_dt_L_type_Ca_current_d_gate__d = (d_inf - d) / tau_d;
  //double d_dt_L_type_Ca_current_f_gate__f = (f_inf - f) / tau_f;
  //double d_dt_L_type_Ca_current_f2_gate__f2 = (f2_inf - f2) / tau_f2;
  //double d_dt_L_type_Ca_current_fCass_gate__fCass = (fCass_inf - fCass) / tau_fCass;
  //double d_dt_transient_outward_current_s_gate__s = (s_inf - s) / tau_s;
  //double d_dt_transient_outward_current_r_gate__r = (r_inf - r) / tau_r;

  double d_dt_V = -(IK1+Ito+IKr+IKs+ICaL+INaK+INa+IbNa+INaCa+IbCa+IpK+IpCa+Istim);
  double d_dt_R_prime = ((-k2)*Ca_ss*R_prime) + (k4*(1.0-R_prime));
  double d_dt_Ca_i  = Ca_i_bufc*(((((Ileak-Iup)*V_sr)/V_c)+Ixfer)-((((IbCa+IpCa)-(2*INaCa))*Cm)/(2*VcF)));
  double d_dt_Ca_SR = Ca_sr_bufsr*(Iup-(Irel+Ileak));
  double d_dt_Ca_ss = Ca_ss_bufss*((((-ICaL*Cm)/(2*V_ss*F))+((Irel*V_sr)/V_ss))-((Ixfer*V_c)/V_ss));
  double d_dt_Na_i = ((-(INa+IbNa+(3*INaK)+(3*INaCa)))/(VcF))*Cm;
  double d_dt_K_i  = ((-((IK1+Ito+IKr+IKs+IpK+Istim)-(2*INaK)))/(VcF))*Cm;

  // Isac
# ifdef SAC
  double Gs, Es, lambda, lambda_max, Isac;
  Gs =         0.075;
  Es =         -50;
  lambda =     1.05;
  lambda_max = 1.1;
  if(lambda > 1.0)
    Isac = Gs*(lambda -1.0)*(V - Es);
    //Isac = Gs*((lambda -1.0)/(lambda_max-1.0))*(V - Es);
  else
    Isac = 0.0;
  d_dt_V -= Isac;
# endif

  values[ 0] = d_dt_V;

  values[ 1] = xr1_inf + (Xr1 - xr1_inf) * exp(-dt / tau_xr1);
  values[ 2] = xr2_inf + (Xr2 - xr2_inf) * exp(-dt / tau_xr2);
  values[ 3] = xs_inf + (Xs - xs_inf) * exp(-dt / tau_xs);
  values[ 4] = m_inf + (m - m_inf) * exp(-dt / tau_m);
  values[ 5] = h_inf + (h - h_inf) * exp(-dt / tau_h);
  values[ 6] = j_inf + (j - j_inf) * exp(-dt / tau_j);
  values[ 7] = d_inf + (d - d_inf) * exp(-dt / tau_d);
  values[ 8] = f_inf + (f - f_inf) * exp(-dt / tau_f);
  values[ 9] = f2_inf + (f2 - f2_inf) * exp(-dt / tau_f2);
  values[10] = fCass_inf + (fCass - fCass_inf) * exp(-dt / tau_fCass);
  values[11] = s_inf + (s - s_inf) * exp(-dt / tau_s);
  values[12] = r_inf + (r - r_inf) * exp(-dt / tau_r);

  values[13] = d_dt_Ca_i;
  values[14] = d_dt_Ca_SR;
  values[15] = d_dt_Ca_ss;
  values[16] = d_dt_R_prime;
  values[17] = d_dt_Na_i;
  values[18] = d_dt_K_i;
}
