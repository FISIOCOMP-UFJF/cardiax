#include "ten_tusscher_ta.hpp"

// Heaviside function
inline double MYMAX(double x, double y) { return (x < y) ? y : x; }

TenTusscherTa::TenTusscherTa() : CellModel(21), marked(false), elapsedTime(0)
{
  // default type of cell
  type = EPI;

  // state variables names
  const string vn[] = {"V", "Xr1", "Xr2", "Xs", "m", "h", "j", "d", "f",
                       "f2", "fCass", "s", "r", "Ca_i", "Ca_SR", "Ca_ss",
                       "R_prime", "Na_i", "K_i", "Ta", "Aux"};

  for(int i=0; i<20; i++)
    var_names.insert( std::pair<int, std::string>(i,  vn[i]) );

  // variables for RL integration in Explicit Euler
  for(int i=1; i<=12; i++) rlvars.insert(i);
  
  // some test for active tension (Ta)
  rlvars.insert(19); 
  rlvars.insert(20); 
}

void TenTusscherTa::init(double * values) const
{
    /*
 // assert(values != 0);
  values[ 0] = -85.423;  // V; millivolt; -85.423
  values[ 1] = 0.0165;   // Xr1; dimensionless; 0.0165
  values[ 2] = 0.473;    // Xr2; dimensionless; 0.473
  values[ 3] = 0.0174;   // Xs; dimensionless; 0.0174
  values[ 4] = 0.00165;  // m; dimensionless; 0.00165
  values[ 5] = 0.749;    // h; dimensionless; 0.749
  values[ 6] = 0.6788;   // j; dimensionless; 0.6788
  values[ 7] = 3.288e-5; // d; dimensionless; 3.288e-5
  values[ 8] = 0.7026;   // f; dimensionless; 0.7026
  values[ 9] = 0.9526;   // f2; dimensionless; 0.9526
  values[10] = 0.9942;   // fCass; dimensionless; 0.9942
  values[11] = 0.999998; // s; dimensionless; 0.999998
  values[12] = 2.347e-8; // r; dimensionless; 2.347e-8
  values[13] = 0.000153; // Ca_i; millimolar; 0.000153
  values[14] = 4.272;    // Ca_SR; millimolar; 4.272
  values[15] = 0.00042;  // Ca_ss; millimolar; 0.00042
  values[16] = 0.8978;   // R_prime; dimensionless; 0.8978
  values[17] = 10.132;   // Na_i; millimolar; 10.132
  values[18] = 138.52;   // K_i; millimolar; 138.52
  
  values[19] = 0.0;      // Ta; kPa
  values[20] = 0.0;      // AUX
  */

  //
  // Electrophysioloy benchmark initial conditions
  //  
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
  values[19] = 1.0;      // NOT USED
  
  values[19] = 0.0;      // Ta; kPa
  values[20] = 0.0;      // AUX
}

void TenTusscherTa::equation(const double time, const double * statevars, double * values)
{
  // State variables
  const double V = statevars[0];   // Transmembrane potential
  const double Xr1 = statevars[1]; // Rapid time dependent potassium Xr1 gate
  const double Xr2 = statevars[2]; // Rapid time dependent potassium Xr2 gate
  const double Xs = statevars[3];  // Slow time dependent potassium Xs gate
  const double m = statevars[4];   // Fast sodium current m gate
  const double h = statevars[5];   // Fast sodium current h gate
  const double j = statevars[6];   // Fast sodium current j gate
  const double d = statevars[7];   // L-type Ca current d gate
  const double f = statevars[8];   // L-type Ca current f gate
  const double f2 = statevars[9];  // L-type Ca current f2 gate
  const double fCass = statevars[10];   // L-type Ca current fCass gate
  const double s = statevars[11];       // Transient outward current s gate
  const double r = statevars[12];       // Transient outward current r gate
  const double Ca_i = statevars[13];    // Calcium dynamics
  const double Ca_SR = statevars[14];   // Calcium dynamics
  const double Ca_ss = statevars[15];   // Calcium dynamics Ca_ss
  const double R_prime = statevars[16]; // Calcium dynamics 
  const double Na_i = statevars[17];    // Sodium dynamics
  const double K_i = statevars[18];     // Potassium dynamics
  const double Ta = statevars[19];      // Active tension
  const double aux2 = statevars[20];

  // Some constants
  const double R = 8314.472;
  const double T = 310.0;
  const double F = 96485.3415;
  const double RTONF = (R*T)/F;
  const double Cm = 0.185;
  
  // Intracellular volumes
  const double Vc = 0.016404;
  const double Vsr = 0.001094;

  // External concentrations
  const double Ko = 5.4;
  const double Cao = 2.0;
  const double Nao = 140.0;

  // Parameters for currents
  double GK1 = 5.405;
  const double Gkr = 0.096;

  const double GCaL = 0.000175;
  
  const double GNa = 14.838;
  const double Gbna = 0.00029;
  const double GpK = 0.0146;
  const double Gbca = 0.000592;
  const double KmK = 1.0;
  const double P_NaK = 2.724;
  const double P_kna = 0.03;  
  const double K_mNa = 40.0;

  // Parameters for INaCa
  const double ksat = 0.1;
  const double KmCa = 1.38;
  double K_NaCa = 1000.0;
  const double KmNai = 87.5;

  // Parameters for IpCa
  const double GpCa = 0.1238;
  const double KpCa = 0.0005;

  // Parameters for Calcium dynamics
  const double EC = 1.5;
  const double V_rel = 0.102;
  const double k1_prime = 0.15;
  const double max_sr = 2.5;
  const double min_sr = 1.0;
  const double k3 = 0.06;
  double Vmax_up = 0.006375;
  const double K_up = 0.00025;  
  double V_leak = 0.00036;
  const double V_xfer = 0.0038;
  const double k2_prime = 0.045;
  const double k4 = 0.005;
  const double Buf_c = 0.2;
  const double K_buf_c = 0.001;
  const double K_buf_sr = 0.3;
  const double Buf_sr = 10.0;
  const double Buf_ss = 0.4;
  const double K_buf_ss = 0.00025;
  const double V_ss = 5.468e-05;   

  // Cell type specific conductances (Endo, M, Epi, Apex or Base)
  double Gto = 0.294;
  if (type != ENDO)
    Gto = 0.294;
  else
    Gto = 0.073;

 double Gks = 0.245;
  if (type == APEX)
    Gks = 0.36064;
  else if (type == BASE)
    Gks = 0.54096;
  else if (type == EPI)
    Gks = 0.245;
  else if (type == ENDO)
    Gks = 0.245;
  else if (type == MCELL)
    Gks = 0.062;

  /*
  if (type == EPI)
    Gks = 0.245 * 0.4;
  else if (type == ENDO)
    Gks = 0.245 * 0.4;
  else if (type == MCELL)
    Gks = 0.062 * 0.4;
      
  if (type != ENDO)
    Gto = 0.294 * 0.65;
  else
    Gto = 0.073 * 0.65;
  
  GK1     = 5.405 * 0.5;
  V_leak  = 1.32 * 0.00036; 
  Vmax_up = 0.006375 * 0.6014;  
  K_NaCa  = 2000;
      
  */
  // Needed to compute currents
  double E_K  = RTONF * log(Ko/K_i);
  double E_Ks = RTONF * log((Ko + (P_kna * Nao)) / (K_i + (P_kna * Na_i)));
  double E_Na = RTONF * log(Nao / Na_i);
  double E_Ca = (0.5 * RTONF) * log(Cao / Ca_i);

  // Inward rectifier potassium current 
  double beta_K1 = ((3.0 * exp(0.0002 * ((V - E_K) + 100.0))) + exp(0.1 * ((V - E_K) - 10.0))) / (1.0 + exp((-0.5) * (V - E_K)));
  double alpha_K1 = 0.1 / (1.0 + exp(0.06 * ((V - E_K) - 200.0)));
  double xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
  double IK1 = GK1 * xK1_inf * (V - E_K);

  // Compute currents
  double Ito = Gto * r * s * (V - E_K);
  double IKr = Gkr * Xr1 * Xr2 * (V - E_K) * sqrt(Ko / 5.4);
  double IKs = Gks * pow(Xs, 2) * (V - E_Ks);
  double ICaL = (V < 15.0-1.0e-5 || V > 15.0+1.0e-5) ? ((((GCaL*d*f*f2*fCass*4.0*(V-15.0)*pow(F,2)) / (R * T)) * ((0.25 * Ca_ss * exp((2.0 * (V - 15.0) * F) / (R*T))) - Cao)) / (exp((2.0 * (V-15.0) * F) / (R * T)) - 1.0)) : GCaL * d * f * f2 * fCass * 2.0 * F * (0.25 * Ca_ss  - Cao);
  double INaK = ((((P_NaK * Ko) / (Ko + KmK)) * Na_i) / (Na_i + K_mNa)) / (1.0 + (0.1245 * exp(((-0.1) * V * F) / (R * T))) + (0.0353 * exp(((-V) * F) / (R * T))));
  double INa = GNa * pow(m, 3) * h * j * (V - E_Na);
  double IbNa = Gbna * (V - E_Na);
  double IbCa = Gbca * (V - E_Ca);
  double IpK = (GpK * (V - E_K)) / (1.0 + exp((25.0 - V) / 5.98));
  double IpCa = (GpCa * Ca_i) / (Ca_i + KpCa);
  
  double Istim = i_stim; 
  
  //double Istim = /*i_stim*/stim(time);
//   if(Istim!=0.0)printf(":::::::%f\n", Istim);


  //double i_stimAmplitude = -5.200000000000000e+01;
  //double i_stimStart = 0.000000000000000e+00;
  //double i_stimEnd = 5.000000000000000e+04;
  //double i_stimPeriod = 1.000000000000000e+03;
  //double i_stimPulseDuration = 1.000000000000000e+00;
  //if(((time>=i_stimStart)&&(time<=i_stimEnd)&&(((time-i_stimStart)-(floor(((time-i_stimStart)/i_stimPeriod))*i_stimPeriod))<=i_stimPulseDuration)))
  //  Istim = i_stimAmplitude;
  //else
  //  Istim=0;
    
  // Sodium Calcium exchanger current
  const double alpha = 2.5;
  const double gamma = 0.35;
  double INaCa = (K_NaCa * ((exp((gamma * V * F) / (R * T)) * pow(Na_i, 3) * Cao) - (exp(((gamma - 1.0) * V * F) / (R * T)) * pow(Nao, 3) * Ca_i * alpha))) / ((pow(KmNai, 3) + pow(Nao, 3)) * (KmCa + Cao) * (1.0 + (ksat * exp(((gamma - 1.0) * V * F) / (R * T)))));

  // Rapid time dependent potassium current (Xr1 gate)
  double xr1_inf = 1.0 / (1.0 + exp(((-26.0) - V) / 7.0));
  double alpha_xr1 = 450.0 / (1.0 + exp(((-45.0) - V) / 10.0));
  double beta_xr1 = 6.0 / (1.0 + exp((V + 30.0) / 11.5));
  double tau_xr1 = 1.0 * alpha_xr1 * beta_xr1;

  // Rapid time dependent potassium current (Xr2 gate)
  double xr2_inf = 1.0 / (1.0 + exp((V + 88.0) / 24.0));
  double alpha_xr2 = 3.0 / (1.0 + exp(((-60.0) - V) / 20.0));
  double beta_xr2 = 1.12 / (1.0 + exp((V - 60.0) / 20.0));
  double tau_xr2 = 1.0 * alpha_xr2 * beta_xr2;

  // Slow time dependent potassium current (Xs gate)
  double xs_inf = 1.0 / (1.0 + exp(((-5.0) - V) / 14.0));
  double alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0));
  double beta_xs = 1.0 / (1.0 + exp((V - 35.0) / 15.0));
  double tau_xs = (1.0 * alpha_xs * beta_xs) + 80.0;

  // Fast sodium current (m gate)
  double m_inf = 1.0 / pow(1.0 + exp(((-56.86) - V) / 9.03), 2);
  double alpha_m = 1.0 / (1.0 + exp(((-60.0) - V) / 5.0));
  double beta_m = (0.1 / (1.0 + exp((V + 35.0) / 5.0))) + (0.1 / (1.0 + exp((V - 50.0) / 200.0)));
  double tau_m = 1.0 * alpha_m * beta_m;

  // Fast sodium current (h gate)
  double h_inf = 1.0 / pow(1.0 + exp((V + 71.55) / 7.43), 2);
  double alpha_h = (V < (-40.0)) ? (0.057 * exp((-(V + 80.0)) / 6.8)) : 0.0;
  double beta_h = (V < (-40.0)) ? ((2.7 * exp(0.079 * V)) + (310000.0 * exp(0.3485 * V))) : (0.77 / (0.13 * (1.0 + exp((V + 10.66) / (-11.1)))));
  double tau_h = 1.0 / (alpha_h + beta_h);

  // Fast sodium current (j gate)
  double j_inf = 1.0 / pow(1.0 + exp((V + 71.55) / 7.43), 2);
  double alpha_j = (V < (-40.0)) ? ((((((-25428.0) * exp(0.2444 * V)) - (6.948e-06 * exp((-0.04391) * V))) * (V + 37.78)) / 1.0) / (1.0 + exp(0.311 * (V + 79.23)))) : 0.0;
  double beta_j = (V < (-40.0)) ? ((0.02424 * exp((-0.01052) * V)) / (1.0 + exp((-0.1378) * (V + 40.14)))) : ((0.6 * exp(0.057 * V)) / (1.0 + exp((-0.1) * (V + 32.0))));
  double tau_j = 1.0 / (alpha_j + beta_j);

  // L-type Calcium current (d gate)
  double d_inf = 1.0 / (1.0 + exp(((-8.0) - V) / 7.5));
  double alpha_d = (1.4 / (1.0 + exp(((-35.0) - V) / 13.0))) + 0.25;
  double beta_d = 1.4 / (1.0 + exp((V + 5.0) / 5.0));
  double gamma_d = 1.0 / (1.0 + exp((50.0 - V) / 20.0));
  double tau_d = (1.0 * alpha_d * beta_d) + gamma_d;

  // L-type Calcium current (f, f2 and fCass gates)
  double f_inf = 1.0 / (1.0 + exp((V + 20.0) / 7.0));
  double tau_f = (1102.5 * exp((-pow(V + 27.0, 2)) / 225.0)) + (200.0 / (1.0 + exp((13.0 - V) / 10.0))) + (180.0 / (1.0 + exp((V + 30.0) / 10.0))) + 20.0;
  double f2_inf = (0.67 / (1.0 + exp((V + 35.0) / 7.0))) + 0.33;
  double tau_f2 = (562.0 * exp((-pow(V + 27.0, 2)) / 240.0)) + (31.0 / (1.0 + exp((25.0 - V) / 10.0))) + (80.0 / (1.0 + exp((V + 30.0) / 10.0)));
  double fCass_inf = (0.6 / (1.0 + pow(Ca_ss / 0.05, 2))) + 0.4;
  double tau_fCass = (80.0 / (1.0 + pow(Ca_ss / 0.05, 2))) + 2.0;

  // Transient outward current r gate
  double r_inf = 1.0 / (1.0 + exp((20.0 - V) / 6.0));
  double tau_r = (9.5 * exp((-pow(V + 40.0, 2)) / 1800.0)) + 0.8;
    
  // Transient outward current s gate (endo, M, epi)
  double s_inf;
  double tau_s;
  if (type != ENDO)
  {
    s_inf = 1.0 / (1.0 + exp((V + 20.0) / 5.0));
    tau_s = (85.0 * exp((-pow(V + 45.0, 2)) / 320.0)) + (5.0 / (1.0 + exp((V - 20.0) / 5.0))) + 3.0;  
  }
  else // ENDO (or BASE or APEX)
  {
    s_inf = 1.0 / (1.0 + exp((V + 28.0) / 5.0));
    tau_s = (1000.0 * exp((-pow(V + 67.0, 2)) / 1000.0)) + 8.0;
  }
  
  // Calcium dynamics
  double kcasr = max_sr - ((max_sr - min_sr) / (1.0 + pow(EC / Ca_SR, 2)));
  double k1 = k1_prime / kcasr;
  double var_calcium_dynamics__O = (k1 * pow(Ca_ss, 2) * R_prime) / (k3 + (k1 * pow(Ca_ss, 2)));
  double Irel = V_rel * var_calcium_dynamics__O * (Ca_SR - Ca_ss);
  double Iup = Vmax_up / (1.0 + (pow(K_up, 2) / pow(Ca_i, 2)));
  double Ileak = V_leak * (Ca_SR - Ca_i); 
  double Ixfer = V_xfer * (Ca_ss - Ca_i);
  double k2 = k2_prime * kcasr;
  double Ca_i_bufc = 1.0 / (1.0 + ((Buf_c * K_buf_c) / pow(Ca_i + K_buf_c, 2)));
  double Ca_sr_bufsr = 1.0 / (1.0 + ((Buf_sr * K_buf_sr) / pow(Ca_SR + K_buf_sr, 2)));
  double Ca_ss_bufss = 1.0 / (1.0 + ((Buf_ss * K_buf_ss) / pow(Ca_ss + K_buf_ss, 2)));
  
  // Rate calculations
  double d_dt_V = -(IK1 + Ito + IKr + IKs + ICaL + INaK + INa 
		  + IbNa + INaCa + IbCa + IpK + IpCa + Istim);

  /*
  double d_dt_Xr1 = (xr1_inf - Xr1) / tau_xr1;
  double d_dt_Xr2 = (xr2_inf - Xr2) / tau_xr2;
  double d_dt_Xs = (xs_inf - Xs) / tau_xs;
  double d_dt_m = (m_inf - m) / tau_m;
  double d_dt_h = (h_inf - h) / tau_h;
  double d_dt_j = (j_inf - j) / tau_j;
  double d_dt_d = (d_inf - d) / tau_d;
  double d_dt_f = (f_inf - f) / tau_f;
  double d_dt_f2 = (f2_inf - f2) / tau_f2;
  double d_dt_fCass = (fCass_inf - fCass) / tau_fCass;
  double d_dt_s = (s_inf - s) / tau_s;
  double d_dt_r = (r_inf - r) / tau_r;
  */

  double d_dt_R_prime = ((-k2) * Ca_ss * R_prime) + (k4 * (1.0 - R_prime));
  double d_dt_Ca_i = Ca_i_bufc*(((((Ileak-Iup)*Vsr)/Vc)+Ixfer)-((((IbCa + IpCa)-(2.0*INaCa))*Cm)/(2.0*Vc*F)));
  double d_dt_Ca_SR = Ca_sr_bufsr*(Iup-(Irel+Ileak));
  double d_dt_Ca_ss = Ca_ss_bufss*(((((-ICaL)*Cm)/(2.0*V_ss*F))+((Irel*Vsr)/V_ss))-((Ixfer*Vc)/V_ss));
  double d_dt_Na_i = ((-(INa + IbNa + (3.0 * INaK) + (3.0 * INaCa))) / (Vc * F)) * Cm;
  double d_dt_K_i = ((-((IK1 + Ito + IKr + IKs + IpK + Istim) - (2.0 * INaK))) / (Vc * F)) * Cm;

  // Compute active tension
  //const double kTa = 0.5; //0.5; //13.7;
  //const double Vr = -85.423;
  //double d_dt_Ta = active_stress_twitch(V) * (kTa * (V - Vr) - Ta);
  //double d_dt_Ta = active_stress_twitch(V) * (kTa * (V-Vr) - Ta);
  
  double d_dt_aux;
  double d_dt_Ta;  
  
  // NOVO TESTE    
  double lambda = 1.0;
  double a7 = 0.7;
  double a6 = 6.0;
  double a4 = 650.0;  // ms
  double td = 90.0;   // ms
  double b = 1.0;     // ms
  double ld = -500.0; //
  double T0 = 80; //180;    // kPa
  double tr0 = 80.0;  
  double Fiso = MYMAX(0.0, tanh(a6-(lambda-a7)));  
  double te = b*(lambda - ld);
  double tr = tr0 + a4*(1-Fiso);    
  double Ftwitch = 0.0;
    
  //if (fabs(Istim) != 0.0)
  //    elapsedTime = time;
  //cout << "V = " << V << "\tAUX2 = " << aux2 << endl;
  
  double delay = 0.0; //10.0;
  
  if(V > 0 && aux2 <= 0 )
  {
      elapsedTime = time + delay;
      d_dt_aux = time + delay;
      marked = true;
  }
  else
      d_dt_aux = aux2;
  
  double ts;
  if (time >= elapsedTime && d_dt_aux != 0)
    ts = time - d_dt_aux;
  else
    ts = 0.0;
  
  if( ts < te ) Ftwitch = pow(tanh(ts/tr),2) * pow(tanh((te-ts)/td),2);
  
  d_dt_Ta = T0 * Fiso * Ftwitch;
  
  // KERCKHOFFS
  
  
  /*
  // CADIM
  const double kTa = 500.0;
  const double err = 5e-3;
  const double Vthreshold   = -80.0;//0.05;
  const double auxThreshold = 0.05;
  
  double myCond   = rateswitchNP(V, aux2, Vthreshold, auxThreshold);
  double d_dt_Ta  = myCond * kTa*(((V+85.914)/(14.273+85.914)) - Ta/(kTa*0.1));  
  double d_dt_aux = V<Vthreshold?-aux2:err;
   */

  // Update values
  values[ 0] = d_dt_V;
  
//  values[ 1] = d_dt_Xr1;
//  values[ 2] = d_dt_Xr2;
//  values[ 3] = d_dt_Xs;
//  values[ 4] = d_dt_m;
//  values[ 5] = d_dt_h;
//  values[ 6] = d_dt_j;
//  values[ 7] = d_dt_d;
//  values[ 8] = d_dt_f;
//  values[ 9] = d_dt_f2;
//  values[10] = d_dt_fCass;
//  values[11] = d_dt_s;
//  values[12] = d_dt_r;
  
  const double dt = 0.05;
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
 
  values[19] = d_dt_Ta;  
  values[20] = d_dt_aux;
}

inline double TenTusscherTa::rateswitchNP(double v, double aux, double Vthreshold, double auxThreshold)
{
    if (v < Vthreshold)
      return 0.1;
    else
      if(aux<auxThreshold) return 0.0;
    
    return 0.01;
}

inline double TenTusscherTa::active_stress_twitch(double v)
{
  // My adjustments
  const double e0   = 1.0;
  const double einf = 0.1;
  const double psi  = 0.3;
  //const double vbar = -85.423;
  const double vs = -30.0;
  return e0 + (einf-e0) * exp(-exp(-psi*(v-vs))); 
  
  // Eriksson
  //const double e0   = 1.0;
  //const double einf = 0.1;
  //const double psi  = 0.1;
  //const double vbar = -80.0;
  //return e0 + (einf-e0)*exp(-exp(-psi*(v-vbar))); 

  // Kuhl
  //const double e0   = 0.1;
  //const double einf = 1.0;
  //const double psi  = 0.1;
  //const double vbar = 0.0;
  //return e0 + (einf-e0) * exp(-exp(-psi*(v-vbar))); 
  
  // Nash-Panfilov
  /*if(abs(v) < abs(-30))
    return (0.01);
  else 
    return 1.0;
  */
}
