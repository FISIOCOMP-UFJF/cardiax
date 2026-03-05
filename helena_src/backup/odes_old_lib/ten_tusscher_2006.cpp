#include "ten_Tusscher_2006.h"


// Setting initial conditions to the ten Tusscher 2006 model
void setIC_ode_cpu(Real *sv, int cellID, int NEQ = 19)
{
   // Vector containing the initial conditions of the cell
   Real IC[NEQ];

   // V; Units: millivolt; Initial value: -85.423
   IC[ 0] = -85.423;
   // Xr1; Units: dimensionless; Initial value: 0.0165
   IC[ 1] = 0.0165;
   // Xr2; Units: dimensionless; Initial value: 0.473
   IC[ 2] = 0.473;
   // Xs; Units: dimensionless; Initial value: 0.0174
   IC[ 3] = 0.0174;
   // m; Units: dimensionless; Initial value: 0.00165
   IC[ 4] = 0.00165;
   // h; Units: dimensionless; Initial value: 0.749
   IC[ 5] = 0.749;
   // j; Units: dimensionless; Initial value: 0.6788
   IC[ 6] = 0.6788;
   // d; Units: dimensionless; Initial value: 3.288e-5
   IC[ 7] = 3.288e-5;
   // f; Units: dimensionless; Initial value: 0.7026
   IC[ 8] = 0.7026;
   // f2; Units: dimensionless; Initial value: 0.9526
   IC[ 9] = 0.9526;
   // fCass; Units: dimensionless; Initial value: 0.9942
   IC[10] = 0.9942;
   // s; Units: dimensionless; Initial value: 0.999998
   IC[11] = 0.999998;
   // r; Units: dimensionless; Initial value: 2.347e-8
   IC[12] = 2.347e-8;
   // Ca_i; Units: millimolar; Initial value: 0.000153
   IC[13] = 0.000153;
   // Ca_SR; Units: millimolar; Initial value: 4.272
   IC[14] = 4.272;
   // Ca_ss; Units: millimolar; Initial value: 0.00042
   IC[15] = 0.00042;
   // R_prime; Units: dimensionless; Initial value: 0.8978
   IC[16] = 0.8978;
   // Na_i; Units: millimolar; Initial value: 10.132
   IC[17] = 10.132;
   // K_i; Units: millimolar; Initial value: 138.52
   IC[18] = 138.52;

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = IC[i];
}

void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, int NEQ, int col)
{
   solve_Forward_Euler_cpu(time, sv, dt, NEQ, cellID, col);
}

// The Forward Euler Method
void solve_Forward_Euler_cpu(Real time, Real *sv, Real dt, int NEQ, int cellID, int col)
{
   // sv contains the membrane state variables (NEQ) of all cells in the 'tissue' (ni*nj)
   // cellID defines where the state variable starts in sv

   Real rY[NEQ], rDY[NEQ];

   for(int i = 0; i < NEQ; i++)
      rY[i] = sv[cellID+i];

   RHS_ten_Tusscher_2006_cpu(time, rY, rDY, col);

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = dt*rDY[i] + rY[i];
}

// The ten Tusscher 2006 model
void RHS_ten_Tusscher_2006_cpu(Real time, Real *rY, Real *rDY, int col)
{
   Real var_membrane__V = rY[0];
   Real var_rapid_time_dependent_potassium_current_Xr1_gate__Xr1 = rY[1];
   Real var_rapid_time_dependent_potassium_current_Xr2_gate__Xr2 = rY[2];
   Real var_slow_time_dependent_potassium_current_Xs_gate__Xs = rY[3];
   Real var_fast_sodium_current_m_gate__m = rY[4];
   Real var_fast_sodium_current_h_gate__h = rY[5];
   Real var_fast_sodium_current_j_gate__j = rY[6];
   Real var_L_type_Ca_current_d_gate__d = rY[7];
   Real var_L_type_Ca_current_f_gate__f = rY[8];
   Real var_L_type_Ca_current_f2_gate__f2 = rY[9];
   Real var_L_type_Ca_current_fCass_gate__fCass = rY[10];
   Real var_transient_outward_current_s_gate__s = rY[11];
   Real var_transient_outward_current_r_gate__r = rY[12];
   Real var_calcium_dynamics__Ca_i = rY[13];
   Real var_calcium_dynamics__Ca_SR = rY[14];
   Real var_calcium_dynamics__Ca_ss = rY[15];
   Real var_calcium_dynamics__R_prime = rY[16];
   Real var_sodium_dynamics__Na_i = rY[17];
   Real var_potassium_dynamics__K_i = rY[18];

   // Inputs:
   // Time units: ms
   Real var_environment__time = time;

   // Mathematics
   const Real var_membrane__R = 8314.472;
   const Real var_membrane__T = 310.0;
   const Real var_membrane__F = 96485.3415;
   const Real var_membrane__Cm = 0.185;
   const Real var_membrane__V_c = 0.016404;
   Real var_membrane__time = var_environment__time;
   Real var_inward_rectifier_potassium_current__V = var_membrane__V;
   Real var_reversal_potentials__K_i = var_potassium_dynamics__K_i;
   Real var_reversal_potentials__R = var_membrane__R;
   Real var_reversal_potentials__T = var_membrane__T;
   Real var_reversal_potentials__F = var_membrane__F;
   const Real var_potassium_dynamics__K_o = 5.4;
   Real var_reversal_potentials__K_o = var_potassium_dynamics__K_o;
   Real var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
   Real var_inward_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
   Real var_inward_rectifier_potassium_current__beta_K1 = ((3.0 * exp(0.0002 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) + 100.0))) + exp(0.1 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) - 10.0))) / (1.0 + exp((-0.5) * (var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K)));
   Real var_inward_rectifier_potassium_current__alpha_K1 = 0.1 / (1.0 + exp(0.06 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) - 200.0)));
   Real var_inward_rectifier_potassium_current__xK1_inf = var_inward_rectifier_potassium_current__alpha_K1 / (var_inward_rectifier_potassium_current__alpha_K1 + var_inward_rectifier_potassium_current__beta_K1);
   const Real var_inward_rectifier_potassium_current__g_K1 = 5.405;
   Real var_inward_rectifier_potassium_current__i_K1 = var_inward_rectifier_potassium_current__g_K1 * var_inward_rectifier_potassium_current__xK1_inf * (var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K);
   Real var_membrane__i_K1 = var_inward_rectifier_potassium_current__i_K1;
   Real var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
   Real var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
   const Real var_transient_outward_current__g_to = 0.294;
   Real var_transient_outward_current__V = var_membrane__V;
   Real var_transient_outward_current__E_K = var_reversal_potentials__E_K;
   Real var_transient_outward_current__i_to = var_transient_outward_current__g_to * var_transient_outward_current__r * var_transient_outward_current__s * (var_transient_outward_current__V - var_transient_outward_current__E_K);
   Real var_membrane__i_to = var_transient_outward_current__i_to;
   Real var_rapid_time_dependent_potassium_current__Xr1 = var_rapid_time_dependent_potassium_current_Xr1_gate__Xr1;
   Real var_rapid_time_dependent_potassium_current__V = var_membrane__V;
   Real var_rapid_time_dependent_potassium_current__Xr2 = var_rapid_time_dependent_potassium_current_Xr2_gate__Xr2;
   Real var_rapid_time_dependent_potassium_current__K_o = var_potassium_dynamics__K_o;
   Real var_rapid_time_dependent_potassium_current__E_K = var_reversal_potentials__E_K;
   const Real var_rapid_time_dependent_potassium_current__g_Kr = 0.153;
   Real var_rapid_time_dependent_potassium_current__i_Kr = var_rapid_time_dependent_potassium_current__g_Kr * var_rapid_time_dependent_potassium_current__Xr1 * var_rapid_time_dependent_potassium_current__Xr2 * (var_rapid_time_dependent_potassium_current__V - var_rapid_time_dependent_potassium_current__E_K) * sqrt(var_rapid_time_dependent_potassium_current__K_o / 5.4);
   Real var_membrane__i_Kr = var_rapid_time_dependent_potassium_current__i_Kr;
   const Real var_slow_time_dependent_potassium_current__g_Ks = 0.098;
   const Real var_sodium_dynamics__Na_o = 140.0;
   Real var_reversal_potentials__Na_o = var_sodium_dynamics__Na_o;
   Real var_reversal_potentials__Na_i = var_sodium_dynamics__Na_i;
   const Real var_reversal_potentials__P_kna = 0.03;
   Real var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
   Real var_slow_time_dependent_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
   Real var_slow_time_dependent_potassium_current__V = var_membrane__V;
   Real var_slow_time_dependent_potassium_current__Xs = var_slow_time_dependent_potassium_current_Xs_gate__Xs;
   Real var_slow_time_dependent_potassium_current__i_Ks = var_slow_time_dependent_potassium_current__g_Ks * pow(var_slow_time_dependent_potassium_current__Xs, 2) * (var_slow_time_dependent_potassium_current__V - var_slow_time_dependent_potassium_current__E_Ks);
   Real var_membrane__i_Ks = var_slow_time_dependent_potassium_current__i_Ks;
   Real var_L_type_Ca_current__Ca_ss = var_calcium_dynamics__Ca_ss;
   const Real var_L_type_Ca_current__g_CaL = 3.98e-05;
   Real var_L_type_Ca_current__f = var_L_type_Ca_current_f_gate__f;
   Real var_L_type_Ca_current__d = var_L_type_Ca_current_d_gate__d;
   Real var_L_type_Ca_current__F = var_membrane__F;
   Real var_L_type_Ca_current__f2 = var_L_type_Ca_current_f2_gate__f2;
   Real var_L_type_Ca_current__fCass = var_L_type_Ca_current_fCass_gate__fCass;
   Real var_L_type_Ca_current__V = var_membrane__V;
   Real var_L_type_Ca_current__T = var_membrane__T;
   const Real var_calcium_dynamics__Ca_o = 2.0;
   Real var_L_type_Ca_current__Ca_o = var_calcium_dynamics__Ca_o;
   Real var_L_type_Ca_current__R = var_membrane__R;
   Real var_L_type_Ca_current__i_CaL = (var_L_type_Ca_current__V < 15.0-1.0e-5 || var_L_type_Ca_current__V > 15.0+1.0e-5) ? ((((var_L_type_Ca_current__g_CaL * var_L_type_Ca_current__d * var_L_type_Ca_current__f * var_L_type_Ca_current__f2 * var_L_type_Ca_current__fCass * 4.0 * (var_L_type_Ca_current__V - 15.0) * pow(var_L_type_Ca_current__F, 2)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((0.25 * var_L_type_Ca_current__Ca_ss * exp((2.0 * (var_L_type_Ca_current__V - 15.0) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__Ca_o)) / (exp((2.0 * (var_L_type_Ca_current__V - 15.0) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0)) : var_L_type_Ca_current__g_CaL * var_L_type_Ca_current__d * var_L_type_Ca_current__f * var_L_type_Ca_current__f2 * var_L_type_Ca_current__fCass * 2.0 * var_L_type_Ca_current__F * (0.25 * var_L_type_Ca_current__Ca_ss  - var_L_type_Ca_current__Ca_o);
   Real var_membrane__i_CaL = var_L_type_Ca_current__i_CaL;
   Real var_sodium_potassium_pump_current__Na_i = var_sodium_dynamics__Na_i;
   Real var_sodium_potassium_pump_current__R = var_membrane__R;
   Real var_sodium_potassium_pump_current__T = var_membrane__T;
   Real var_sodium_potassium_pump_current__V = var_membrane__V;
   const Real var_sodium_potassium_pump_current__K_mk = 1.0;
   const Real var_sodium_potassium_pump_current__P_NaK = 2.724;
   const Real var_sodium_potassium_pump_current__K_mNa = 40.0;
   Real var_sodium_potassium_pump_current__F = var_membrane__F;
   Real var_sodium_potassium_pump_current__K_o = var_potassium_dynamics__K_o;
   Real var_sodium_potassium_pump_current__i_NaK = ((((var_sodium_potassium_pump_current__P_NaK * var_sodium_potassium_pump_current__K_o) / (var_sodium_potassium_pump_current__K_o + var_sodium_potassium_pump_current__K_mk)) * var_sodium_potassium_pump_current__Na_i) / (var_sodium_potassium_pump_current__Na_i + var_sodium_potassium_pump_current__K_mNa)) / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump_current__V * var_sodium_potassium_pump_current__F) / (var_sodium_potassium_pump_current__R * var_sodium_potassium_pump_current__T))) + (0.0353 * exp(((-var_sodium_potassium_pump_current__V) * var_sodium_potassium_pump_current__F) / (var_sodium_potassium_pump_current__R * var_sodium_potassium_pump_current__T))));
   Real var_membrane__i_NaK = var_sodium_potassium_pump_current__i_NaK;
   Real var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
   Real var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
   const Real var_fast_sodium_current__g_Na = 14.838;
   Real var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
   Real var_fast_sodium_current__V = var_membrane__V;
   Real var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
   Real var_fast_sodium_current__E_Na = var_reversal_potentials__E_Na;
   Real var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
   Real var_membrane__i_Na = var_fast_sodium_current__i_Na;
   Real var_sodium_background_current__V = var_membrane__V;
   Real var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
   const Real var_sodium_background_current__g_bna = 0.00029;
   Real var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
   Real var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
   const Real var_sodium_calcium_exchanger_current__alpha = 2.5;
   const Real var_sodium_calcium_exchanger_current__gamma = 0.35;
   const Real var_sodium_calcium_exchanger_current__K_sat = 0.1;
   const Real var_sodium_calcium_exchanger_current__Km_Ca = 1.38;
   const Real var_sodium_calcium_exchanger_current__K_NaCa = 1000.0;
   Real var_sodium_calcium_exchanger_current__F = var_membrane__F;
   Real var_sodium_calcium_exchanger_current__Ca_i = var_calcium_dynamics__Ca_i;
   Real var_sodium_calcium_exchanger_current__Ca_o = var_calcium_dynamics__Ca_o;
   Real var_sodium_calcium_exchanger_current__V = var_membrane__V;
   Real var_sodium_calcium_exchanger_current__R = var_membrane__R;
   const Real var_sodium_calcium_exchanger_current__Km_Nai = 87.5;
   Real var_sodium_calcium_exchanger_current__Na_o = var_sodium_dynamics__Na_o;
   Real var_sodium_calcium_exchanger_current__Na_i = var_sodium_dynamics__Na_i;
   Real var_sodium_calcium_exchanger_current__T = var_membrane__T;
   Real var_sodium_calcium_exchanger_current__i_NaCa = (var_sodium_calcium_exchanger_current__K_NaCa * ((exp((var_sodium_calcium_exchanger_current__gamma * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)) * pow(var_sodium_calcium_exchanger_current__Na_i, 3) * var_sodium_calcium_exchanger_current__Ca_o) - (exp(((var_sodium_calcium_exchanger_current__gamma - 1.0) * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)) * pow(var_sodium_calcium_exchanger_current__Na_o, 3) * var_sodium_calcium_exchanger_current__Ca_i * var_sodium_calcium_exchanger_current__alpha))) / ((pow(var_sodium_calcium_exchanger_current__Km_Nai, 3) + pow(var_sodium_calcium_exchanger_current__Na_o, 3)) * (var_sodium_calcium_exchanger_current__Km_Ca + var_sodium_calcium_exchanger_current__Ca_o) * (1.0 + (var_sodium_calcium_exchanger_current__K_sat * exp(((var_sodium_calcium_exchanger_current__gamma - 1.0) * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)))));
   Real var_membrane__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa;
   Real var_reversal_potentials__Ca_o = var_calcium_dynamics__Ca_o;
   Real var_reversal_potentials__Ca_i = var_calcium_dynamics__Ca_i;
   Real var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
   Real var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
   const Real var_calcium_background_current__g_bca = 0.000592;
   Real var_calcium_background_current__V = var_membrane__V;
   Real var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
   Real var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
   const Real var_potassium_pump_current__g_pK = 0.0146;
   Real var_potassium_pump_current__V = var_membrane__V;
   Real var_potassium_pump_current__E_K = var_reversal_potentials__E_K;
   Real var_potassium_pump_current__i_p_K = (var_potassium_pump_current__g_pK * (var_potassium_pump_current__V - var_potassium_pump_current__E_K)) / (1.0 + exp((25.0 - var_potassium_pump_current__V) / 5.98));
   Real var_membrane__i_p_K = var_potassium_pump_current__i_p_K;
   const Real var_calcium_pump_current__K_pCa = 0.0005;
   Real var_calcium_pump_current__Ca_i = var_calcium_dynamics__Ca_i;
   const Real var_calcium_pump_current__g_pCa = 0.1238;
   Real var_calcium_pump_current__i_p_Ca = (var_calcium_pump_current__g_pCa * var_calcium_pump_current__Ca_i) / (var_calcium_pump_current__Ca_i + var_calcium_pump_current__K_pCa);
   Real var_membrane__i_p_Ca = var_calcium_pump_current__i_p_Ca;
   const Real var_membrane__stim_end = 100000.0;
   //const real var_membrane__stim_amplitude =  -52.0;
   const Real var_membrane__stim_amplitude =  0.0;
   const Real var_membrane__stim_duration = 1.0;
   const Real var_membrane__stim_period = 1000.0;
   const Real var_membrane__stim_start = 100.0 + (Real)col;
   Real var_membrane__i_Stim = ((var_membrane__time >= var_membrane__stim_start) && (var_membrane__time <= var_membrane__stim_end) && (((var_membrane__time - var_membrane__stim_start) - (floor((var_membrane__time - var_membrane__stim_start) / var_membrane__stim_period) * var_membrane__stim_period)) <= var_membrane__stim_duration)) ? var_membrane__stim_amplitude : 0.0;
   Real var_rapid_time_dependent_potassium_current_Xr1_gate__V = var_rapid_time_dependent_potassium_current__V;
   Real var_rapid_time_dependent_potassium_current_Xr1_gate__xr1_inf = 1.0 / (1.0 + exp(((-26.0) - var_rapid_time_dependent_potassium_current_Xr1_gate__V) / 7.0));
   Real var_rapid_time_dependent_potassium_current_Xr1_gate__alpha_xr1 = 450.0 / (1.0 + exp(((-45.0) - var_rapid_time_dependent_potassium_current_Xr1_gate__V) / 10.0));
   Real var_rapid_time_dependent_potassium_current_Xr1_gate__beta_xr1 = 6.0 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr1_gate__V + 30.0) / 11.5));
   Real var_rapid_time_dependent_potassium_current_Xr1_gate__tau_xr1 = 1.0 * var_rapid_time_dependent_potassium_current_Xr1_gate__alpha_xr1 * var_rapid_time_dependent_potassium_current_Xr1_gate__beta_xr1;
   Real var_rapid_time_dependent_potassium_current_Xr2_gate__V = var_rapid_time_dependent_potassium_current__V;
   Real var_rapid_time_dependent_potassium_current_Xr2_gate__xr2_inf = 1.0 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr2_gate__V + 88.0) / 24.0));
   Real var_rapid_time_dependent_potassium_current_Xr2_gate__alpha_xr2 = 3.0 / (1.0 + exp(((-60.0) - var_rapid_time_dependent_potassium_current_Xr2_gate__V) / 20.0));
   Real var_rapid_time_dependent_potassium_current_Xr2_gate__beta_xr2 = 1.12 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr2_gate__V - 60.0) / 20.0));
   Real var_rapid_time_dependent_potassium_current_Xr2_gate__tau_xr2 = 1.0 * var_rapid_time_dependent_potassium_current_Xr2_gate__alpha_xr2 * var_rapid_time_dependent_potassium_current_Xr2_gate__beta_xr2;
   Real var_slow_time_dependent_potassium_current_Xs_gate__V = var_slow_time_dependent_potassium_current__V;
   Real var_slow_time_dependent_potassium_current_Xs_gate__xs_inf = 1.0 / (1.0 + exp(((-5.0) - var_slow_time_dependent_potassium_current_Xs_gate__V) / 14.0));
   Real var_slow_time_dependent_potassium_current_Xs_gate__alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - var_slow_time_dependent_potassium_current_Xs_gate__V) / 6.0));
   Real var_slow_time_dependent_potassium_current_Xs_gate__beta_xs = 1.0 / (1.0 + exp((var_slow_time_dependent_potassium_current_Xs_gate__V - 35.0) / 15.0));
   Real var_slow_time_dependent_potassium_current_Xs_gate__tau_xs = (1.0 * var_slow_time_dependent_potassium_current_Xs_gate__alpha_xs * var_slow_time_dependent_potassium_current_Xs_gate__beta_xs) + 80.0;
   Real var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
   Real var_fast_sodium_current_m_gate__m_inf = 1.0 / pow(1.0 + exp(((-56.86) - var_fast_sodium_current_m_gate__V) / 9.03), 2);
   Real var_fast_sodium_current_m_gate__alpha_m = 1.0 / (1.0 + exp(((-60.0) - var_fast_sodium_current_m_gate__V) / 5.0));
   Real var_fast_sodium_current_m_gate__beta_m = (0.1 / (1.0 + exp((var_fast_sodium_current_m_gate__V + 35.0) / 5.0))) + (0.1 / (1.0 + exp((var_fast_sodium_current_m_gate__V - 50.0) / 200.0)));
   Real var_fast_sodium_current_m_gate__tau_m = 1.0 * var_fast_sodium_current_m_gate__alpha_m * var_fast_sodium_current_m_gate__beta_m;
   Real var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
   Real var_fast_sodium_current_h_gate__h_inf = 1.0 / pow(1.0 + exp((var_fast_sodium_current_h_gate__V + 71.55) / 7.43), 2);
   Real var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.057 * exp((-(var_fast_sodium_current_h_gate__V + 80.0)) / 6.8)) : 0.0;
   Real var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((2.7 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.3485 * var_fast_sodium_current_h_gate__V))) : (0.77 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1)))));
   Real var_fast_sodium_current_h_gate__tau_h = 1.0 / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h);
   Real var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
   Real var_fast_sodium_current_j_gate__j_inf = 1.0 / pow(1.0 + exp((var_fast_sodium_current_j_gate__V + 71.55) / 7.43), 2);
   Real var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((((((-25428.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (6.948e-06 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / 1.0) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0;
   Real var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.02424 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.6 * exp(0.057 * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))));
   Real var_fast_sodium_current_j_gate__tau_j = 1.0 / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j);
   Real var_L_type_Ca_current_d_gate__V = var_L_type_Ca_current__V;
   Real var_L_type_Ca_current_d_gate__d_inf = 1.0 / (1.0 + exp(((-8.0) - var_L_type_Ca_current_d_gate__V) / 7.5));
   Real var_L_type_Ca_current_d_gate__alpha_d = (1.4 / (1.0 + exp(((-35.0) - var_L_type_Ca_current_d_gate__V) / 13.0))) + 0.25;
   Real var_L_type_Ca_current_d_gate__beta_d = 1.4 / (1.0 + exp((var_L_type_Ca_current_d_gate__V + 5.0) / 5.0));
   Real var_L_type_Ca_current_d_gate__gamma_d = 1.0 / (1.0 + exp((50.0 - var_L_type_Ca_current_d_gate__V) / 20.0));
   Real var_L_type_Ca_current_d_gate__tau_d = (1.0 * var_L_type_Ca_current_d_gate__alpha_d * var_L_type_Ca_current_d_gate__beta_d) + var_L_type_Ca_current_d_gate__gamma_d;
   Real var_L_type_Ca_current_f_gate__V = var_L_type_Ca_current__V;
   Real var_L_type_Ca_current_f_gate__f_inf = 1.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 20.0) / 7.0));
   Real var_L_type_Ca_current_f_gate__tau_f = (1102.5 * exp((-pow(var_L_type_Ca_current_f_gate__V + 27.0, 2)) / 225.0)) + (200.0 / (1.0 + exp((13.0 - var_L_type_Ca_current_f_gate__V) / 10.0))) + (180.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 30.0) / 10.0))) + 20.0;
   Real var_L_type_Ca_current_f2_gate__V = var_L_type_Ca_current__V;
   Real var_L_type_Ca_current_f2_gate__f2_inf = (0.67 / (1.0 + exp((var_L_type_Ca_current_f2_gate__V + 35.0) / 7.0))) + 0.33;
   Real var_L_type_Ca_current_f2_gate__tau_f2 = (562.0 * exp((-pow(var_L_type_Ca_current_f2_gate__V + 27.0, 2)) / 240.0)) + (31.0 / (1.0 + exp((25.0 - var_L_type_Ca_current_f2_gate__V) / 10.0))) + (80.0 / (1.0 + exp((var_L_type_Ca_current_f2_gate__V + 30.0) / 10.0)));
   Real var_L_type_Ca_current_fCass_gate__Ca_ss = var_L_type_Ca_current__Ca_ss;
   Real var_L_type_Ca_current_fCass_gate__fCass_inf = (0.6 / (1.0 + pow(var_L_type_Ca_current_fCass_gate__Ca_ss / 0.05, 2))) + 0.4;
   Real var_L_type_Ca_current_fCass_gate__tau_fCass = (80.0 / (1.0 + pow(var_L_type_Ca_current_fCass_gate__Ca_ss / 0.05, 2))) + 2.0;
   Real var_transient_outward_current_s_gate__V = var_transient_outward_current__V;
   Real var_transient_outward_current_s_gate__s_inf = 1.0 / (1.0 + exp((var_transient_outward_current_s_gate__V + 20.0) / 5.0));
   Real var_transient_outward_current_s_gate__tau_s = (85.0 * exp((-pow(var_transient_outward_current_s_gate__V + 45.0, 2)) / 320.0)) + (5.0 / (1.0 + exp((var_transient_outward_current_s_gate__V - 20.0) / 5.0))) + 3.0;
   Real var_transient_outward_current_r_gate__V = var_transient_outward_current__V;
   Real var_transient_outward_current_r_gate__r_inf = 1.0 / (1.0 + exp((20.0 - var_transient_outward_current_r_gate__V) / 6.0));
   Real var_transient_outward_current_r_gate__tau_r = (9.5 * exp((-pow(var_transient_outward_current_r_gate__V + 40.0, 2)) / 1800.0)) + 0.8;
   const Real var_calcium_dynamics__V_rel = 0.102;
   const Real var_calcium_dynamics__k1_prime = 0.15;
   const Real var_calcium_dynamics__max_sr = 2.5;
   const Real var_calcium_dynamics__EC = 1.5;
   const Real var_calcium_dynamics__min_sr = 1.0;
   Real var_calcium_dynamics__kcasr = var_calcium_dynamics__max_sr - ((var_calcium_dynamics__max_sr - var_calcium_dynamics__min_sr) / (1.0 + pow(var_calcium_dynamics__EC / var_calcium_dynamics__Ca_SR, 2)));
   Real var_calcium_dynamics__k1 = var_calcium_dynamics__k1_prime / var_calcium_dynamics__kcasr;
   const Real var_calcium_dynamics__k3 = 0.06;
   Real var_calcium_dynamics__O = (var_calcium_dynamics__k1 * pow(var_calcium_dynamics__Ca_ss, 2) * var_calcium_dynamics__R_prime) / (var_calcium_dynamics__k3 + (var_calcium_dynamics__k1 * pow(var_calcium_dynamics__Ca_ss, 2)));
   Real var_calcium_dynamics__i_rel = var_calcium_dynamics__V_rel * var_calcium_dynamics__O * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_ss);
   const Real var_calcium_dynamics__Vmax_up = 0.006375;
   const Real var_calcium_dynamics__K_up = 0.00025;
   Real var_calcium_dynamics__i_up = var_calcium_dynamics__Vmax_up / (1.0 + (pow(var_calcium_dynamics__K_up, 2) / pow(var_calcium_dynamics__Ca_i, 2)));
   const Real var_calcium_dynamics__V_leak = 0.00036;
   Real var_calcium_dynamics__i_leak = var_calcium_dynamics__V_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i);
   const Real var_calcium_dynamics__V_xfer = 0.0038;
   Real var_calcium_dynamics__i_xfer = var_calcium_dynamics__V_xfer * (var_calcium_dynamics__Ca_ss - var_calcium_dynamics__Ca_i);
   const Real var_calcium_dynamics__k2_prime = 0.045;
   Real var_calcium_dynamics__k2 = var_calcium_dynamics__k2_prime * var_calcium_dynamics__kcasr;
   const Real var_calcium_dynamics__k4 = 0.005;
   const Real var_calcium_dynamics__Buf_c = 0.2;
   const Real var_calcium_dynamics__K_buf_c = 0.001;
   Real var_calcium_dynamics__Ca_i_bufc = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_c * var_calcium_dynamics__K_buf_c) / pow(var_calcium_dynamics__Ca_i + var_calcium_dynamics__K_buf_c, 2)));
   const Real var_calcium_dynamics__K_buf_sr = 0.3;
   const Real var_calcium_dynamics__Buf_sr = 10.0;
   Real var_calcium_dynamics__Ca_sr_bufsr = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_sr * var_calcium_dynamics__K_buf_sr) / pow(var_calcium_dynamics__Ca_SR + var_calcium_dynamics__K_buf_sr, 2)));
   const Real var_calcium_dynamics__Buf_ss = 0.4;
   const Real var_calcium_dynamics__K_buf_ss = 0.00025;
   Real var_calcium_dynamics__Ca_ss_bufss = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_ss * var_calcium_dynamics__K_buf_ss) / pow(var_calcium_dynamics__Ca_ss + var_calcium_dynamics__K_buf_ss, 2)));
   const Real var_calcium_dynamics__V_sr = 0.001094;
   const Real var_calcium_dynamics__V_ss = 5.468e-05;
   Real var_calcium_dynamics__V_c = var_membrane__V_c;
   Real var_calcium_dynamics__F = var_membrane__F;
   Real var_calcium_dynamics__Cm = var_membrane__Cm;
   Real var_calcium_dynamics__i_CaL = var_L_type_Ca_current__i_CaL;
   Real var_calcium_dynamics__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa;
   Real var_calcium_dynamics__i_p_Ca = var_calcium_pump_current__i_p_Ca;
   Real var_calcium_dynamics__i_b_Ca = var_calcium_background_current__i_b_Ca;
   Real var_sodium_dynamics__F = var_membrane__F;
   Real var_sodium_dynamics__Cm = var_membrane__Cm;
   Real var_sodium_dynamics__V_c = var_membrane__V_c;
   Real var_sodium_dynamics__i_Na = var_fast_sodium_current__i_Na;
   Real var_sodium_dynamics__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa;
   Real var_sodium_dynamics__i_NaK = var_sodium_potassium_pump_current__i_NaK;
   Real var_sodium_dynamics__i_b_Na = var_sodium_background_current__i_b_Na;
   Real var_potassium_dynamics__F = var_membrane__F;
   Real var_potassium_dynamics__Cm = var_membrane__Cm;
   Real var_potassium_dynamics__V_c = var_membrane__V_c;
   Real var_potassium_dynamics__i_K1 = var_inward_rectifier_potassium_current__i_K1;
   Real var_potassium_dynamics__i_to = var_transient_outward_current__i_to;
   Real var_potassium_dynamics__i_NaK = var_sodium_potassium_pump_current__i_NaK;
   Real var_potassium_dynamics__i_Kr = var_rapid_time_dependent_potassium_current__i_Kr;
   Real var_potassium_dynamics__i_Ks = var_slow_time_dependent_potassium_current__i_Ks;
   Real var_potassium_dynamics__i_p_K = var_potassium_pump_current__i_p_K;
   Real var_potassium_dynamics__i_Stim = var_membrane__i_Stim;
   Real d_dt_membrane__V = -(var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_CaL + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_NaCa + var_membrane__i_b_Ca + var_membrane__i_p_K + var_membrane__i_p_Ca + var_membrane__i_Stim);
   Real d_dt_rapid_time_dependent_potassium_current_Xr1_gate__Xr1 = (var_rapid_time_dependent_potassium_current_Xr1_gate__xr1_inf - var_rapid_time_dependent_potassium_current_Xr1_gate__Xr1) / var_rapid_time_dependent_potassium_current_Xr1_gate__tau_xr1;
   Real d_dt_rapid_time_dependent_potassium_current_Xr2_gate__Xr2 = (var_rapid_time_dependent_potassium_current_Xr2_gate__xr2_inf - var_rapid_time_dependent_potassium_current_Xr2_gate__Xr2) / var_rapid_time_dependent_potassium_current_Xr2_gate__tau_xr2;
   Real d_dt_slow_time_dependent_potassium_current_Xs_gate__Xs = (var_slow_time_dependent_potassium_current_Xs_gate__xs_inf - var_slow_time_dependent_potassium_current_Xs_gate__Xs) / var_slow_time_dependent_potassium_current_Xs_gate__tau_xs;
   Real d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__m_inf - var_fast_sodium_current_m_gate__m) / var_fast_sodium_current_m_gate__tau_m;
   Real d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__h_inf - var_fast_sodium_current_h_gate__h) / var_fast_sodium_current_h_gate__tau_h;
   Real d_dt_fast_sodium_current_j_gate__j = (var_fast_sodium_current_j_gate__j_inf - var_fast_sodium_current_j_gate__j) / var_fast_sodium_current_j_gate__tau_j;
   Real d_dt_L_type_Ca_current_d_gate__d = (var_L_type_Ca_current_d_gate__d_inf - var_L_type_Ca_current_d_gate__d) / var_L_type_Ca_current_d_gate__tau_d;
   Real d_dt_L_type_Ca_current_f_gate__f = (var_L_type_Ca_current_f_gate__f_inf - var_L_type_Ca_current_f_gate__f) / var_L_type_Ca_current_f_gate__tau_f;
   Real d_dt_L_type_Ca_current_f2_gate__f2 = (var_L_type_Ca_current_f2_gate__f2_inf - var_L_type_Ca_current_f2_gate__f2) / var_L_type_Ca_current_f2_gate__tau_f2;
   Real d_dt_L_type_Ca_current_fCass_gate__fCass = (var_L_type_Ca_current_fCass_gate__fCass_inf - var_L_type_Ca_current_fCass_gate__fCass) / var_L_type_Ca_current_fCass_gate__tau_fCass;
   Real d_dt_transient_outward_current_s_gate__s = (var_transient_outward_current_s_gate__s_inf - var_transient_outward_current_s_gate__s) / var_transient_outward_current_s_gate__tau_s;
   Real d_dt_transient_outward_current_r_gate__r = (var_transient_outward_current_r_gate__r_inf - var_transient_outward_current_r_gate__r) / var_transient_outward_current_r_gate__tau_r;
   Real d_dt_calcium_dynamics__R_prime = ((-var_calcium_dynamics__k2) * var_calcium_dynamics__Ca_ss * var_calcium_dynamics__R_prime) + (var_calcium_dynamics__k4 * (1.0 - var_calcium_dynamics__R_prime));
   Real d_dt_calcium_dynamics__Ca_i = var_calcium_dynamics__Ca_i_bufc * (((((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_c) + var_calcium_dynamics__i_xfer) - ((((var_calcium_dynamics__i_b_Ca + var_calcium_dynamics__i_p_Ca) - (2.0 * var_calcium_dynamics__i_NaCa)) * var_calcium_dynamics__Cm) / (2.0 * var_calcium_dynamics__V_c * var_calcium_dynamics__F)));
   Real d_dt_calcium_dynamics__Ca_SR = var_calcium_dynamics__Ca_sr_bufsr * (var_calcium_dynamics__i_up - (var_calcium_dynamics__i_rel + var_calcium_dynamics__i_leak));
   Real d_dt_calcium_dynamics__Ca_ss = var_calcium_dynamics__Ca_ss_bufss * (((((-var_calcium_dynamics__i_CaL) * var_calcium_dynamics__Cm) / (2.0 * var_calcium_dynamics__V_ss * var_calcium_dynamics__F)) + ((var_calcium_dynamics__i_rel * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_ss)) - ((var_calcium_dynamics__i_xfer * var_calcium_dynamics__V_c) / var_calcium_dynamics__V_ss));
   Real d_dt_sodium_dynamics__Na_i = ((-(var_sodium_dynamics__i_Na + var_sodium_dynamics__i_b_Na + (3.0 * var_sodium_dynamics__i_NaK) + (3.0 * var_sodium_dynamics__i_NaCa))) / (var_sodium_dynamics__V_c * var_sodium_dynamics__F)) * var_sodium_dynamics__Cm;
   Real d_dt_potassium_dynamics__K_i = ((-((var_potassium_dynamics__i_K1 + var_potassium_dynamics__i_to + var_potassium_dynamics__i_Kr + var_potassium_dynamics__i_Ks + var_potassium_dynamics__i_p_K + var_potassium_dynamics__i_Stim) - (2.0 * var_potassium_dynamics__i_NaK))) / (var_potassium_dynamics__V_c * var_potassium_dynamics__F)) * var_potassium_dynamics__Cm;

   rDY[ 0] = d_dt_membrane__V;
   rDY[ 1] = d_dt_rapid_time_dependent_potassium_current_Xr1_gate__Xr1;
   rDY[ 2] = d_dt_rapid_time_dependent_potassium_current_Xr2_gate__Xr2;
   rDY[ 3] = d_dt_slow_time_dependent_potassium_current_Xs_gate__Xs;
   rDY[ 4] = d_dt_fast_sodium_current_m_gate__m;
   rDY[ 5] = d_dt_fast_sodium_current_h_gate__h;
   rDY[ 6] = d_dt_fast_sodium_current_j_gate__j;
   rDY[ 7] = d_dt_L_type_Ca_current_d_gate__d;
   rDY[ 8] = d_dt_L_type_Ca_current_f_gate__f;
   rDY[ 9] = d_dt_L_type_Ca_current_f2_gate__f2;
   rDY[10] = d_dt_L_type_Ca_current_fCass_gate__fCass;
   rDY[11] = d_dt_transient_outward_current_s_gate__s;
   rDY[12] = d_dt_transient_outward_current_r_gate__r;
   rDY[13] = d_dt_calcium_dynamics__Ca_i;
   rDY[14] = d_dt_calcium_dynamics__Ca_SR;
   rDY[15] = d_dt_calcium_dynamics__Ca_ss;
   rDY[16] = d_dt_calcium_dynamics__R_prime;
   rDY[17] = d_dt_sodium_dynamics__Na_i;
   rDY[18] = d_dt_potassium_dynamics__K_i;
}
