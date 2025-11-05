/*
 * Luo_Rudy_1991.cpp
 * Created on: 28/10/2009
 */

#include "Luo_Rudy_1991.h"

// Setting initial conditions to the Luo-Rudy 1991 model
void setIC_ode_cpu(Real *sv, int cellID, int NEQ = 8)
{

   // Vector containing the initial conditions of the cell
   Real IC[NEQ];

   // V; Units: millivolt; Initial value: -83.853
   IC[0] = -83.853;
   // m; Units: dimensionless; Initial value: 0.00187018
   IC[1] = 0.00187018;
   // h; Units: dimensionless; Initial value: 0.9804713
   IC[2] = 0.9804713;
   // j; Units: dimensionless; Initial value: 0.98767124
   IC[3] = 0.98767124;
   // d; Units: dimensionless; Initial value: 0.00316354
   IC[4] = 0.00316354;
   // f; Units: dimensionless; Initial value: 0.99427859
   IC[5] = 0.99427859;
   // X; Units: dimensionless; Initial value: 0.16647703
   IC[6] = 0.16647703;
   // Cai; Units: millimolar; Initial value: 0.0002
   IC[7] = 0.0002;

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = IC[i];
}

// Solving the Luo-Rudy 1991 model for each cell in the tissue matrix ni x nj
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

   RHS_Luo_Rudy_1991_cpu(time, rY, rDY, col);

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = dt*rDY[i] + rY[i];
}

// The Luo-Rudy 1991 model
void RHS_Luo_Rudy_1991_cpu(Real time, Real *rY, Real *rDY, int col)
{
   Real var_membrane__V = rY[0];
   // Units: millivolt; Initial value: -83.853
   Real var_fast_sodium_current_m_gate__m = rY[1];
   // Units: dimensionless; Initial value: 0.00187018
   Real var_fast_sodium_current_h_gate__h = rY[2];
   // Units: dimensionless; Initial value: 0.9804713
   Real var_fast_sodium_current_j_gate__j = rY[3];
   // Units: dimensionless; Initial value: 0.98767124
   Real var_slow_inward_current_d_gate__d = rY[4];
   // Units: dimensionless; Initial value: 0.00316354
   Real var_slow_inward_current_f_gate__f = rY[5];
   // Units: dimensionless; Initial value: 0.99427859
   Real var_time_dependent_potassium_current_X_gate__X = rY[6];
   // Units: dimensionless; Initial value: 0.16647703
   Real var_intracellular_calcium_concentration__Cai = rY[7];
   // Units: millimolar; Initial value: 0.0002


   // Inputs:
   // Time units: ms
   Real var_environment__time = time;

   // Mathematics
   const Real var_membrane__R = 8314.0;
   const Real var_membrane__T = 310.0;
   const Real var_membrane__F = 96484.6;
   const Real var_membrane__C = 1.0;
   const Real var_membrane__stim_end = 100000000000.0;
   //const real var_membrane__stim_amplitude =  -25.5;
   const Real var_membrane__stim_amplitude =  0.0;
   Real var_membrane__time = var_environment__time;
   const Real var_membrane__stim_duration = 2.0;
   const Real var_membrane__stim_period = 1000.0;
   const Real var_membrane__stim_start = 100.0 + (Real)col;
   Real var_membrane__I_stim = ((var_membrane__time >= var_membrane__stim_start) && (var_membrane__time <= var_membrane__stim_end) && (((var_membrane__time - var_membrane__stim_start) - (floor((var_membrane__time - var_membrane__stim_start) / var_membrane__stim_period) * var_membrane__stim_period)) <= var_membrane__stim_duration)) ? var_membrane__stim_amplitude : 0.0;
   Real var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
   Real var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
   const Real var_fast_sodium_current__g_Na = 23.0;
   Real var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
   Real var_fast_sodium_current__V = var_membrane__V;
   Real var_fast_sodium_current__R = var_membrane__R;
   Real var_fast_sodium_current__F = var_membrane__F;
   const Real var_ionic_concentrations__Nao = 140.0;
   Real var_fast_sodium_current__Nao = var_ionic_concentrations__Nao;
   const Real var_ionic_concentrations__Nai = 18.0;
   Real var_fast_sodium_current__Nai = var_ionic_concentrations__Nai;
   Real var_fast_sodium_current__T = var_membrane__T;
   Real var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Nao / var_fast_sodium_current__Nai);
   Real var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
   Real var_membrane__i_Na = var_fast_sodium_current__i_Na;
   Real var_slow_inward_current__d = var_slow_inward_current_d_gate__d;
   Real var_slow_inward_current__f = var_slow_inward_current_f_gate__f;
   Real var_slow_inward_current__V = var_membrane__V;
   Real var_slow_inward_current__Cai = var_intracellular_calcium_concentration__Cai;
   Real var_slow_inward_current__E_si = 7.7 - (13.0287 * log(var_slow_inward_current__Cai / 1.0));
   Real var_slow_inward_current__i_si = 0.09 * var_slow_inward_current__d * var_slow_inward_current__f * (var_slow_inward_current__V - var_slow_inward_current__E_si);
   Real var_membrane__i_si = var_slow_inward_current__i_si;
   Real var_time_dependent_potassium_current__V = var_membrane__V;
   Real var_time_dependent_potassium_current_Xi_gate__V = var_time_dependent_potassium_current__V;
   Real var_time_dependent_potassium_current_Xi_gate__Xi = (var_time_dependent_potassium_current_Xi_gate__V > (-100.0)) ? ((var_time_dependent_potassium_current_Xi_gate__V != -77.0)?((2.837 * (exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 77.0)) - 1.0)) / ((var_time_dependent_potassium_current_Xi_gate__V + 77.0) * exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 35.0)))):((2.837 * (exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 77.0+1.0e-6)) - 1.0)) / ((var_time_dependent_potassium_current_Xi_gate__V + 77.0+1.0e-6) * exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 35.0)))) ): 1.0;
   Real var_time_dependent_potassium_current__Xi = var_time_dependent_potassium_current_Xi_gate__Xi;
   const Real var_ionic_concentrations__Ko = 5.4;
   Real var_time_dependent_potassium_current__Ko = var_ionic_concentrations__Ko;
   Real var_time_dependent_potassium_current__g_K = 0.282 * sqrt(var_time_dependent_potassium_current__Ko / 5.4);
   Real var_time_dependent_potassium_current__X = var_time_dependent_potassium_current_X_gate__X;
   const Real var_time_dependent_potassium_current__PR_NaK = 0.01833;
   Real var_time_dependent_potassium_current__F = var_membrane__F;
   Real var_time_dependent_potassium_current__Nao = var_ionic_concentrations__Nao;
   const Real var_ionic_concentrations__Ki = 145.0;
   Real var_time_dependent_potassium_current__Ki = var_ionic_concentrations__Ki;
   Real var_time_dependent_potassium_current__Nai = var_ionic_concentrations__Nai;
   Real var_time_dependent_potassium_current__T = var_membrane__T;
   Real var_time_dependent_potassium_current__R = var_membrane__R;
   Real var_time_dependent_potassium_current__E_K = ((var_time_dependent_potassium_current__R * var_time_dependent_potassium_current__T) / var_time_dependent_potassium_current__F) * log((var_time_dependent_potassium_current__Ko + (var_time_dependent_potassium_current__PR_NaK * var_time_dependent_potassium_current__Nao)) / (var_time_dependent_potassium_current__Ki + (var_time_dependent_potassium_current__PR_NaK * var_time_dependent_potassium_current__Nai)));
   Real var_time_dependent_potassium_current__i_K = var_time_dependent_potassium_current__g_K * var_time_dependent_potassium_current__X * var_time_dependent_potassium_current__Xi * (var_time_dependent_potassium_current__V - var_time_dependent_potassium_current__E_K);
   Real var_membrane__i_K = var_time_dependent_potassium_current__i_K;
   Real var_time_independent_potassium_current__V = var_membrane__V;
   Real var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V;
   Real var_time_independent_potassium_current__Ki = var_ionic_concentrations__Ki;
   Real var_time_independent_potassium_current__R = var_membrane__R;
   Real var_time_independent_potassium_current__F = var_membrane__F;
   Real var_time_independent_potassium_current__Ko = var_ionic_concentrations__Ko;
   Real var_time_independent_potassium_current__T = var_membrane__T;
   Real var_time_independent_potassium_current__E_K1 = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__Ko / var_time_independent_potassium_current__Ki);
   Real var_time_independent_potassium_current_K1_gate__E_K1 = var_time_independent_potassium_current__E_K1;
   Real var_time_independent_potassium_current_K1_gate__beta_K1 = ((0.49124 * exp(0.08032 * ((var_time_independent_potassium_current_K1_gate__V + 5.476) - var_time_independent_potassium_current_K1_gate__E_K1))) + (1.0 * exp(0.06175 * (var_time_independent_potassium_current_K1_gate__V - (var_time_independent_potassium_current_K1_gate__E_K1 + 594.31))))) / (1.0 + exp((-0.5143) * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K1) + 4.753)));
   Real var_time_independent_potassium_current_K1_gate__alpha_K1 = 1.02 / (1.0 + exp(0.2385 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K1) - 59.215)));
   Real var_time_independent_potassium_current_K1_gate__K1_infinity = var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + var_time_independent_potassium_current_K1_gate__beta_K1);
   Real var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity;
   Real var_time_independent_potassium_current__g_K1 = 0.6047 * sqrt(var_time_independent_potassium_current__Ko / 5.4);
   Real var_time_independent_potassium_current__i_K1 = var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K1);
   Real var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
   const Real var_plateau_potassium_current__g_Kp = 0.0183;
   Real var_plateau_potassium_current__V = var_membrane__V;
   Real var_plateau_potassium_current__Kp = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current__V) / 5.98));
   Real var_plateau_potassium_current__E_K1 = var_time_independent_potassium_current__E_K1;
   Real var_plateau_potassium_current__E_Kp = var_plateau_potassium_current__E_K1;
   Real var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_Kp);
   Real var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
   const Real var_background_current__E_b =  -59.87;
   const Real var_background_current__g_b = 0.03921;
   Real var_background_current__V = var_membrane__V;
   Real var_background_current__i_b = var_background_current__g_b * (var_background_current__V - var_background_current__E_b);
   Real var_membrane__i_b = var_background_current__i_b;
   Real var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
   Real var_fast_sodium_current_m_gate__alpha_m = (var_fast_sodium_current_m_gate__V != -47.13) ? ((0.32 * (var_fast_sodium_current_m_gate__V + 47.13)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13)))) : ((0.32 * (var_fast_sodium_current_m_gate__V + 47.13 + 1.0e-06)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13 + 1.0e-06))));
   Real var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0);
   Real var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
   Real var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.135 * exp((80.0 + var_fast_sodium_current_h_gate__V) / (-6.8))) : 0.0;
   Real var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3.56 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1)))));
   Real var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
   Real var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? (((((-127140.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0;
   Real var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.1212 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.3 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))));
   Real var_slow_inward_current_d_gate__V = var_slow_inward_current__V;
   Real var_slow_inward_current_d_gate__alpha_d = (0.095 * exp((-0.01) * (var_slow_inward_current_d_gate__V - 5.0))) / (1.0 + exp((-0.072) * (var_slow_inward_current_d_gate__V - 5.0)));
   Real var_slow_inward_current_d_gate__beta_d = (0.07 * exp((-0.017) * (var_slow_inward_current_d_gate__V + 44.0))) / (1.0 + exp(0.05 * (var_slow_inward_current_d_gate__V + 44.0)));
   Real var_slow_inward_current_f_gate__V = var_slow_inward_current__V;
   Real var_slow_inward_current_f_gate__alpha_f = (0.012 * exp((-0.008) * (var_slow_inward_current_f_gate__V + 28.0))) / (1.0 + exp(0.15 * (var_slow_inward_current_f_gate__V + 28.0)));
   Real var_slow_inward_current_f_gate__beta_f = (0.0065 * exp((-0.02) * (var_slow_inward_current_f_gate__V + 30.0))) / (1.0 + exp((-0.2) * (var_slow_inward_current_f_gate__V + 30.0)));
   Real var_time_dependent_potassium_current_X_gate__V = var_time_dependent_potassium_current__V;
   Real var_time_dependent_potassium_current_X_gate__alpha_X = (0.0005 * exp(0.083 * (var_time_dependent_potassium_current_X_gate__V + 50.0))) / (1.0 + exp(0.057 * (var_time_dependent_potassium_current_X_gate__V + 50.0)));
   Real var_time_dependent_potassium_current_X_gate__beta_X = (0.0013 * exp((-0.06) * (var_time_dependent_potassium_current_X_gate__V + 20.0))) / (1.0 + exp((-0.04) * (var_time_dependent_potassium_current_X_gate__V + 20.0)));
   Real var_intracellular_calcium_concentration__i_si = var_slow_inward_current__i_si;
   Real d_dt_membrane__V = ((-1.0) / var_membrane__C) * (var_membrane__I_stim + var_membrane__i_Na + var_membrane__i_si + var_membrane__i_K + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_b);

   Real d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__alpha_m * (1.0 - var_fast_sodium_current_m_gate__m)) - (var_fast_sodium_current_m_gate__beta_m * var_fast_sodium_current_m_gate__m);
   Real d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - var_fast_sodium_current_h_gate__h)) - (var_fast_sodium_current_h_gate__beta_h * var_fast_sodium_current_h_gate__h);
   Real d_dt_fast_sodium_current_j_gate__j = (var_fast_sodium_current_j_gate__alpha_j * (1.0 - var_fast_sodium_current_j_gate__j)) - (var_fast_sodium_current_j_gate__beta_j * var_fast_sodium_current_j_gate__j);
   Real d_dt_slow_inward_current_d_gate__d = (var_slow_inward_current_d_gate__alpha_d * (1.0 - var_slow_inward_current_d_gate__d)) - (var_slow_inward_current_d_gate__beta_d * var_slow_inward_current_d_gate__d);
   Real d_dt_slow_inward_current_f_gate__f = (var_slow_inward_current_f_gate__alpha_f * (1.0 - var_slow_inward_current_f_gate__f)) - (var_slow_inward_current_f_gate__beta_f * var_slow_inward_current_f_gate__f);
   Real d_dt_time_dependent_potassium_current_X_gate__X = (var_time_dependent_potassium_current_X_gate__alpha_X * (1.0 - var_time_dependent_potassium_current_X_gate__X)) - (var_time_dependent_potassium_current_X_gate__beta_X * var_time_dependent_potassium_current_X_gate__X);
   Real d_dt_intracellular_calcium_concentration__Cai = (((-0.0001) / 1.0) * var_intracellular_calcium_concentration__i_si) + (0.07 * (0.0001 - var_intracellular_calcium_concentration__Cai));

   rDY[0] = d_dt_membrane__V;
   rDY[1] = d_dt_fast_sodium_current_m_gate__m;
   rDY[2] = d_dt_fast_sodium_current_h_gate__h;
   rDY[3] = d_dt_fast_sodium_current_j_gate__j;
   rDY[4] = d_dt_slow_inward_current_d_gate__d;
   rDY[5] = d_dt_slow_inward_current_f_gate__f;
   rDY[6] = d_dt_time_dependent_potassium_current_X_gate__X;
   rDY[7] = d_dt_intracellular_calcium_concentration__Cai;
}

