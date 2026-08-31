#include "luo_rudy.hpp"

LuoRudy::LuoRudy() : CellModel(8)
{
  var_names.insert( std::pair<int, std::string>(0, "V") );
  var_names.insert( std::pair<int, std::string>(1, "m") );
  var_names.insert( std::pair<int, std::string>(2, "h") );
  var_names.insert( std::pair<int, std::string>(3, "j") );
  var_names.insert( std::pair<int, std::string>(4, "d") );
  var_names.insert( std::pair<int, std::string>(5, "f") );
  var_names.insert( std::pair<int, std::string>(6, "X") );
  var_names.insert( std::pair<int, std::string>(7, "Cai") );
}

void LuoRudy::init(double * values) const
{
  assert(values != 0);

  // initial conditions
  values[0] = -83.853;    // V;   Units: millivolt;     Initial value: -83.853
  values[1] = 0.00187018; // m;   Units: dimensionless; Initial value: 0.00187018
  values[2] = 0.9804713;  // h;   Units: dimensionless; Initial value: 0.9804713
  values[3] = 0.98767124; // j;   Units: dimensionless; Initial value: 0.98767124
  values[4] = 0.00316354; // d;   Units: dimensionless; Initial value: 0.00316354
  values[5] = 0.99427859; // f;   Units: dimensionless; Initial value: 0.99427859
  values[6] = 0.16647703; // X;   Units: dimensionless; Initial value: 0.16647703
  values[7] = 0.0002;     // Cai; Units: millimolar;    Initial value: 0.0002
}

void LuoRudy::equation(const double time, const double * statevars, double * values)
{
  // State variables
  const double V = statevars[0]; // Membrane variable
  const double m = statevars[1]; // Fast sodium current m gate
  const double h = statevars[2]; // Fast sodium current h gate
  const double j = statevars[3]; // Fast sodium current j gate
  const double d = statevars[4]; // Slow inward current d gate
  const double f = statevars[5]; // Slow inward current f gate
  const double X = statevars[6]; // Time dependent potassium current X gate
  const double Cai = statevars[7]; // Intracellular calcium concentration

  // Inputs:
  // Time units: ms
  // double var_environment__time = time;

  // Some constants
  const double R = 8314.0;
  const double T = 310.0;
  const double F = 96484.6;
  const double C = 1.0;

  // Stimulus
  double var_membrane__I_stim = i_stim;

  // Calculations
  double var_fast_sodium_current__j = j;
  double var_fast_sodium_current__h = h;
  const double var_fast_sodium_current__g_Na = 23.0;
  double var_fast_sodium_current__m = m;
  double var_fast_sodium_current__V = V;
  double var_fast_sodium_current__R = R;
  double var_fast_sodium_current__F = F;
  const double var_ionic_concentrations__Nao = 140.0;
  double var_fast_sodium_current__Nao = var_ionic_concentrations__Nao;
  const double var_ionic_concentrations__Nai = 18.0;
  double var_fast_sodium_current__Nai = var_ionic_concentrations__Nai;
  double var_fast_sodium_current__T = T;
  double var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Nao / var_fast_sodium_current__Nai);
  double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
  double var_membrane__i_Na = var_fast_sodium_current__i_Na;
  double var_slow_inward_current__d = d;
  double var_slow_inward_current__f = f;
  double var_slow_inward_current__V = V;
  double var_slow_inward_current__Cai = Cai;
  double var_slow_inward_current__E_si = 7.7 - (13.0287 * log(var_slow_inward_current__Cai / 1.0));
  double var_slow_inward_current__i_si = 0.09 * var_slow_inward_current__d * var_slow_inward_current__f * (var_slow_inward_current__V - var_slow_inward_current__E_si);
  double var_membrane__i_si = var_slow_inward_current__i_si;
  double var_time_dependent_potassium_current__V = V;
  double var_time_dependent_potassium_current_Xi_gate__V = var_time_dependent_potassium_current__V;
  double var_time_dependent_potassium_current_Xi_gate__Xi = (var_time_dependent_potassium_current_Xi_gate__V > (-100.0)) ? ((var_time_dependent_potassium_current_Xi_gate__V != -77.0)?((2.837 * (exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 77.0)) - 1.0)) / ((var_time_dependent_potassium_current_Xi_gate__V + 77.0) * exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 35.0)))):((2.837 * (exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 77.0+1.0e-6)) - 1.0)) / ((var_time_dependent_potassium_current_Xi_gate__V + 77.0+1.0e-6) * exp(0.04 * (var_time_dependent_potassium_current_Xi_gate__V + 35.0)))) ): 1.0;
  double var_time_dependent_potassium_current__Xi = var_time_dependent_potassium_current_Xi_gate__Xi;
  const double var_ionic_concentrations__Ko = 5.4;
  double var_time_dependent_potassium_current__Ko = var_ionic_concentrations__Ko;
  double var_time_dependent_potassium_current__g_K = 0.282 * sqrt(var_time_dependent_potassium_current__Ko / 5.4);
  double var_time_dependent_potassium_current__X = X;
  const double var_time_dependent_potassium_current__PR_NaK = 0.01833;
  double var_time_dependent_potassium_current__F = F;
  double var_time_dependent_potassium_current__Nao = var_ionic_concentrations__Nao;
  const double var_ionic_concentrations__Ki = 145.0;
  double var_time_dependent_potassium_current__Ki = var_ionic_concentrations__Ki;
  double var_time_dependent_potassium_current__Nai = var_ionic_concentrations__Nai;
  double var_time_dependent_potassium_current__T = T;
  double var_time_dependent_potassium_current__R = R;
  double var_time_dependent_potassium_current__E_K = ((var_time_dependent_potassium_current__R * var_time_dependent_potassium_current__T) / var_time_dependent_potassium_current__F) * log((var_time_dependent_potassium_current__Ko + (var_time_dependent_potassium_current__PR_NaK * var_time_dependent_potassium_current__Nao)) / (var_time_dependent_potassium_current__Ki + (var_time_dependent_potassium_current__PR_NaK * var_time_dependent_potassium_current__Nai)));
  double var_time_dependent_potassium_current__i_K = var_time_dependent_potassium_current__g_K * var_time_dependent_potassium_current__X * var_time_dependent_potassium_current__Xi * (var_time_dependent_potassium_current__V - var_time_dependent_potassium_current__E_K);
  double var_membrane__i_K = var_time_dependent_potassium_current__i_K;
  double var_time_independent_potassium_current__V = V;
  double var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V;
  double var_time_independent_potassium_current__Ki = var_ionic_concentrations__Ki;
  double var_time_independent_potassium_current__R = R;
  double var_time_independent_potassium_current__F = F;
  double var_time_independent_potassium_current__Ko = var_ionic_concentrations__Ko;
  double var_time_independent_potassium_current__T = T;
  double var_time_independent_potassium_current__E_K1 = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__Ko / var_time_independent_potassium_current__Ki);
  double var_time_independent_potassium_current_K1_gate__E_K1 = var_time_independent_potassium_current__E_K1;
  double var_time_independent_potassium_current_K1_gate__beta_K1 = ((0.49124 * exp(0.08032 * ((var_time_independent_potassium_current_K1_gate__V + 5.476) - var_time_independent_potassium_current_K1_gate__E_K1))) + (1.0 * exp(0.06175 * (var_time_independent_potassium_current_K1_gate__V - (var_time_independent_potassium_current_K1_gate__E_K1 + 594.31))))) / (1.0 + exp((-0.5143) * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K1) + 4.753)));
  double var_time_independent_potassium_current_K1_gate__alpha_K1 = 1.02 / (1.0 + exp(0.2385 * ((var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K1) - 59.215)));
  double var_time_independent_potassium_current_K1_gate__K1_infinity = var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + var_time_independent_potassium_current_K1_gate__beta_K1);
  double var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity;
  double var_time_independent_potassium_current__g_K1 = 0.6047 * sqrt(var_time_independent_potassium_current__Ko / 5.4);
  double var_time_independent_potassium_current__i_K1 = var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K1);
  double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
  const double var_plateau_potassium_current__g_Kp = 0.0183;
  double var_plateau_potassium_current__V = V;
  double var_plateau_potassium_current__Kp = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current__V) / 5.98));
  double var_plateau_potassium_current__E_K1 = var_time_independent_potassium_current__E_K1;
  double var_plateau_potassium_current__E_Kp = var_plateau_potassium_current__E_K1;
  double var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_Kp);
  double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
  const double var_background_current__E_b =  -59.87;
  const double var_background_current__g_b = 0.03921;
  double var_background_current__V = V;
  double var_background_current__i_b = var_background_current__g_b * (var_background_current__V - var_background_current__E_b);
  double var_membrane__i_b = var_background_current__i_b;
  double var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
  double var_fast_sodium_current_m_gate__alpha_m = (var_fast_sodium_current_m_gate__V != -47.13) ? ((0.32 * (var_fast_sodium_current_m_gate__V + 47.13)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13)))) : ((0.32 * (var_fast_sodium_current_m_gate__V + 47.13 + 1.0e-06)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13 + 1.0e-06))));
  double var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0);
  double var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
  double var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.135 * exp((80.0 + var_fast_sodium_current_h_gate__V) / (-6.8))) : 0.0;
  double var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3.56 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1)))));
  double var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
  double var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? (((((-127140.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0;
  double var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.1212 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.3 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))));
  double var_slow_inward_current_d_gate__V = var_slow_inward_current__V;
  double var_slow_inward_current_d_gate__alpha_d = (0.095 * exp((-0.01) * (var_slow_inward_current_d_gate__V - 5.0))) / (1.0 + exp((-0.072) * (var_slow_inward_current_d_gate__V - 5.0)));
  double var_slow_inward_current_d_gate__beta_d = (0.07 * exp((-0.017) * (var_slow_inward_current_d_gate__V + 44.0))) / (1.0 + exp(0.05 * (var_slow_inward_current_d_gate__V + 44.0)));
  double var_slow_inward_current_f_gate__V = var_slow_inward_current__V;
  double var_slow_inward_current_f_gate__alpha_f = (0.012 * exp((-0.008) * (var_slow_inward_current_f_gate__V + 28.0))) / (1.0 + exp(0.15 * (var_slow_inward_current_f_gate__V + 28.0)));
  double var_slow_inward_current_f_gate__beta_f = (0.0065 * exp((-0.02) * (var_slow_inward_current_f_gate__V + 30.0))) / (1.0 + exp((-0.2) * (var_slow_inward_current_f_gate__V + 30.0)));
  double var_time_dependent_potassium_current_X_gate__V = var_time_dependent_potassium_current__V;
  double var_time_dependent_potassium_current_X_gate__alpha_X = (0.0005 * exp(0.083 * (var_time_dependent_potassium_current_X_gate__V + 50.0))) / (1.0 + exp(0.057 * (var_time_dependent_potassium_current_X_gate__V + 50.0)));
  double var_time_dependent_potassium_current_X_gate__beta_X = (0.0013 * exp((-0.06) * (var_time_dependent_potassium_current_X_gate__V + 20.0))) / (1.0 + exp((-0.04) * (var_time_dependent_potassium_current_X_gate__V + 20.0)));
  double var_intracellular_calcium_concentration__i_si = var_slow_inward_current__i_si;
  double d_dt_membrane__V = ((-1.0) / C) * (var_membrane__I_stim + var_membrane__i_Na + var_membrane__i_si + var_membrane__i_K + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_b);

  double d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__alpha_m * (1.0 - m)) - (var_fast_sodium_current_m_gate__beta_m * m);
  double d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - h)) - (var_fast_sodium_current_h_gate__beta_h * h);
  double d_dt_fast_sodium_current_j_gate__j = (var_fast_sodium_current_j_gate__alpha_j * (1.0 - j)) - (var_fast_sodium_current_j_gate__beta_j * j);
  double d_dt_slow_inward_current_d_gate__d = (var_slow_inward_current_d_gate__alpha_d * (1.0 - d)) - (var_slow_inward_current_d_gate__beta_d * d);
  double d_dt_slow_inward_current_f_gate__f = (var_slow_inward_current_f_gate__alpha_f * (1.0 - f)) - (var_slow_inward_current_f_gate__beta_f * f);
  double d_dt_time_dependent_potassium_current_X_gate__X = (var_time_dependent_potassium_current_X_gate__alpha_X * (1.0 - X)) - (var_time_dependent_potassium_current_X_gate__beta_X * X);
  double d_dt_intracellular_calcium_concentration__Cai = (((-0.0001) / 1.0) * var_intracellular_calcium_concentration__i_si) + (0.07 * (0.0001 - Cai));

  values[0] = d_dt_membrane__V;
  values[1] = d_dt_fast_sodium_current_m_gate__m;
  values[2] = d_dt_fast_sodium_current_h_gate__h;
  values[3] = d_dt_fast_sodium_current_j_gate__j;
  values[4] = d_dt_slow_inward_current_d_gate__d;
  values[5] = d_dt_slow_inward_current_f_gate__f;
  values[6] = d_dt_time_dependent_potassium_current_X_gate__X;
  values[7] = d_dt_intracellular_calcium_concentration__Cai;
}

