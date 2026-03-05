#include "courtemanche.h"

void setIC_ode_cpu(Real *sv, int cellID, int NEQ = 20)
{
   // Vector containing the initial conditions of the cell
   Real IC[NEQ];
 
   IC[0]  = -81.18;
   IC[1]  = 2.908e-3;
   IC[2]  = 9.649e-1;
   IC[3]  = 9.775e-1;
   IC[4]  = 3.043e-2;
   IC[5]  = 9.992e-1;
   IC[6]  = 4.966e-3;
   IC[7]  = 9.986e-1;
   IC[8]  = 3.296e-5;
   IC[9]  = 1.869e-2;

   IC[10] = 1.367e-4;
   IC[11] = 9.996e-1;
   IC[12] = 7.755e-1;
   IC[13] = 2.35e-112;
   IC[14] = 1;
   IC[15] = 0.9992;
   IC[16] = 1.117e+01;
   IC[17] = 1.013e-4;
   IC[18] = 1.39e+02;
   IC[19] = 1.488;

   IC[20] = 1.488;
   //IC[21] = ;

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = IC[i];
}

void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, int NEQ, int col)
{
   solve_Forward_Euler_cpu(time, sv, dt, NEQ, cellID, col);
}

void solve_Forward_Euler_cpu(Real time, Real *sv, Real dt, int NEQ, int cellID, int col)
{
   Real rY[NEQ], rDY[NEQ];

   for(int i = 0; i < NEQ; i++)
      rY[i] = sv[cellID+i];

   RHS_Courtemanche_cpu(time, rY, rDY, col);

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = dt*rDY[i] + rY[i];
}

void RHS_Courtemanche_cpu(Real time, Real *rY, Real *rDY, int col)
{
  Real var_environment__time = time;


  // Time units: millisecond
  Real var_membrane__V = rY[0];
  // Units: millivolt; Initial value: -81.18
  Real var_fast_sodium_current_m_gate__m = rY[1];
  // Units: dimensionless; Initial value: 2.908e-3
  Real var_fast_sodium_current_h_gate__h = rY[2];
  // Units: dimensionless; Initial value: 9.649e-1
  Real var_fast_sodium_current_j_gate__j = rY[3];
  // Units: dimensionless; Initial value: 9.775e-1
  Real var_transient_outward_K_current_oa_gate__oa = rY[4];
  // Units: dimensionless; Initial value: 3.043e-2
  Real var_transient_outward_K_current_oi_gate__oi = rY[5];
  // Units: dimensionless; Initial value: 9.992e-1
  Real var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua = rY[6];
  // Units: dimensionless; Initial value: 4.966e-3
  Real var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui = rY[7];
  // Units: dimensionless; Initial value: 9.986e-1
  Real var_rapid_delayed_rectifier_K_current_xr_gate__xr = rY[8];
  // Units: dimensionless; Initial value: 3.296e-5
  Real var_slow_delayed_rectifier_K_current_xs_gate__xs = rY[9];
  // Units: dimensionless; Initial value: 1.869e-2
  Real var_L_type_Ca_channel_d_gate__d = rY[10];
  // Units: dimensionless; Initial value: 1.367e-4
  Real var_L_type_Ca_channel_f_gate__f = rY[11];
  // Units: dimensionless; Initial value: 9.996e-1
  Real var_L_type_Ca_channel_f_Ca_gate__f_Ca = rY[12];
  // Units: dimensionless; Initial value: 7.755e-1
  Real var_Ca_release_current_from_JSR_u_gate__u = rY[13];
  // Units: dimensionless; Initial value: 2.35e-112
  Real var_Ca_release_current_from_JSR_v_gate__v = rY[14];
  // Units: dimensionless; Initial value: 1
  Real var_Ca_release_current_from_JSR_w_gate__w = rY[15];
  // Units: dimensionless; Initial value: 0.9992
  Real var_intracellular_ion_concentrations__Na_i = rY[16];
  // Units: millimolar; Initial value: 1.117e+01
  Real var_intracellular_ion_concentrations__Ca_i = rY[17];
  // Units: millimolar; Initial value: 1.013e-4
  Real var_intracellular_ion_concentrations__K_i = rY[18];
  // Units: millimolar; Initial value: 1.39e+02
  Real var_intracellular_ion_concentrations__Ca_rel = rY[19];
  // Units: millimolar; Initial value: 1.488
  Real var_intracellular_ion_concentrations__Ca_up = rY[20];
  // Units: millimolar; Initial value: 1.488
        
  // Mathematics
  const Real var_fast_sodium_current__E_Na = 26.712831924 * log(140.0 / var_intracellular_ion_concentrations__Na_i);
  const Real var_fast_sodium_current__i_Na = 100.0 * 7.8 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * var_fast_sodium_current_j_gate__j * (var_membrane__V - var_fast_sodium_current__E_Na);
  const Real var_time_independent_potassium_current__E_K = 26.712831924 * log(5.4 / var_intracellular_ion_concentrations__K_i);
  const Real var_time_independent_potassium_current__i_K1 = (100.0 * 0.09 * (var_membrane__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp(0.07 * (var_membrane__V + 80.0)));
  const Real var_transient_outward_K_current__i_to = 100.0 * 0.1652 * pow(var_transient_outward_K_current_oa_gate__oa, 3.0) * var_transient_outward_K_current_oi_gate__oi * (var_membrane__V - var_time_independent_potassium_current__E_K);
  const Real var_ultrarapid_delayed_rectifier_K_current__i_Kur = 100.0 * (0.005 + (0.05 / (1.0 + exp((var_membrane__V - 15.0) *  -0.0769230769231)))) * pow(var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua, 3.0) * var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui * (var_membrane__V - var_time_independent_potassium_current__E_K);
  const Real var_rapid_delayed_rectifier_K_current__i_Kr = (100.0 * 0.029411765 * var_rapid_delayed_rectifier_K_current_xr_gate__xr * (var_membrane__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp((var_membrane__V + 15.0) * 0.0446428571429));
  const Real var_slow_delayed_rectifier_K_current__i_Ks = 100.0 * 0.12941176 * pow(var_slow_delayed_rectifier_K_current_xs_gate__xs, 2.0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
  const Real var_L_type_Ca_channel__i_Ca_L = 100.0 * 0.12375 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * (var_membrane__V - 65.0);
  const Real var_sarcolemmal_calcium_pump_current__i_CaP = (100.0 * 0.275 * var_intracellular_ion_concentrations__Ca_i) / (0.0005 + var_intracellular_ion_concentrations__Ca_i);
  const Real var_sodium_potassium_pump__i_NaK = (((100.0 * 0.59933874 * pow(1.0 + (0.1245 * exp(( -0.1 * 96.4867 * var_membrane__V) * 0.000387982927199)) + (0.0365 * 1.00091030495 * exp(( -96.4867 * var_membrane__V) * 0.000387982927199)),  -1.0) * 1.0) / (1.0 + pow(10.0 / var_intracellular_ion_concentrations__Na_i, 1.5))) * 5.4) * 0.144927536232;
  const Real var_Na_Ca_exchanger_current__i_NaCa = (100.0 * 1600.0 * ((exp((0.35 * 96.4867 * var_membrane__V) * 0.000387982927199) * pow(var_intracellular_ion_concentrations__Na_i, 3.0) * 1.8) - (exp(( -0.65 * 96.4867 * var_membrane__V) * 0.000387982927199) * 2744000.0 * var_intracellular_ion_concentrations__Ca_i))) / (3413921.875 * 3.18 * (1.0 + (0.1 * exp(( -0.65 * var_membrane__V * 96.4867) * 0.000387982927199))));
  const Real var_background_currents__i_B_Na = 100.0 * 0.0006744375 * (var_membrane__V - var_fast_sodium_current__E_Na);
  const Real var_background_currents__i_B_Ca = 100.0 * 0.001131 * (var_membrane__V - (13.356415962 * log(1.8 / var_intracellular_ion_concentrations__Ca_i)));
  const Real var_fast_sodium_current_m_gate__alpha_m = (var_membrane__V ==  -47.13) ? 3.2 : ((0.32 * (var_membrane__V + 47.13)) / (1.0 - exp( -0.1 * (var_membrane__V + 47.13))));
  const Real var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_membrane__V) * 0.0909090909091);
  const Real var_fast_sodium_current_h_gate__alpha_h = (var_membrane__V <  -40.0) ? (0.135 * exp((var_membrane__V + 80.0) *  -0.147058823529)) : 0.0;
  const Real var_fast_sodium_current_h_gate__beta_h = (var_membrane__V <  -40.0) ? ((3.56 * exp(0.079 * var_membrane__V)) + (310000.0 * exp(0.35 * var_membrane__V))) : (1.0 / (0.13 * (1.0 + exp((var_membrane__V + 10.66) *  -0.0900900900901))));
  const Real var_fast_sodium_current_j_gate__alpha_j = (var_membrane__V <  -40.0) ? (((( -127140.0 * exp(0.2444 * var_membrane__V)) - (3.474e-05 * exp( -0.04391 * var_membrane__V))) * (var_membrane__V + 37.78)) / (1.0 + exp(0.311 * (var_membrane__V + 79.23)))) : 0.0;
  const Real var_fast_sodium_current_j_gate__beta_j = (var_membrane__V <  -40.0) ? ((0.1212 * exp( -0.01052 * var_membrane__V)) / (1.0 + exp( -0.1378 * (var_membrane__V + 40.14)))) : ((0.3 * exp( -2.535e-07 * var_membrane__V)) / (1.0 + exp( -0.1 * (var_membrane__V + 32.0))));
  const Real var_Ca_release_current_from_JSR__i_rel = 30.0 * pow(var_Ca_release_current_from_JSR_u_gate__u, 2.0) * var_Ca_release_current_from_JSR_v_gate__v * var_Ca_release_current_from_JSR_w_gate__w * (var_intracellular_ion_concentrations__Ca_rel - var_intracellular_ion_concentrations__Ca_i);
  const Real var_Ca_release_current_from_JSR__Fn = 1000.0 * ((1e-15 * 96.48 * var_Ca_release_current_from_JSR__i_rel) - (5.18206136182e-18 * ((0.5 * var_L_type_Ca_channel__i_Ca_L) - (0.2 * var_Na_Ca_exchanger_current__i_NaCa))));
  const Real var_transfer_current_from_NSR_to_JSR__i_tr = (var_intracellular_ion_concentrations__Ca_up - var_intracellular_ion_concentrations__Ca_rel) * 0.00555555555556;
  const Real var_Ca_uptake_current_by_the_NSR__i_up = 0.005 / (1.0 + (0.00092 / var_intracellular_ion_concentrations__Ca_i));
  const Real var_Ca_leak_current_by_the_NSR__i_up_leak = (0.005 * var_intracellular_ion_concentrations__Ca_up) * 0.0666666666667;
  const Real d_dt_membrane__V = (-(var_fast_sodium_current__i_Na + var_time_independent_potassium_current__i_K1 + var_transient_outward_K_current__i_to + var_ultrarapid_delayed_rectifier_K_current__i_Kur + var_rapid_delayed_rectifier_K_current__i_Kr + var_slow_delayed_rectifier_K_current__i_Ks + var_background_currents__i_B_Na + var_background_currents__i_B_Ca + var_sodium_potassium_pump__i_NaK + var_sarcolemmal_calcium_pump_current__i_CaP + var_Na_Ca_exchanger_current__i_NaCa + var_L_type_Ca_channel__i_Ca_L + (((var_environment__time >= 50.0) && (var_environment__time <= 50000.0) && (((var_environment__time - 50.0) - (floor((var_environment__time - 50.0) * 0.001) * 1000.0)) <= 2.0)) ?  -2000.0 : 0.0))) * 0.01;
  const Real d_dt_fast_sodium_current_m_gate__m = ((var_fast_sodium_current_m_gate__alpha_m / (var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m)) - var_fast_sodium_current_m_gate__m) / (1.0 / (var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m));
  const Real d_dt_fast_sodium_current_h_gate__h = ((var_fast_sodium_current_h_gate__alpha_h / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h)) - var_fast_sodium_current_h_gate__h) / (1.0 / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h));
  const Real d_dt_fast_sodium_current_j_gate__j = ((var_fast_sodium_current_j_gate__alpha_j / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j)) - var_fast_sodium_current_j_gate__j) / (1.0 / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j));
  const Real d_dt_transient_outward_K_current_oa_gate__oa = (pow(1.0 + exp(((var_membrane__V -  -10.0) + 10.47) *  -0.0570125427594),  -1.0) - var_transient_outward_K_current_oa_gate__oa) / (pow((0.65 * pow(exp((var_membrane__V -  -10.0) *  -0.117647058824) + exp(((var_membrane__V -  -10.0) - 40.0) *  -0.0169491525424),  -1.0)) + (0.65 * pow(2.5 + exp(((var_membrane__V -  -10.0) + 72.0) * 0.0588235294118),  -1.0)),  -1.0) * 0.333333333333);
  const Real d_dt_transient_outward_K_current_oi_gate__oi = (pow(1.0 + exp(((var_membrane__V -  -10.0) + 33.1) * 0.188679245283),  -1.0) - var_transient_outward_K_current_oi_gate__oi) / (pow(pow(18.53 + (1.0 * exp(((var_membrane__V -  -10.0) + 103.7) * 0.0913242009132)),  -1.0) + pow(35.56 + (1.0 * exp(((var_membrane__V -  -10.0) - 8.74) *  -0.134408602151)),  -1.0),  -1.0) * 0.333333333333);
  const Real d_dt_ultrarapid_delayed_rectifier_K_current_ua_gate__ua = (pow(1.0 + exp(((var_membrane__V -  -10.0) + 20.3) *  -0.104166666667),  -1.0) - var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua) / (pow((0.65 * pow(exp((var_membrane__V -  -10.0) *  -0.117647058824) + exp(((var_membrane__V -  -10.0) - 40.0) *  -0.0169491525424),  -1.0)) + (0.65 * pow(2.5 + exp(((var_membrane__V -  -10.0) + 72.0) * 0.0588235294118),  -1.0)),  -1.0) * 0.333333333333);
  const Real d_dt_ultrarapid_delayed_rectifier_K_current_ui_gate__ui = (pow(1.0 + exp(((var_membrane__V -  -10.0) - 109.45) * 0.0363901018923),  -1.0) - var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui) / (pow(pow(21.0 + (1.0 * exp(((var_membrane__V -  -10.0) - 195.0) *  -0.0357142857143)),  -1.0) + (1.0 / exp(((var_membrane__V -  -10.0) - 168.0) *  -0.0625)),  -1.0) * 0.333333333333);
  const Real d_dt_rapid_delayed_rectifier_K_current_xr_gate__xr = (pow(1.0 + exp((var_membrane__V + 14.1) *  -0.153846153846),  -1.0) - var_rapid_delayed_rectifier_K_current_xr_gate__xr) / pow(((fabs(var_membrane__V + 14.1) < 1e-10) ? 0.0015 : ((0.0003 * (var_membrane__V + 14.1)) / (1.0 - exp((var_membrane__V + 14.1) *  -0.2)))) + ((fabs(var_membrane__V - 3.3328) < 1e-10) ? 0.00037836118 : ((7.3898e-05 * (var_membrane__V - 3.3328)) / (exp((var_membrane__V - 3.3328) * 0.195171458126) - 1.0))),  -1.0);
  const Real d_dt_slow_delayed_rectifier_K_current_xs_gate__xs = (pow(1.0 + exp((var_membrane__V - 19.9) *  -0.0787401574803),  -0.5) - var_slow_delayed_rectifier_K_current_xs_gate__xs) / (0.5 * pow(((fabs(var_membrane__V - 19.9) < 1e-10) ? 0.00068 : ((4e-05 * (var_membrane__V - 19.9)) / (1.0 - exp((var_membrane__V - 19.9) *  -0.0588235294118)))) + ((fabs(var_membrane__V - 19.9) < 1e-10) ? 0.000315 : ((3.5e-05 * (var_membrane__V - 19.9)) / (exp((var_membrane__V - 19.9) * 0.111111111111) - 1.0))),  -1.0));
  const Real d_dt_L_type_Ca_channel_d_gate__d = (pow(1.0 + exp((var_membrane__V + 10.0) *  -0.125),  -1.0) - var_L_type_Ca_channel_d_gate__d) / ((fabs(var_membrane__V + 10.0) < 1e-10) ? (4.579 / (1.0 + exp((var_membrane__V + 10.0) *  -0.160256410256))) : ((1.0 - exp((var_membrane__V + 10.0) *  -0.160256410256)) / (0.035 * (var_membrane__V + 10.0) * (1.0 + exp((var_membrane__V + 10.0) *  -0.160256410256)))));
  const Real d_dt_L_type_Ca_channel_f_gate__f = ((exp((-(var_membrane__V + 28.0)) * 0.144927536232) / (1.0 + exp((-(var_membrane__V + 28.0)) * 0.144927536232))) - var_L_type_Ca_channel_f_gate__f) / (9.0 * pow((0.0197 * exp( -0.00113569 * pow(var_membrane__V + 10.0, 2.0))) + 0.02,  -1.0));
  const Real d_dt_L_type_Ca_channel_f_Ca_gate__f_Ca = (pow(1.0 + (var_intracellular_ion_concentrations__Ca_i * 2857.14285714),  -1.0) - var_L_type_Ca_channel_f_Ca_gate__f_Ca) * 0.5;
  const Real d_dt_Ca_release_current_from_JSR_u_gate__u = (pow(1.0 + exp((-(var_Ca_release_current_from_JSR__Fn - 3.4175e-13)) * 7.31528895391e+14),  -1.0) - var_Ca_release_current_from_JSR_u_gate__u) * 0.125;
  const Real d_dt_Ca_release_current_from_JSR_v_gate__v = ((1.0 - pow(1.0 + exp((-(var_Ca_release_current_from_JSR__Fn - 6.835e-14)) * 7.31528895391e+14),  -1.0)) - var_Ca_release_current_from_JSR_v_gate__v) / (1.91 + (2.09 * pow(1.0 + exp((-(var_Ca_release_current_from_JSR__Fn - 3.4175e-13)) * 7.31528895391e+14),  -1.0)));
  const Real d_dt_Ca_release_current_from_JSR_w_gate__w = ((1.0 - pow(1.0 + exp((-(var_membrane__V - 40.0)) * 0.0588235294118),  -1.0)) - var_Ca_release_current_from_JSR_w_gate__w) / ((fabs(var_membrane__V - 7.9) < 1e-10) ? 0.923076923077 : ((6.0 * (1.0 - exp((-(var_membrane__V - 7.9)) * 0.2))) / ((1.0 + (0.3 * exp((-(var_membrane__V - 7.9)) * 0.2))) * 1.0 * (var_membrane__V - 7.9))));
  const Real d_dt_intracellular_ion_concentrations__Na_i = (( -3.0 * var_sodium_potassium_pump__i_NaK) - ((3.0 * var_Na_Ca_exchanger_current__i_NaCa) + var_background_currents__i_B_Na + var_fast_sodium_current__i_Na)) * 7.58276465002e-07;
  const Real d_dt_intracellular_ion_concentrations__K_i = ((2.0 * var_sodium_potassium_pump__i_NaK) - (var_time_independent_potassium_current__i_K1 + var_transient_outward_K_current__i_to + var_ultrarapid_delayed_rectifier_K_current__i_Kur + var_rapid_delayed_rectifier_K_current__i_Kr + var_slow_delayed_rectifier_K_current__i_Ks + (100.0 * 0.0 * (var_membrane__V - var_time_independent_potassium_current__E_K)))) * 7.58276465002e-07;
  const Real d_dt_intracellular_ion_concentrations__Ca_i = ((((2.0 * var_Na_Ca_exchanger_current__i_NaCa) - (var_sarcolemmal_calcium_pump_current__i_CaP + var_L_type_Ca_channel__i_Ca_L + var_background_currents__i_B_Ca)) * 3.79138232501e-07) + (((1109.52 * (var_Ca_leak_current_by_the_NSR__i_up_leak - var_Ca_uptake_current_by_the_NSR__i_up)) + (var_Ca_release_current_from_JSR__i_rel * 96.48)) * 7.31635937957e-05)) / (1.0 + (3.5e-05 / pow(var_intracellular_ion_concentrations__Ca_i + 0.0005, 2.0)) + (0.000119 / pow(var_intracellular_ion_concentrations__Ca_i + 0.00238, 2.0)));
  const Real d_dt_intracellular_ion_concentrations__Ca_up = var_Ca_uptake_current_by_the_NSR__i_up - (var_Ca_leak_current_by_the_NSR__i_up_leak + ((var_transfer_current_from_NSR_to_JSR__i_tr * 96.48) * 0.000901290648208));
  const Real d_dt_intracellular_ion_concentrations__Ca_rel = (var_transfer_current_from_NSR_to_JSR__i_tr - var_Ca_release_current_from_JSR__i_rel) * pow(1.0 + (8.0 / pow(var_intracellular_ion_concentrations__Ca_rel + 0.8, 2.0)),  -1.0);
        

  rDY[0] = d_dt_membrane__V;
  rDY[1] = d_dt_fast_sodium_current_m_gate__m;
  rDY[2] = d_dt_fast_sodium_current_h_gate__h;
  rDY[3] = d_dt_fast_sodium_current_j_gate__j;
  rDY[4] = d_dt_transient_outward_K_current_oa_gate__oa;
  rDY[5] = d_dt_transient_outward_K_current_oi_gate__oi;
  rDY[6] = d_dt_ultrarapid_delayed_rectifier_K_current_ua_gate__ua;
  rDY[7] = d_dt_ultrarapid_delayed_rectifier_K_current_ui_gate__ui;
  rDY[8] = d_dt_rapid_delayed_rectifier_K_current_xr_gate__xr;
  rDY[9] = d_dt_slow_delayed_rectifier_K_current_xs_gate__xs;
  rDY[10] = d_dt_L_type_Ca_channel_d_gate__d;
  rDY[11] = d_dt_L_type_Ca_channel_f_gate__f;
  rDY[12] = d_dt_L_type_Ca_channel_f_Ca_gate__f_Ca;
  rDY[13] = d_dt_Ca_release_current_from_JSR_u_gate__u;
  rDY[14] = d_dt_Ca_release_current_from_JSR_v_gate__v;
  rDY[15] = d_dt_Ca_release_current_from_JSR_w_gate__w;
  rDY[16] = d_dt_intracellular_ion_concentrations__Na_i;
  rDY[17] = d_dt_intracellular_ion_concentrations__Ca_i;
  rDY[18] = d_dt_intracellular_ion_concentrations__K_i;
  rDY[19] = d_dt_intracellular_ion_concentrations__Ca_rel;
  rDY[20] = d_dt_intracellular_ion_concentrations__Ca_up;


  /*
  // Inputs:
  // Time units: millisecond
  double var_membrane__V = rY[0];
  // Units: millivolt; Initial value: -81.18
  double var_fast_sodium_current_m_gate__m = rY[1];
  // Units: dimensionless; Initial value: 2.908e-3
  double var_fast_sodium_current_h_gate__h = rY[2];
  // Units: dimensionless; Initial value: 9.649e-1
  double var_fast_sodium_current_j_gate__j = rY[3];
  // Units: dimensionless; Initial value: 9.775e-1
  double var_transient_outward_K_current_oa_gate__oa = rY[4];
  // Units: dimensionless; Initial value: 3.043e-2
  double var_transient_outward_K_current_oi_gate__oi = rY[5];
  // Units: dimensionless; Initial value: 9.992e-1
  double var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua = rY[6];
  // Units: dimensionless; Initial value: 4.966e-3
  double var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui = rY[7];
  // Units: dimensionless; Initial value: 9.986e-1
  double var_rapid_delayed_rectifier_K_current_xr_gate__xr = rY[8];
  // Units: dimensionless; Initial value: 3.296e-5
  double var_slow_delayed_rectifier_K_current_xs_gate__xs = rY[9];
  // Units: dimensionless; Initial value: 1.869e-2
  double var_L_type_Ca_channel_d_gate__d = rY[10];
  // Units: dimensionless; Initial value: 1.367e-4
  double var_L_type_Ca_channel_f_gate__f = rY[11];
  // Units: dimensionless; Initial value: 9.996e-1
  double var_L_type_Ca_channel_f_Ca_gate__f_Ca = rY[12];
  // Units: dimensionless; Initial value: 7.755e-1
  double var_Ca_release_current_from_JSR_u_gate__u = rY[13];
  // Units: dimensionless; Initial value: 2.35e-112
  double var_Ca_release_current_from_JSR_v_gate__v = rY[14];
  // Units: dimensionless; Initial value: 1
  double var_Ca_release_current_from_JSR_w_gate__w = rY[15];
  // Units: dimensionless; Initial value: 0.9992
  double var_intracellular_ion_concentrations__Na_i = rY[16];
  // Units: millimolar; Initial value: 1.117e+01
  double var_intracellular_ion_concentrations__Ca_i = rY[17];
  // Units: millimolar; Initial value: 1.013e-4
  double var_intracellular_ion_concentrations__K_i = rY[18];
  // Units: millimolar; Initial value: 1.39e+02
  double var_intracellular_ion_concentrations__Ca_rel = rY[19];
  // Units: millimolar; Initial value: 1.488
  double var_intracellular_ion_concentrations__Ca_up = rY[20];
  // Units: millimolar; Initial value: 1.488
        
  // Mathematics

	

  const double var_membrane__R = 8.3143;
  const double var_membrane__T = 310.0;
  const double var_membrane__F = 96.4867;
  const double var_membrane__Cm = 100.0;
  const double var_membrane__stim_end = 50000.0;
  const double var_membrane__stim_amplitude =  -2000.0;
  const double var_membrane__time = var_environment__time;
  const double var_membrane__stim_duration = 2.0;
  const double var_membrane__stim_period = 1000.0;
  const double var_membrane__stim_start = 50.0;
  const double var_membrane__i_st = ((var_membrane__time >= var_membrane__stim_start) && (var_membrane__time <= var_membrane__stim_end) && (((var_membrane__time - var_membrane__stim_start) - (floor((var_membrane__time - var_membrane__stim_start) / var_membrane__stim_period) * var_membrane__stim_period)) <= var_membrane__stim_duration)) ? var_membrane__stim_amplitude : 0.0;
  const double var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
  const double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
  const double var_fast_sodium_current__g_Na = 7.8;
  const double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
  const double var_fast_sodium_current__Cm = var_membrane__Cm;
  const double var_fast_sodium_current__V = var_membrane__V;
  const double var_fast_sodium_current__R = var_membrane__R;
  const double var_standard_ionic_concentrations__Na_o = 140.0;
  const double var_fast_sodium_current__Na_o = var_standard_ionic_concentrations__Na_o;
  const double var_fast_sodium_current__F = var_membrane__F;
  const double var_fast_sodium_current__T = var_membrane__T;
  const double var_fast_sodium_current__Na_i = var_intracellular_ion_concentrations__Na_i;
  const double var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Na_o / var_fast_sodium_current__Na_i);
  const double var_fast_sodium_current__i_Na = var_fast_sodium_current__Cm * var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
  const double var_membrane__i_Na = var_fast_sodium_current__i_Na;
  const double var_time_independent_potassium_current__g_K1 = 0.09;
  const double var_time_independent_potassium_current__Cm = var_membrane__Cm;
  const double var_time_independent_potassium_current__V = var_membrane__V;
  const double var_standard_ionic_concentrations__K_o = 5.4;
  const double var_time_independent_potassium_current__K_o = var_standard_ionic_concentrations__K_o;
  const double var_time_independent_potassium_current__R = var_membrane__R;
  const double var_time_independent_potassium_current__F = var_membrane__F;
  const double var_time_independent_potassium_current__K_i = var_intracellular_ion_concentrations__K_i;
  const double var_time_independent_potassium_current__T = var_membrane__T;
  const double var_time_independent_potassium_current__E_K = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__K_o / var_time_independent_potassium_current__K_i);
  const double var_time_independent_potassium_current__i_K1 = (var_time_independent_potassium_current__Cm * var_time_independent_potassium_current__g_K1 * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp(0.07 * (var_time_independent_potassium_current__V + 80.0)));
  const double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
  const double var_transient_outward_K_current__oi = var_transient_outward_K_current_oi_gate__oi;
  const double var_transient_outward_K_current__Cm = var_membrane__Cm;
  const double var_transient_outward_K_current__V = var_membrane__V;
  const double var_transient_outward_K_current__oa = var_transient_outward_K_current_oa_gate__oa;
  const double var_transient_outward_K_current__g_to = 0.1652;
  const double var_transient_outward_K_current__E_K = var_time_independent_potassium_current__E_K;
  const double var_transient_outward_K_current__i_to = var_transient_outward_K_current__Cm * var_transient_outward_K_current__g_to * pow(var_transient_outward_K_current__oa, 3.0) * var_transient_outward_K_current__oi * (var_transient_outward_K_current__V - var_transient_outward_K_current__E_K);
  const double var_membrane__i_to = var_transient_outward_K_current__i_to;
  const double var_ultrarapid_delayed_rectifier_K_current__V = var_membrane__V;
  const double var_ultrarapid_delayed_rectifier_K_current__g_Kur = 0.005 + (0.05 / (1.0 + exp((var_ultrarapid_delayed_rectifier_K_current__V - 15.0) / (-13.0))));
  const double var_ultrarapid_delayed_rectifier_K_current__Cm = var_membrane__Cm;
  const double var_ultrarapid_delayed_rectifier_K_current__ua = var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua;
  const double var_ultrarapid_delayed_rectifier_K_current__E_K = var_time_independent_potassium_current__E_K;
  const double var_ultrarapid_delayed_rectifier_K_current__ui = var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui;
  const double var_ultrarapid_delayed_rectifier_K_current__i_Kur = var_ultrarapid_delayed_rectifier_K_current__Cm * var_ultrarapid_delayed_rectifier_K_current__g_Kur * pow(var_ultrarapid_delayed_rectifier_K_current__ua, 3.0) * var_ultrarapid_delayed_rectifier_K_current__ui * (var_ultrarapid_delayed_rectifier_K_current__V - var_ultrarapid_delayed_rectifier_K_current__E_K);
  const double var_membrane__i_Kur = var_ultrarapid_delayed_rectifier_K_current__i_Kur;
  const double var_rapid_delayed_rectifier_K_current__g_Kr = 0.029411765;
  const double var_rapid_delayed_rectifier_K_current__E_K = var_time_independent_potassium_current__E_K;
  const double var_rapid_delayed_rectifier_K_current__V = var_membrane__V;
  const double var_rapid_delayed_rectifier_K_current__xr = var_rapid_delayed_rectifier_K_current_xr_gate__xr;
  const double var_rapid_delayed_rectifier_K_current__Cm = var_membrane__Cm;
  const double var_rapid_delayed_rectifier_K_current__i_Kr = (var_rapid_delayed_rectifier_K_current__Cm * var_rapid_delayed_rectifier_K_current__g_Kr * var_rapid_delayed_rectifier_K_current__xr * (var_rapid_delayed_rectifier_K_current__V - var_rapid_delayed_rectifier_K_current__E_K)) / (1.0 + exp((var_rapid_delayed_rectifier_K_current__V + 15.0) / 22.4));
  const double var_membrane__i_Kr = var_rapid_delayed_rectifier_K_current__i_Kr;
  const double var_slow_delayed_rectifier_K_current__E_K = var_time_independent_potassium_current__E_K;
  const double var_slow_delayed_rectifier_K_current__xs = var_slow_delayed_rectifier_K_current_xs_gate__xs;
  const double var_slow_delayed_rectifier_K_current__Cm = var_membrane__Cm;
  const double var_slow_delayed_rectifier_K_current__g_Ks = 0.12941176;
  const double var_slow_delayed_rectifier_K_current__V = var_membrane__V;
  const double var_slow_delayed_rectifier_K_current__i_Ks = var_slow_delayed_rectifier_K_current__Cm * var_slow_delayed_rectifier_K_current__g_Ks * pow(var_slow_delayed_rectifier_K_current__xs, 2.0) * (var_slow_delayed_rectifier_K_current__V - var_slow_delayed_rectifier_K_current__E_K);
  const double var_membrane__i_Ks = var_slow_delayed_rectifier_K_current__i_Ks;
  const double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
  const double var_L_type_Ca_channel__f_Ca = var_L_type_Ca_channel_f_Ca_gate__f_Ca;
  const double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
  const double var_L_type_Ca_channel__Cm = var_membrane__Cm;
  const double var_L_type_Ca_channel__g_Ca_L = 0.12375;
  const double var_L_type_Ca_channel__V = var_membrane__V;
  const double var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__Cm * var_L_type_Ca_channel__g_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * (var_L_type_Ca_channel__V - 65.0);
  const double var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
  const double var_sarcolemmal_calcium_pump_current__i_CaP_max = 0.275;
  const double var_sarcolemmal_calcium_pump_current__Cm = var_membrane__Cm;
  const double var_sarcolemmal_calcium_pump_current__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_sarcolemmal_calcium_pump_current__i_CaP = (var_sarcolemmal_calcium_pump_current__Cm * var_sarcolemmal_calcium_pump_current__i_CaP_max * var_sarcolemmal_calcium_pump_current__Ca_i) / (0.0005 + var_sarcolemmal_calcium_pump_current__Ca_i);
  const double var_membrane__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP;
  const double var_sodium_potassium_pump__Km_Na_i = 10.0;
  const double var_sodium_potassium_pump__i_NaK_max = 0.59933874;
  const double var_sodium_potassium_pump__Km_K_o = 1.5;
  const double var_sodium_potassium_pump__Cm = var_membrane__Cm;
  const double var_sodium_potassium_pump__K_o = var_standard_ionic_concentrations__K_o;
  const double var_sodium_potassium_pump__Na_i = var_intracellular_ion_concentrations__Na_i;
  const double var_sodium_potassium_pump__V = var_membrane__V;
  const double var_sodium_potassium_pump__F = var_membrane__F;
  const double var_sodium_potassium_pump__T = var_membrane__T;
  const double var_sodium_potassium_pump__Na_o = var_standard_ionic_concentrations__Na_o;
  const double var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Na_o / 67.3) - 1.0);
  const double var_sodium_potassium_pump__R = var_membrane__R;
  const double var_sodium_potassium_pump__f_NaK = pow(1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__F * var_sodium_potassium_pump__V) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__F) * var_sodium_potassium_pump__V) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))), -1.0);
  const double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__Cm * var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__f_NaK * 1.0) / (1.0 + pow(var_sodium_potassium_pump__Km_Na_i / var_sodium_potassium_pump__Na_i, 1.5))) * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_o + var_sodium_potassium_pump__Km_K_o);
  const double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
  const double var_Na_Ca_exchanger_current__F = var_membrane__F;
  const double var_Na_Ca_exchanger_current__K_mNa = 87.5;
  const double var_Na_Ca_exchanger_current__I_NaCa_max = 1600.0;
  const double var_Na_Ca_exchanger_current__T = var_membrane__T;
  const double var_Na_Ca_exchanger_current__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_Na_Ca_exchanger_current__K_sat = 0.1;
  const double var_Na_Ca_exchanger_current__gamma = 0.35;
  const double var_standard_ionic_concentrations__Ca_o = 1.8;
  const double var_Na_Ca_exchanger_current__Ca_o = var_standard_ionic_concentrations__Ca_o;
  const double var_Na_Ca_exchanger_current__Na_o = var_standard_ionic_concentrations__Na_o;
  const double var_Na_Ca_exchanger_current__R = var_membrane__R;
  const double var_Na_Ca_exchanger_current__Cm = var_membrane__Cm;
  const double var_Na_Ca_exchanger_current__V = var_membrane__V;
  const double var_Na_Ca_exchanger_current__K_mCa = 1.38;
  const double var_Na_Ca_exchanger_current__Na_i = var_intracellular_ion_concentrations__Na_i;
  const double var_Na_Ca_exchanger_current__i_NaCa = (var_Na_Ca_exchanger_current__Cm * var_Na_Ca_exchanger_current__I_NaCa_max * ((exp((var_Na_Ca_exchanger_current__gamma * var_Na_Ca_exchanger_current__F * var_Na_Ca_exchanger_current__V) / (var_Na_Ca_exchanger_current__R * var_Na_Ca_exchanger_current__T)) * pow(var_Na_Ca_exchanger_current__Na_i, 3.0) * var_Na_Ca_exchanger_current__Ca_o) - (exp(((var_Na_Ca_exchanger_current__gamma - 1.0) * var_Na_Ca_exchanger_current__F * var_Na_Ca_exchanger_current__V) / (var_Na_Ca_exchanger_current__R * var_Na_Ca_exchanger_current__T)) * pow(var_Na_Ca_exchanger_current__Na_o, 3.0) * var_Na_Ca_exchanger_current__Ca_i))) / ((pow(var_Na_Ca_exchanger_current__K_mNa, 3.0) + pow(var_Na_Ca_exchanger_current__Na_o, 3.0)) * (var_Na_Ca_exchanger_current__K_mCa + var_Na_Ca_exchanger_current__Ca_o) * (1.0 + (var_Na_Ca_exchanger_current__K_sat * exp(((var_Na_Ca_exchanger_current__gamma - 1.0) * var_Na_Ca_exchanger_current__V * var_Na_Ca_exchanger_current__F) / (var_Na_Ca_exchanger_current__R * var_Na_Ca_exchanger_current__T)))));
  const double var_membrane__i_NaCa = var_Na_Ca_exchanger_current__i_NaCa;
  const double var_background_currents__Cm = var_membrane__Cm;
  const double var_background_currents__E_Na = var_fast_sodium_current__E_Na;
  const double var_background_currents__g_B_Na = 0.0006744375;
  const double var_background_currents__V = var_membrane__V;
  const double var_background_currents__i_B_Na = var_background_currents__Cm * var_background_currents__g_B_Na * (var_background_currents__V - var_background_currents__E_Na);
  const double var_membrane__i_B_Na = var_background_currents__i_B_Na;
  const double var_background_currents__R = var_membrane__R;
  const double var_background_currents__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_background_currents__Ca_o = var_standard_ionic_concentrations__Ca_o;
  const double var_background_currents__F = var_membrane__F;
  const double var_background_currents__T = var_membrane__T;
  const double var_background_currents__E_Ca = ((var_background_currents__R * var_background_currents__T) / (2.0 * var_background_currents__F)) * log(var_background_currents__Ca_o / var_background_currents__Ca_i);
  const double var_background_currents__g_B_Ca = 0.001131;
  const double var_background_currents__i_B_Ca = var_background_currents__Cm * var_background_currents__g_B_Ca * (var_background_currents__V - var_background_currents__E_Ca);
  const double var_membrane__i_B_Ca = var_background_currents__i_B_Ca;
  const double var_fast_sodium_current__time = var_environment__time;
  const double var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
  const double var_fast_sodium_current_m_gate__alpha_m = (var_fast_sodium_current_m_gate__V == (-47.13)) ? 3.2 : ((0.32 * (var_fast_sodium_current_m_gate__V + 47.13)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13))));
  const double var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0);
  const double var_fast_sodium_current_m_gate__m_inf = var_fast_sodium_current_m_gate__alpha_m / (var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m);
  const double var_fast_sodium_current_m_gate__tau_m = 1.0 / (var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m);
  const double var_fast_sodium_current_m_gate__time = var_fast_sodium_current__time;
  const double var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
  const double var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.135 * exp((var_fast_sodium_current_h_gate__V + 80.0) / (-6.8))) : 0.0;
  const double var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3.56 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1)))));
  const double var_fast_sodium_current_h_gate__h_inf = var_fast_sodium_current_h_gate__alpha_h / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h);
  const double var_fast_sodium_current_h_gate__tau_h = 1.0 / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h);
  const double var_fast_sodium_current_h_gate__time = var_fast_sodium_current__time;
  const double var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
  const double var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? (((((-127140.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0;
  const double var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.1212 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.3 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))));
  const double var_fast_sodium_current_j_gate__j_inf = var_fast_sodium_current_j_gate__alpha_j / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j);
  const double var_fast_sodium_current_j_gate__tau_j = 1.0 / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j);
  const double var_fast_sodium_current_j_gate__time = var_fast_sodium_current__time;
  const double var_transient_outward_K_current__K_Q10 = 3.0;
  const double var_transient_outward_K_current__time = var_environment__time;
  const double var_transient_outward_K_current_oa_gate__V = var_transient_outward_K_current__V;
  const double var_transient_outward_K_current_oa_gate__alpha_oa = 0.65 * pow(exp((var_transient_outward_K_current_oa_gate__V - (-10.0)) / (-8.5)) + exp(((var_transient_outward_K_current_oa_gate__V - (-10.0)) - 40.0) / (-59.0)), -1.0);
  const double var_transient_outward_K_current_oa_gate__beta_oa = 0.65 * pow(2.5 + exp(((var_transient_outward_K_current_oa_gate__V - (-10.0)) + 72.0) / 17.0), -1.0);
  const double var_transient_outward_K_current_oa_gate__K_Q10 = var_transient_outward_K_current__K_Q10;
  const double var_transient_outward_K_current_oa_gate__tau_oa = pow(var_transient_outward_K_current_oa_gate__alpha_oa + var_transient_outward_K_current_oa_gate__beta_oa, -1.0) / var_transient_outward_K_current_oa_gate__K_Q10;
  const double var_transient_outward_K_current_oa_gate__oa_infinity = pow(1.0 + exp(((var_transient_outward_K_current_oa_gate__V - (-10.0)) + 10.47) / (-17.54)), -1.0);
  const double var_transient_outward_K_current_oa_gate__time = var_transient_outward_K_current__time;
  const double var_transient_outward_K_current_oi_gate__V = var_transient_outward_K_current__V;
  const double var_transient_outward_K_current_oi_gate__alpha_oi = pow(18.53 + (1.0 * exp(((var_transient_outward_K_current_oi_gate__V - (-10.0)) + 103.7) / 10.95)), -1.0);
  const double var_transient_outward_K_current_oi_gate__beta_oi = pow(35.56 + (1.0 * exp(((var_transient_outward_K_current_oi_gate__V - (-10.0)) - 8.74) / (-7.44))), -1.0);
  const double var_transient_outward_K_current_oi_gate__K_Q10 = var_transient_outward_K_current__K_Q10;
  const double var_transient_outward_K_current_oi_gate__tau_oi = pow(var_transient_outward_K_current_oi_gate__alpha_oi + var_transient_outward_K_current_oi_gate__beta_oi, -1.0) / var_transient_outward_K_current_oi_gate__K_Q10;
  const double var_transient_outward_K_current_oi_gate__oi_infinity = pow(1.0 + exp(((var_transient_outward_K_current_oi_gate__V - (-10.0)) + 33.1) / 5.3), -1.0);
  const double var_transient_outward_K_current_oi_gate__time = var_transient_outward_K_current__time;
  const double var_ultrarapid_delayed_rectifier_K_current__time = var_environment__time;
  const double var_ultrarapid_delayed_rectifier_K_current__K_Q10 = var_transient_outward_K_current__K_Q10;
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__V = var_ultrarapid_delayed_rectifier_K_current__V;
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__alpha_ua = 0.65 * pow(exp((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) / (-8.5)) + exp(((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) - 40.0) / (-59.0)), -1.0);
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__beta_ua = 0.65 * pow(2.5 + exp(((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) + 72.0) / 17.0), -1.0);
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__K_Q10 = var_ultrarapid_delayed_rectifier_K_current__K_Q10;
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__tau_ua = pow(var_ultrarapid_delayed_rectifier_K_current_ua_gate__alpha_ua + var_ultrarapid_delayed_rectifier_K_current_ua_gate__beta_ua, -1.0) / var_ultrarapid_delayed_rectifier_K_current_ua_gate__K_Q10;
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua_infinity = pow(1.0 + exp(((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) + 20.3) / (-9.6)), -1.0);
  const double var_ultrarapid_delayed_rectifier_K_current_ua_gate__time = var_ultrarapid_delayed_rectifier_K_current__time;
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__V = var_ultrarapid_delayed_rectifier_K_current__V;
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__alpha_ui = pow(21.0 + (1.0 * exp(((var_ultrarapid_delayed_rectifier_K_current_ui_gate__V - (-10.0)) - 195.0) / (-28.0))), -1.0);
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__beta_ui = 1.0 / exp(((var_ultrarapid_delayed_rectifier_K_current_ui_gate__V - (-10.0)) - 168.0) / (-16.0));
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__K_Q10 = var_ultrarapid_delayed_rectifier_K_current__K_Q10;
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__tau_ui = pow(var_ultrarapid_delayed_rectifier_K_current_ui_gate__alpha_ui + var_ultrarapid_delayed_rectifier_K_current_ui_gate__beta_ui, -1.0) / var_ultrarapid_delayed_rectifier_K_current_ui_gate__K_Q10;
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui_infinity = pow(1.0 + exp(((var_ultrarapid_delayed_rectifier_K_current_ui_gate__V - (-10.0)) - 109.45) / 27.48), -1.0);
  const double var_ultrarapid_delayed_rectifier_K_current_ui_gate__time = var_ultrarapid_delayed_rectifier_K_current__time;
  const double var_rapid_delayed_rectifier_K_current__time = var_environment__time;
  const double var_rapid_delayed_rectifier_K_current_xr_gate__V = var_rapid_delayed_rectifier_K_current__V;
  const double var_rapid_delayed_rectifier_K_current_xr_gate__alpha_xr = (fabs(var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1) < 1e-10) ? 0.0015 : ((0.0003 * (var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1)) / (1.0 - exp((var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1) / (-5.0))));
  const double var_rapid_delayed_rectifier_K_current_xr_gate__beta_xr = (fabs(var_rapid_delayed_rectifier_K_current_xr_gate__V - 3.3328) < 1e-10) ? 0.00037836118 : ((7.3898e-05 * (var_rapid_delayed_rectifier_K_current_xr_gate__V - 3.3328)) / (exp((var_rapid_delayed_rectifier_K_current_xr_gate__V - 3.3328) / 5.1237) - 1.0));
  const double var_rapid_delayed_rectifier_K_current_xr_gate__tau_xr = pow(var_rapid_delayed_rectifier_K_current_xr_gate__alpha_xr + var_rapid_delayed_rectifier_K_current_xr_gate__beta_xr, -1.0);
  const double var_rapid_delayed_rectifier_K_current_xr_gate__xr_infinity = pow(1.0 + exp((var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1) / (-6.5)), -1.0);
  const double var_rapid_delayed_rectifier_K_current_xr_gate__time = var_rapid_delayed_rectifier_K_current__time;
  const double var_slow_delayed_rectifier_K_current__time = var_environment__time;
  const double var_slow_delayed_rectifier_K_current_xs_gate__V = var_slow_delayed_rectifier_K_current__V;
  const double var_slow_delayed_rectifier_K_current_xs_gate__alpha_xs = (fabs(var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) < 1e-10) ? 0.00068 : ((4e-05 * (var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9)) / (1.0 - exp((var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) / (-17.0))));
  const double var_slow_delayed_rectifier_K_current_xs_gate__beta_xs = (fabs(var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) < 1e-10) ? 0.000315 : ((3.5e-05 * (var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9)) / (exp((var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) / 9.0) - 1.0));
  const double var_slow_delayed_rectifier_K_current_xs_gate__tau_xs = 0.5 * pow(var_slow_delayed_rectifier_K_current_xs_gate__alpha_xs + var_slow_delayed_rectifier_K_current_xs_gate__beta_xs, -1.0);
  const double var_slow_delayed_rectifier_K_current_xs_gate__xs_infinity = pow(1.0 + exp((var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) / (-12.7)), -0.5);
  const double var_slow_delayed_rectifier_K_current_xs_gate__time = var_slow_delayed_rectifier_K_current__time;
  const double var_L_type_Ca_channel__time = var_environment__time;
  const double var_L_type_Ca_channel__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_L_type_Ca_channel_d_gate__V = var_L_type_Ca_channel__V;
  const double var_L_type_Ca_channel_d_gate__d_infinity = pow(1.0 + exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-8.0)), -1.0);
  const double var_L_type_Ca_channel_d_gate__tau_d = (fabs(var_L_type_Ca_channel_d_gate__V + 10.0) < 1e-10) ? (4.579 / (1.0 + exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-6.24)))) : ((1.0 - exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-6.24))) / (0.035 * (var_L_type_Ca_channel_d_gate__V + 10.0) * (1.0 + exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-6.24)))));
  const double var_L_type_Ca_channel_d_gate__time = var_L_type_Ca_channel__time;
  const double var_L_type_Ca_channel_f_gate__V = var_L_type_Ca_channel__V;
  const double var_L_type_Ca_channel_f_gate__f_infinity = exp((-(var_L_type_Ca_channel_f_gate__V + 28.0)) / 6.9) / (1.0 + exp((-(var_L_type_Ca_channel_f_gate__V + 28.0)) / 6.9));
  const double var_L_type_Ca_channel_f_gate__tau_f = 9.0 * pow((0.0197 * exp((-pow(0.0337, 2.0)) * pow(var_L_type_Ca_channel_f_gate__V + 10.0, 2.0))) + 0.02, -1.0);
  const double var_L_type_Ca_channel_f_gate__time = var_L_type_Ca_channel__time;
  const double var_L_type_Ca_channel_f_Ca_gate__Ca_i = var_L_type_Ca_channel__Ca_i;
  const double var_L_type_Ca_channel_f_Ca_gate__f_Ca_infinity = pow(1.0 + (var_L_type_Ca_channel_f_Ca_gate__Ca_i / 0.00035), -1.0);
  const double var_L_type_Ca_channel_f_Ca_gate__tau_f_Ca = 2.0;
  const double var_L_type_Ca_channel_f_Ca_gate__time = var_L_type_Ca_channel__time;
  const double var_background_currents__g_B_K = 0.0;
  const double var_background_currents__E_K = var_time_independent_potassium_current__E_K;
  const double var_background_currents__i_B_K = var_background_currents__Cm * var_background_currents__g_B_K * (var_background_currents__V - var_background_currents__E_K);
  const double var_Ca_release_current_from_JSR__K_rel = 30.0;
  const double var_Ca_release_current_from_JSR__u = var_Ca_release_current_from_JSR_u_gate__u;
  const double var_Ca_release_current_from_JSR__w = var_Ca_release_current_from_JSR_w_gate__w;
  const double var_Ca_release_current_from_JSR__Ca_rel = var_intracellular_ion_concentrations__Ca_rel;
  const double var_Ca_release_current_from_JSR__v = var_Ca_release_current_from_JSR_v_gate__v;
  const double var_Ca_release_current_from_JSR__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_Ca_release_current_from_JSR__i_rel = var_Ca_release_current_from_JSR__K_rel * pow(var_Ca_release_current_from_JSR__u, 2.0) * var_Ca_release_current_from_JSR__v * var_Ca_release_current_from_JSR__w * (var_Ca_release_current_from_JSR__Ca_rel - var_Ca_release_current_from_JSR__Ca_i);
  const double var_Ca_release_current_from_JSR__F = var_membrane__F;
  const double var_intracellular_ion_concentrations__V_cell = 20100.0;
  const double var_intracellular_ion_concentrations__V_rel = 0.0048 * var_intracellular_ion_concentrations__V_cell;
  const double var_Ca_release_current_from_JSR__V_rel = var_intracellular_ion_concentrations__V_rel;
  const double var_Ca_release_current_from_JSR__i_NaCa = var_Na_Ca_exchanger_current__i_NaCa;
  const double var_Ca_release_current_from_JSR__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
  const double var_Ca_release_current_from_JSR__Fn = 1000.0 * ((1e-15 * var_Ca_release_current_from_JSR__V_rel * var_Ca_release_current_from_JSR__i_rel) - ((1e-15 / (2.0 * var_Ca_release_current_from_JSR__F)) * ((0.5 * var_Ca_release_current_from_JSR__i_Ca_L) - (0.2 * var_Ca_release_current_from_JSR__i_NaCa))));
  const double var_Ca_release_current_from_JSR__time = var_environment__time;
  const double var_Ca_release_current_from_JSR__V = var_membrane__V;
  const double var_Ca_release_current_from_JSR_u_gate__tau_u = 8.0;
  const double var_Ca_release_current_from_JSR_u_gate__Fn = var_Ca_release_current_from_JSR__Fn;
  const double var_Ca_release_current_from_JSR_u_gate__u_infinity = pow(1.0 + exp((-(var_Ca_release_current_from_JSR_u_gate__Fn - 3.4175e-13)) / 1.367e-15), -1.0);
  const double var_Ca_release_current_from_JSR_u_gate__time = var_Ca_release_current_from_JSR__time;
  const double var_Ca_release_current_from_JSR_v_gate__Fn = var_Ca_release_current_from_JSR__Fn;
  const double var_Ca_release_current_from_JSR_v_gate__tau_v = 1.91 + (2.09 * pow(1.0 + exp((-(var_Ca_release_current_from_JSR_v_gate__Fn - 3.4175e-13)) / 1.367e-15), -1.0));
  const double var_Ca_release_current_from_JSR_v_gate__v_infinity = 1.0 - pow(1.0 + exp((-(var_Ca_release_current_from_JSR_v_gate__Fn - 6.835e-14)) / 1.367e-15), -1.0);
  const double var_Ca_release_current_from_JSR_v_gate__time = var_Ca_release_current_from_JSR__time;
  const double var_Ca_release_current_from_JSR_w_gate__V = var_Ca_release_current_from_JSR__V;
  const double var_Ca_release_current_from_JSR_w_gate__tau_w = (fabs(var_Ca_release_current_from_JSR_w_gate__V - 7.9) < 1e-10) ? ((6.0 * 0.2) / 1.3) : ((6.0 * (1.0 - exp((-(var_Ca_release_current_from_JSR_w_gate__V - 7.9)) / 5.0))) / ((1.0 + (0.3 * exp((-(var_Ca_release_current_from_JSR_w_gate__V - 7.9)) / 5.0))) * 1.0 * (var_Ca_release_current_from_JSR_w_gate__V - 7.9)));
  const double var_Ca_release_current_from_JSR_w_gate__w_infinity = 1.0 - pow(1.0 + exp((-(var_Ca_release_current_from_JSR_w_gate__V - 40.0)) / 17.0), -1.0);
  const double var_Ca_release_current_from_JSR_w_gate__time = var_Ca_release_current_from_JSR__time;
  const double var_transfer_current_from_NSR_to_JSR__Ca_up = var_intracellular_ion_concentrations__Ca_up;
  const double var_transfer_current_from_NSR_to_JSR__Ca_rel = var_intracellular_ion_concentrations__Ca_rel;
  const double var_transfer_current_from_NSR_to_JSR__tau_tr = 180.0;
  const double var_transfer_current_from_NSR_to_JSR__i_tr = (var_transfer_current_from_NSR_to_JSR__Ca_up - var_transfer_current_from_NSR_to_JSR__Ca_rel) / var_transfer_current_from_NSR_to_JSR__tau_tr;
  const double var_Ca_uptake_current_by_the_NSR__I_up_max = 0.005;
  const double var_Ca_uptake_current_by_the_NSR__K_up = 0.00092;
  const double var_Ca_uptake_current_by_the_NSR__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_Ca_uptake_current_by_the_NSR__i_up = var_Ca_uptake_current_by_the_NSR__I_up_max / (1.0 + (var_Ca_uptake_current_by_the_NSR__K_up / var_Ca_uptake_current_by_the_NSR__Ca_i));
  const double var_Ca_leak_current_by_the_NSR__Ca_up_max = 15.0;
  const double var_Ca_leak_current_by_the_NSR__I_up_max = var_Ca_uptake_current_by_the_NSR__I_up_max;
  const double var_Ca_leak_current_by_the_NSR__Ca_up = var_intracellular_ion_concentrations__Ca_up;
  const double var_Ca_leak_current_by_the_NSR__i_up_leak = (var_Ca_leak_current_by_the_NSR__I_up_max * var_Ca_leak_current_by_the_NSR__Ca_up) / var_Ca_leak_current_by_the_NSR__Ca_up_max;
  const double var_Ca_buffers__CMDN_max = 0.05;
  const double var_Ca_buffers__TRPN_max = 0.07;
  const double var_Ca_buffers__CSQN_max = 10.0;
  const double var_Ca_buffers__Km_CMDN = 0.00238;
  const double var_Ca_buffers__Km_TRPN = 0.0005;
  const double var_Ca_buffers__Km_CSQN = 0.8;
  const double var_Ca_buffers__Ca_i = var_intracellular_ion_concentrations__Ca_i;
  const double var_Ca_buffers__Ca_rel = var_intracellular_ion_concentrations__Ca_rel;
  const double var_intracellular_ion_concentrations__V_i = var_intracellular_ion_concentrations__V_cell * 0.68;
  const double var_intracellular_ion_concentrations__V_up = 0.0552 * var_intracellular_ion_concentrations__V_cell;
  const double var_intracellular_ion_concentrations__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
  const double var_intracellular_ion_concentrations__i_rel = var_Ca_release_current_from_JSR__i_rel;
  const double var_intracellular_ion_concentrations__i_B_Ca = var_background_currents__i_B_Ca;
  const double var_intracellular_ion_concentrations__F = var_membrane__F;
  const double var_intracellular_ion_concentrations__i_NaCa = var_Na_Ca_exchanger_current__i_NaCa;
  const double var_intracellular_ion_concentrations__i_up_leak = var_Ca_leak_current_by_the_NSR__i_up_leak;
  const double var_intracellular_ion_concentrations__i_up = var_Ca_uptake_current_by_the_NSR__i_up;
  const double var_intracellular_ion_concentrations__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP;
  const double var_intracellular_ion_concentrations__B1 = (((2.0 * var_intracellular_ion_concentrations__i_NaCa) - (var_intracellular_ion_concentrations__i_CaP + var_intracellular_ion_concentrations__i_Ca_L + var_intracellular_ion_concentrations__i_B_Ca)) / (2.0 * var_intracellular_ion_concentrations__V_i * var_intracellular_ion_concentrations__F)) + (((var_intracellular_ion_concentrations__V_up * (var_intracellular_ion_concentrations__i_up_leak - var_intracellular_ion_concentrations__i_up)) + (var_intracellular_ion_concentrations__i_rel * var_intracellular_ion_concentrations__V_rel)) / var_intracellular_ion_concentrations__V_i);
  const double var_intracellular_ion_concentrations__Km_TRPN = var_Ca_buffers__Km_TRPN;
  const double var_intracellular_ion_concentrations__CMDN_max = var_Ca_buffers__CMDN_max;
  const double var_intracellular_ion_concentrations__Km_CMDN = var_Ca_buffers__Km_CMDN;
  const double var_intracellular_ion_concentrations__TRPN_max = var_Ca_buffers__TRPN_max;
  const double var_intracellular_ion_concentrations__B2 = 1.0 + ((var_intracellular_ion_concentrations__TRPN_max * var_intracellular_ion_concentrations__Km_TRPN) / pow(var_intracellular_ion_concentrations__Ca_i + var_intracellular_ion_concentrations__Km_TRPN, 2.0)) + ((var_intracellular_ion_concentrations__CMDN_max * var_intracellular_ion_concentrations__Km_CMDN) / pow(var_intracellular_ion_concentrations__Ca_i + var_intracellular_ion_concentrations__Km_CMDN, 2.0));
  const double var_intracellular_ion_concentrations__time = var_environment__time;
  const double var_intracellular_ion_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK;
  const double var_intracellular_ion_concentrations__i_B_Na = var_background_currents__i_B_Na;
  const double var_intracellular_ion_concentrations__i_Na = var_fast_sodium_current__i_Na;
  const double var_intracellular_ion_concentrations__i_K1 = var_time_independent_potassium_current__i_K1;
  const double var_intracellular_ion_concentrations__i_to = var_transient_outward_K_current__i_to;
  const double var_intracellular_ion_concentrations__i_Kur = var_ultrarapid_delayed_rectifier_K_current__i_Kur;
  const double var_intracellular_ion_concentrations__i_Kr = var_rapid_delayed_rectifier_K_current__i_Kr;
  const double var_intracellular_ion_concentrations__i_Ks = var_slow_delayed_rectifier_K_current__i_Ks;
  const double var_intracellular_ion_concentrations__i_B_K = var_background_currents__i_B_K;
  const double var_intracellular_ion_concentrations__i_tr = var_transfer_current_from_NSR_to_JSR__i_tr;
  const double var_intracellular_ion_concentrations__Km_CSQN = var_Ca_buffers__Km_CSQN;
  const double var_intracellular_ion_concentrations__CSQN_max = var_Ca_buffers__CSQN_max;
  const double d_dt_membrane__V = (-(var_membrane__i_Na + var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kur + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_B_Na + var_membrane__i_B_Ca + var_membrane__i_NaK + var_membrane__i_CaP + var_membrane__i_NaCa + var_membrane__i_Ca_L + var_membrane__i_st)) / var_membrane__Cm;
  const double d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__m_inf - var_fast_sodium_current_m_gate__m) / var_fast_sodium_current_m_gate__tau_m;
  const double d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__h_inf - var_fast_sodium_current_h_gate__h) / var_fast_sodium_current_h_gate__tau_h;
  const double d_dt_fast_sodium_current_j_gate__j = (var_fast_sodium_current_j_gate__j_inf - var_fast_sodium_current_j_gate__j) / var_fast_sodium_current_j_gate__tau_j;
  const double d_dt_transient_outward_K_current_oa_gate__oa = (var_transient_outward_K_current_oa_gate__oa_infinity - var_transient_outward_K_current_oa_gate__oa) / var_transient_outward_K_current_oa_gate__tau_oa;
  const double d_dt_transient_outward_K_current_oi_gate__oi = (var_transient_outward_K_current_oi_gate__oi_infinity - var_transient_outward_K_current_oi_gate__oi) / var_transient_outward_K_current_oi_gate__tau_oi;
  const double d_dt_ultrarapid_delayed_rectifier_K_current_ua_gate__ua = (var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua_infinity - var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua) / var_ultrarapid_delayed_rectifier_K_current_ua_gate__tau_ua;
  const double d_dt_ultrarapid_delayed_rectifier_K_current_ui_gate__ui = (var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui_infinity - var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui) / var_ultrarapid_delayed_rectifier_K_current_ui_gate__tau_ui;
  const double d_dt_rapid_delayed_rectifier_K_current_xr_gate__xr = (var_rapid_delayed_rectifier_K_current_xr_gate__xr_infinity - var_rapid_delayed_rectifier_K_current_xr_gate__xr) / var_rapid_delayed_rectifier_K_current_xr_gate__tau_xr;
  const double d_dt_slow_delayed_rectifier_K_current_xs_gate__xs = (var_slow_delayed_rectifier_K_current_xs_gate__xs_infinity - var_slow_delayed_rectifier_K_current_xs_gate__xs) / var_slow_delayed_rectifier_K_current_xs_gate__tau_xs;
  const double d_dt_L_type_Ca_channel_d_gate__d = (var_L_type_Ca_channel_d_gate__d_infinity - var_L_type_Ca_channel_d_gate__d) / var_L_type_Ca_channel_d_gate__tau_d;
  const double d_dt_L_type_Ca_channel_f_gate__f = (var_L_type_Ca_channel_f_gate__f_infinity - var_L_type_Ca_channel_f_gate__f) / var_L_type_Ca_channel_f_gate__tau_f;
  const double d_dt_L_type_Ca_channel_f_Ca_gate__f_Ca = (var_L_type_Ca_channel_f_Ca_gate__f_Ca_infinity - var_L_type_Ca_channel_f_Ca_gate__f_Ca) / var_L_type_Ca_channel_f_Ca_gate__tau_f_Ca;
  const double d_dt_Ca_release_current_from_JSR_u_gate__u = (var_Ca_release_current_from_JSR_u_gate__u_infinity - var_Ca_release_current_from_JSR_u_gate__u) / var_Ca_release_current_from_JSR_u_gate__tau_u;
  const double d_dt_Ca_release_current_from_JSR_v_gate__v = (var_Ca_release_current_from_JSR_v_gate__v_infinity - var_Ca_release_current_from_JSR_v_gate__v) / var_Ca_release_current_from_JSR_v_gate__tau_v;
  const double d_dt_Ca_release_current_from_JSR_w_gate__w = (var_Ca_release_current_from_JSR_w_gate__w_infinity - var_Ca_release_current_from_JSR_w_gate__w) / var_Ca_release_current_from_JSR_w_gate__tau_w;
  const double d_dt_intracellular_ion_concentrations__Na_i = (((-3.0) * var_intracellular_ion_concentrations__i_NaK) - ((3.0 * var_intracellular_ion_concentrations__i_NaCa) + var_intracellular_ion_concentrations__i_B_Na + var_intracellular_ion_concentrations__i_Na)) / (var_intracellular_ion_concentrations__V_i * var_intracellular_ion_concentrations__F);
  const double d_dt_intracellular_ion_concentrations__K_i = ((2.0 * var_intracellular_ion_concentrations__i_NaK) - (var_intracellular_ion_concentrations__i_K1 + var_intracellular_ion_concentrations__i_to + var_intracellular_ion_concentrations__i_Kur + var_intracellular_ion_concentrations__i_Kr + var_intracellular_ion_concentrations__i_Ks + var_intracellular_ion_concentrations__i_B_K)) / (var_intracellular_ion_concentrations__V_i * var_intracellular_ion_concentrations__F);
  const double d_dt_intracellular_ion_concentrations__Ca_i = var_intracellular_ion_concentrations__B1 / var_intracellular_ion_concentrations__B2;
  const double d_dt_intracellular_ion_concentrations__Ca_up = var_intracellular_ion_concentrations__i_up - (var_intracellular_ion_concentrations__i_up_leak + ((var_intracellular_ion_concentrations__i_tr * var_intracellular_ion_concentrations__V_rel) / var_intracellular_ion_concentrations__V_up));
  const double d_dt_intracellular_ion_concentrations__Ca_rel = (var_intracellular_ion_concentrations__i_tr - var_intracellular_ion_concentrations__i_rel) * pow(1.0 + ((var_intracellular_ion_concentrations__CSQN_max * var_intracellular_ion_concentrations__Km_CSQN) / pow(var_intracellular_ion_concentrations__Ca_rel + var_intracellular_ion_concentrations__Km_CSQN, 2.0)), -1.0);
        

  rDY[0] = d_dt_membrane__V;
  rDY[1] = d_dt_fast_sodium_current_m_gate__m;
  rDY[2] = d_dt_fast_sodium_current_h_gate__h;
  rDY[3] = d_dt_fast_sodium_current_j_gate__j;
  rDY[4] = d_dt_transient_outward_K_current_oa_gate__oa;
  rDY[5] = d_dt_transient_outward_K_current_oi_gate__oi;
  rDY[6] = d_dt_ultrarapid_delayed_rectifier_K_current_ua_gate__ua;
  rDY[7] = d_dt_ultrarapid_delayed_rectifier_K_current_ui_gate__ui;
  rDY[8] = d_dt_rapid_delayed_rectifier_K_current_xr_gate__xr;
  rDY[9] = d_dt_slow_delayed_rectifier_K_current_xs_gate__xs;
  rDY[10] = d_dt_L_type_Ca_channel_d_gate__d;
  rDY[11] = d_dt_L_type_Ca_channel_f_gate__f;
  rDY[12] = d_dt_L_type_Ca_channel_f_Ca_gate__f_Ca;
  rDY[13] = d_dt_Ca_release_current_from_JSR_u_gate__u;
  rDY[14] = d_dt_Ca_release_current_from_JSR_v_gate__v;
  rDY[15] = d_dt_Ca_release_current_from_JSR_w_gate__w;
  rDY[16] = d_dt_intracellular_ion_concentrations__Na_i;
  rDY[17] = d_dt_intracellular_ion_concentrations__Ca_i;
  rDY[18] = d_dt_intracellular_ion_concentrations__K_i;
  rDY[19] = d_dt_intracellular_ion_concentrations__Ca_rel;
  rDY[20] = d_dt_intracellular_ion_concentrations__Ca_up;
  */

}

