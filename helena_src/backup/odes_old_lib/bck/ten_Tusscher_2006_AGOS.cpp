#include "ten_Tusscher_2006_AGOS.h"


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
   // R_prime; Units: dimensionless; Initial value: 0.8978
   IC[13] = 0.8978;
   // Ca_i; Units: millimolar; Initial value: 0.000153
   IC[14] = 0.000153;
   // Ca_SR; Units: millimolar; Initial value: 4.272
   IC[15] = 4.272;
   // Ca_ss; Units: millimolar; Initial value: 0.00042
   IC[16] = 0.00042;
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

	const float R = 8.3144720e+03;
	const float T = 3.10e+02;
	const float F = 9.648534149999999e+04;
	const float Na_o = 1.40e+02;
	const float K_o = 5.40e+00;
	const float P_kna = 3.0e-02;
	const float Ca_o = 2.0e+00;
	const float g_K1 = 5.4050e+00;
	const float g_Kr = 1.530e-01;
	const float g_Ks = 9.80e-02;
	const float g_Na = 1.48380e+01;
	const float g_bna = 2.90e-04;
	const float g_CaL = 3.980e-05;
	const float g_bca = 5.920e-04;
	const float g_to = 2.940e-01;
	const float P_NaK = 2.7240e+00;
	const float K_mk = 1.0e+00;
	const float K_mNa = 4.0e+01;
	const float K_NaCa = 1.0e+03;
	const float gamma = 3.50e-01;
	const float alpha = 2.50e+00;
	const float Km_Nai = 8.750e+01;
	const float Km_Ca = 1.380e+00;
	const float K_sat = 1.0e-01;
	const float g_pCa = 1.2380e-01;
	const float K_pCa = 5.0e-04;
	const float g_pK = 1.460e-02;
	const float V_rel = 1.020e-01;
	const float Vmax_up = 6.3750e-03;
	const float K_up = 2.50e-04;
	const float V_leak = 3.60e-04;
	const float V_xfer = 3.80e-03;
	const float k3 = 6.0e-02;
	const float k4 = 5.0e-03;
	const float k1_prime = 1.50e-01;
	const float k2_prime = 4.50e-02;
	const float max_sr = 2.50e+00;
	const float min_sr = 1.0e+00;
	const float EC = 1.50e+00;
	const float Buf_c = 2.0e-01;
	const float K_buf_c = 1.0e-03;
	const float Buf_sr = 1.0e+01;
	const float K_buf_sr = 3.0e-01;
	const float Buf_ss = 4.0e-01;
	const float K_buf_ss = 2.50e-04;
	const float V_sr = 1.0940e-03;
	const float V_c = 1.64040e-02;
	const float Cm = 1.850e-01;
	const float V_ss = 5.4680e-05;


	float V_old_       = rY[0];
	float Xr1_old_     = rY[1];
	float Xr2_old_     = rY[2];
	float Xs_old_      = rY[3];
	float m_old_       = rY[4];
	float h_old_       = rY[5];
	float j_old_       = rY[6];
	float d_old_       = rY[7];
	float f_old_       = rY[8];
	float f2_old_      = rY[9];
	float fCass_old_   = rY[10];
	float s_old_       = rY[11];
	float r_old_       = rY[12];
	float R_prime_old_ = rY[13];
	float Ca_i_old_    = rY[14];
	float Ca_SR_old_   = rY[15];
	float Ca_ss_old_   = rY[16];
	float Na_i_old_    = rY[17];
	float K_i_old_     = rY[18];

	// ---------------------------------------------------------------------------------- //
	float calc_i_Stim = 0.0; //(((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration)))? stim_amplitude: 0.0;	//0
	float calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));	//2
	float calc_E_K = (((R*T)/F)*log((K_o/K_i_old_)));	//3
	float calc_E_Ks = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_)))));	//4
	float calc_E_Ca = (((5.000000000000000e-01*R*T)/F)*log((Ca_o/Ca_i_old_)));	//5
	float calc_xr1_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-2.600000000000000e+01)-V_old_)/7.000000000000000e+00))));	//11
	float calc_alpha_xr1 = (4.500000000000000e+02/(1.000000000000000e+00+exp((((-4.500000000000000e+01)-V_old_)/1.000000000000000e+01))));	//12
	float calc_beta_xr1 = (6.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.150000000000000e+01))));	//13
	float calc_xr2_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+8.800000000000000e+01)/2.400000000000000e+01))));	//16
	float calc_alpha_xr2 = (3.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/2.000000000000000e+01))));	//17
	float calc_beta_xr2 = (1.120000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01))));	//18
	float calc_xs_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/1.400000000000000e+01))));	//22
	float calc_alpha_xs = (1.400000000000000e+03/pow((1.000000000000000e+00+exp(((5.000000000000000e+00-V_old_)/6.000000000000000e+00))),1.0/2.0));	//23
	float calc_beta_xs = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_-3.500000000000000e+01)/1.500000000000000e+01))));	//24
	float calc_m_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp((((-5.686000000000000e+01)-V_old_)/9.029999999999999e+00))),2.000000000000000e+00));	//28
	float calc_alpha_m = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/5.000000000000000e+00))));	//29
	float calc_beta_m = ((1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_+3.500000000000000e+01)/5.000000000000000e+00))))+(1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_-5.000000000000000e+01)/2.000000000000000e+02)))));	//30
	float calc_h_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));	//33
	float calc_alpha_h = ((V_old_<(-4.00e+01)))?((5.70e-02*exp(((-(V_old_+8.00e+01))/6.80e+00)))): 0.0;	//34// 1********************
	float calc_beta_h = ((V_old_<(-4.00e+01)))? (((2.70e+00*exp((7.90e-02*V_old_)))+(3.10e+05*exp((3.4850e-01*V_old_))))):  ((7.70e-01/(1.30e-01*(1.00e+00+exp(((V_old_+1.0660e+01)/(-1.110e+01)))))));	//35// 1********************
	float calc_j_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));	//38
	float calc_alpha_j = ((V_old_<(-4.00e+01)))?(((((((-2.5428000e+04)*exp((2.4440e-01*V_old_)))-(6.9480e-06*exp(((-4.3910e-02)*V_old_))))*(V_old_+3.7780e+01))/1.00e+00)/(1.00e+00+exp((3.110e-01*(V_old_+7.9230e+01)))))): 0.0;	//39// 1********************
	float calc_beta_j = ((V_old_<(-4.00e+01)))? (((2.4240e-02*exp(((-1.0520e-02)*V_old_)))/(1.00e+00+exp(((-1.3780e-01)*(V_old_+4.0140e+01)))))): (((6.00e-01*exp((5.70e-02*V_old_)))/(1.00e+00+exp(((-1.00e-01)*(V_old_+3.20e+01)))))); //40// 1********************
	float calc_i_CaL = ((((g_CaL*d_old_*f_old_*f2_old_*fCass_old_*4.000000000000000e+00*(V_old_-1.500000000000000e+01)*pow(F,2.000000000000000e+00))/(R*T))*((2.500000000000000e-01*Ca_ss_old_*exp(((2.000000000000000e+00*(V_old_-1.500000000000000e+01)*F)/(R*T))))-Ca_o))/(exp(((2.000000000000000e+00*(V_old_-1.500000000000000e+01)*F)/(R*T)))-1.000000000000000e+00));	//44
	float calc_d_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-8.000000000000000e+00)-V_old_)/7.500000000000000e+00))));	//45
	float calc_alpha_d = ((1.400000000000000e+00/(1.000000000000000e+00+exp((((-3.500000000000000e+01)-V_old_)/1.300000000000000e+01))))+2.500000000000000e-01);	//46
	float calc_beta_d = (1.400000000000000e+00/(1.000000000000000e+00+exp(((V_old_+5.000000000000000e+00)/5.000000000000000e+00))));	//47
	float calc_gamma_d = (1.000000000000000e+00/(1.000000000000000e+00+exp(((5.000000000000000e+01-V_old_)/2.000000000000000e+01))));	//48
	float calc_f_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/7.000000000000000e+00))));	//51
	float calc_tau_f = ((1.102500000000000e+03*exp(((-pow((V_old_+2.700000000000000e+01),2.000000000000000e+00))/2.250000000000000e+02)))+(2.000000000000000e+02/(1.000000000000000e+00+exp(((1.300000000000000e+01-V_old_)/1.000000000000000e+01))))+(1.800000000000000e+02/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.000000000000000e+01))))+2.000000000000000e+01);	//52
	float calc_f2_inf = ((6.700000000000000e-01/(1.000000000000000e+00+exp(((V_old_+3.500000000000000e+01)/7.000000000000000e+00))))+3.300000000000000e-01);	//54
	float calc_tau_f2 = ((5.620000000000000e+02*exp(((-pow((V_old_+2.700000000000000e+01),2.000000000000000e+00))/2.400000000000000e+02)))+(3.100000000000000e+01/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/1.000000000000000e+01))))+(8.000000000000000e+01/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.000000000000000e+01)))));	//55
	float calc_fCass_inf = ((6.000000000000000e-01/(1.000000000000000e+00+pow((Ca_ss_old_/5.000000000000000e-02),2.000000000000000e+00)))+4.000000000000000e-01);	//57
	float calc_tau_fCass = ((8.000000000000000e+01/(1.000000000000000e+00+pow((Ca_ss_old_/5.000000000000000e-02),2.000000000000000e+00)))+2.000000000000000e+00);	//58
	float calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/5.000000000000000e+00))));	//62
	float calc_tau_s = ((8.500000000000000e+01*exp(((-pow((V_old_+4.500000000000000e+01),2.000000000000000e+00))/3.200000000000000e+02)))+(5.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_-2.000000000000000e+01)/5.000000000000000e+00))))+3.000000000000000e+00);	//63
	float calc_r_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((2.000000000000000e+01-V_old_)/6.000000000000000e+00))));	//65
	float calc_tau_r = ((9.500000000000000e+00*exp(((-pow((V_old_+4.000000000000000e+01),2.000000000000000e+00))/1.800000000000000e+03)))+8.000000000000000e-01);	//66
	float calc_i_NaK = (((((P_NaK*K_o)/(K_o+K_mk))*Na_i_old_)/(Na_i_old_+K_mNa))/(1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*V_old_*F)/(R*T))))+(3.530000000000000e-02*exp((((-V_old_)*F)/(R*T))))));	//68
	float calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F)/(R*T)))*pow(Na_i_old_,3.000000000000000e+00)*Ca_o)-(exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,3.000000000000000e+00)*Ca_i_old_*alpha)))/((pow(Km_Nai,3.000000000000000e+00)+pow(Na_o,3.000000000000000e+00))*(Km_Ca+Ca_o)*(1.000000000000000e+00+(K_sat*exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))))));	//69
	float calc_i_p_Ca = ((g_pCa*Ca_i_old_)/(Ca_i_old_+K_pCa));	//70
	float calc_i_up = (Vmax_up/(1.000000000000000e+00+(pow(K_up,2.000000000000000e+00)/pow(Ca_i_old_,2.000000000000000e+00))));	//73
	float calc_i_leak = (V_leak*(Ca_SR_old_-Ca_i_old_));	//74
	float calc_i_xfer = (V_xfer*(Ca_ss_old_-Ca_i_old_));	//75
	float calc_kcasr = (max_sr-((max_sr-min_sr)/(1.000000000000000e+00+pow((EC/Ca_SR_old_),2.000000000000000e+00))));	//80
	float calc_Ca_i_bufc = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_c*K_buf_c)/pow((Ca_i_old_+K_buf_c),2.000000000000000e+00))));	//81
	float calc_Ca_sr_bufsr = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_sr*K_buf_sr)/pow((Ca_SR_old_+K_buf_sr),2.000000000000000e+00))));	//82
	float calc_Ca_ss_bufss = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_ss*K_buf_ss)/pow((Ca_ss_old_+K_buf_ss),2.000000000000000e+00))));	//83
	float calc_alpha_K1 = (1.000000000000000e-01/(1.000000000000000e+00+exp((6.000000000000000e-02*((V_old_-calc_E_K)-2.000000000000000e+02)))));	//6
	float calc_beta_K1 = (((3.000000000000000e+00*exp((2.000000000000000e-04*((V_old_-calc_E_K)+1.000000000000000e+02))))+exp((1.000000000000000e-01*((V_old_-calc_E_K)-1.000000000000000e+01))))/(1.000000000000000e+00+exp(((-5.000000000000000e-01)*(V_old_-calc_E_K)))));	//7
	float calc_i_Kr = (g_Kr*Xr1_old_*Xr2_old_*(V_old_-calc_E_K)*pow((K_o/5.400000000000000e+00),1.0/2.0));	//10
	float calc_tau_xr1 = (1.000000000000000e+00*calc_alpha_xr1*calc_beta_xr1);	//14
	float calc_tau_xr2 = (1.000000000000000e+00*calc_alpha_xr2*calc_beta_xr2);	//19
	float calc_i_Ks = (g_Ks*pow(Xs_old_,2.000000000000000e+00)*(V_old_-calc_E_Ks));	//21
	float calc_tau_xs = ((1.000000000000000e+00*calc_alpha_xs*calc_beta_xs)+8.000000000000000e+01);	//25
	float calc_i_Na = (g_Na*pow(m_old_,3.000000000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));	//27
	float calc_tau_m = (1.000000000000000e+00*calc_alpha_m*calc_beta_m);	//31
	float calc_tau_h = (1.000000000000000e+00/(calc_alpha_h+calc_beta_h));	//36
	float calc_tau_j = (1.000000000000000e+00/(calc_alpha_j+calc_beta_j));	//41
	float calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));	//43
	float calc_tau_d = ((1.000000000000000e+00*calc_alpha_d*calc_beta_d)+calc_gamma_d);	//49
	float calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));	//60
	float calc_i_to = (g_to*r_old_*s_old_*(V_old_-calc_E_K));	//61
	float calc_i_p_K = ((g_pK*(V_old_-calc_E_K))/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/5.980000000000000e+00))));	//71
	float calc_xK1_inf = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));	//8
	float calc_k1 = (k1_prime/calc_kcasr);	//78
	float calc_k2 = (k2_prime*calc_kcasr);	//79
	float calc_i_K1 = (g_K1*calc_xK1_inf*(V_old_-calc_E_K));	//9
	float calc_O = ((calc_k1*pow(Ca_ss_old_,2.000000000000000e+00)*R_prime_old_)/(k3+(calc_k1*pow(Ca_ss_old_,2.000000000000000e+00))));	//76
	float calc_i_rel = (V_rel*calc_O*(Ca_SR_old_-Ca_ss_old_));	//72



	// ---------------------------------------------------------------------------------- //
	float dV_new_=((-(calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_CaL+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_NaCa+calc_i_b_Ca+calc_i_p_K+calc_i_p_Ca+calc_i_Stim)));
	float dXr1_new_=(((calc_xr1_inf-Xr1_old_)/calc_tau_xr1));	// 15
	float dXr2_new_=(((calc_xr2_inf-Xr2_old_)/calc_tau_xr2));	// 20
	float dXs_new_=(((calc_xs_inf-Xs_old_)/calc_tau_xs));	// 26
	float dm_new_=(((calc_m_inf-m_old_)/calc_tau_m));	// 32
	float dh_new_=(((calc_h_inf-h_old_)/calc_tau_h));	// 37
	float dj_new_=(((calc_j_inf-j_old_)/calc_tau_j));	// 42
	float dd_new_=(((calc_d_inf-d_old_)/calc_tau_d));	// 50
	float df_new_=(((calc_f_inf-f_old_)/calc_tau_f));	// 53
	float df2_new_=(((calc_f2_inf-f2_old_)/calc_tau_f2));	// 56
	float dfCass_new_=(((calc_fCass_inf-fCass_old_)/calc_tau_fCass));	// 59
	float ds_new_=(((calc_s_inf-s_old_)/calc_tau_s));	// 64
	float dr_new_=(((calc_r_inf-r_old_)/calc_tau_r));	// 67
	float dR_prime_new_=((((-calc_k2)*Ca_ss_old_*R_prime_old_)+(k4*(1.0e+00-R_prime_old_))));	// 77
	float dCa_i_new_=((calc_Ca_i_bufc*(((((calc_i_leak-calc_i_up)*V_sr)/V_c)+calc_i_xfer)-((((calc_i_b_Ca+calc_i_p_Ca)-(2.0e+00*calc_i_NaCa))*Cm)/(2.0e+00*V_c*F)))));	// 84
	float dCa_SR_new_=((calc_Ca_sr_bufsr*(calc_i_up-(calc_i_rel+calc_i_leak))));	// 85
	float dCa_ss_new_=((calc_Ca_ss_bufss*(((((-calc_i_CaL)*Cm)/(2.0e+00*V_ss*F))+((calc_i_rel*V_sr)/V_ss))-((calc_i_xfer*V_c)/V_ss))));	// 86
	float dNa_i_new_=((((-(calc_i_Na+calc_i_b_Na+(3.0e+00*calc_i_NaK)+(3.0e+00*calc_i_NaCa)))/(V_c*F))*Cm));	// 87
	float dK_i_new_=((((-((calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_p_K+calc_i_Stim)-(2.0e+00*calc_i_NaK)))/(V_c*F))*Cm));	// 88
	// ---------------------------------------------------------------------------------- //


	// ------------------------------------------- //
	rDY[0] =  dV_new_;
	rDY[1] =  dXr1_new_;
	rDY[2] =  dXr2_new_;
	rDY[3] =  dXs_new_;
	rDY[4] =  dm_new_;
	rDY[5] =  dh_new_;
	rDY[6] =  dj_new_;
	rDY[7] =  dd_new_;
	rDY[8] =  df_new_;
	rDY[9] =  df2_new_;
	rDY[10] =  dfCass_new_;
	rDY[11] =  ds_new_;
	rDY[12] =  dr_new_;
	rDY[13] =  dR_prime_new_;
	rDY[14] =  dCa_i_new_;
	rDY[15] =  dCa_SR_new_;
	rDY[16] =  dCa_ss_new_;
	rDY[17] =  dNa_i_new_;
	rDY[18] =  dK_i_new_;
	// ------------------------------------------- //

	//	float V_old_       = rY[0];
	//	float Xr1_old_     = rY[1];
	//	float Xr2_old_     = rY[2];
	//	float Xs_old_      = rY[3];
	//	float m_old_       = rY[4];
	//	float h_old_       = rY[5];
	//	float j_old_       = rY[6];
	//	float d_old_       = rY[7];
	//	float f_old_       = rY[8];
	//	float f2_old_      = rY[9];
	//	float fCass_old_   = rY[10];
	//	float s_old_       = rY[11];
	//	float r_old_       = rY[12];
	//	float R_prime_old_ = rY[13];
	//	float Ca_i_old_    = rY[14];
	//	float Ca_SR_old_   = rY[15];
	//	float Ca_ss_old_   = rY[16];
	//	float Na_i_old_    = rY[17];
	//	float K_i_old_     = rY[18];


//	Real var_membrane__V = rY[0];
//	   Real var_rapid_time_dependent_potassium_current_Xr1_gate__Xr1 = rY[1];
//	   Real var_rapid_time_dependent_potassium_current_Xr2_gate__Xr2 = rY[2];
//	   Real var_slow_time_dependent_potassium_current_Xs_gate__Xs = rY[3];
//	   Real var_fast_sodium_current_m_gate__m = rY[4];
//	   Real var_fast_sodium_current_h_gate__h = rY[5];
//	   Real var_fast_sodium_current_j_gate__j = rY[6];
//	   Real var_L_type_Ca_current_d_gate__d = rY[7];
//	   Real var_L_type_Ca_current_f_gate__f = rY[8];
//	   Real var_L_type_Ca_current_f2_gate__f2 = rY[9];
//	   Real var_L_type_Ca_current_fCass_gate__fCass = rY[10];
//	   Real var_transient_outward_current_s_gate__s = rY[11];
//	   Real var_transient_outward_current_r_gate__r = rY[12];
//	   Real var_calcium_dynamics__Ca_i = rY[13];
//	   Real var_calcium_dynamics__Ca_SR = rY[14];
//	   Real var_calcium_dynamics__Ca_ss = rY[15];
//	   Real var_calcium_dynamics__R_prime = rY[16];
//	   Real var_sodium_dynamics__Na_i = rY[17];
//	   Real var_potassium_dynamics__K_i = rY[18];
//
//	   // Inputs:
//	   // Time units: ms
//	   Real var_environment__time = time;
//
//	   // Mathematics
//	   const Real var_membrane__R = 8314.472;
//	   const Real var_membrane__T = 310.0;
//	   const Real var_membrane__F = 96485.3415;
//	   const Real var_membrane__Cm = 0.185;
//	   const Real var_membrane__V_c = 0.016404;
//	   Real var_membrane__time = var_environment__time;
//	   Real var_inward_rectifier_potassium_current__V = var_membrane__V;
//	   Real var_reversal_potentials__K_i = var_potassium_dynamics__K_i;
//	   Real var_reversal_potentials__R = var_membrane__R;
//	   Real var_reversal_potentials__T = var_membrane__T;
//	   Real var_reversal_potentials__F = var_membrane__F;
//	   const Real var_potassium_dynamics__K_o = 5.4;
//	   Real var_reversal_potentials__K_o = var_potassium_dynamics__K_o;
//	   Real var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
//	   Real var_inward_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
//	   Real var_inward_rectifier_potassium_current__beta_K1 = ((3.0 * exp(0.0002 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) + 100.0))) + exp(0.1 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) - 10.0))) / (1.0 + exp((-0.5) * (var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K)));
//	   Real var_inward_rectifier_potassium_current__alpha_K1 = 0.1 / (1.0 + exp(0.06 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) - 200.0)));
//	   Real var_inward_rectifier_potassium_current__xK1_inf = var_inward_rectifier_potassium_current__alpha_K1 / (var_inward_rectifier_potassium_current__alpha_K1 + var_inward_rectifier_potassium_current__beta_K1);
//	   const Real var_inward_rectifier_potassium_current__g_K1 = 5.405;
//	   Real var_inward_rectifier_potassium_current__i_K1 = var_inward_rectifier_potassium_current__g_K1 * var_inward_rectifier_potassium_current__xK1_inf * (var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K);
//	   Real var_membrane__i_K1 = var_inward_rectifier_potassium_current__i_K1;
//	   Real var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
//	   Real var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
//	   const Real var_transient_outward_current__g_to = 0.294;
//	   Real var_transient_outward_current__V = var_membrane__V;
//	   Real var_transient_outward_current__E_K = var_reversal_potentials__E_K;
//	   Real var_transient_outward_current__i_to = var_transient_outward_current__g_to * var_transient_outward_current__r * var_transient_outward_current__s * (var_transient_outward_current__V - var_transient_outward_current__E_K);
//	   Real var_membrane__i_to = var_transient_outward_current__i_to;
//	   Real var_rapid_time_dependent_potassium_current__Xr1 = var_rapid_time_dependent_potassium_current_Xr1_gate__Xr1;
//	   Real var_rapid_time_dependent_potassium_current__V = var_membrane__V;
//	   Real var_rapid_time_dependent_potassium_current__Xr2 = var_rapid_time_dependent_potassium_current_Xr2_gate__Xr2;
//	   Real var_rapid_time_dependent_potassium_current__K_o = var_potassium_dynamics__K_o;
//	   Real var_rapid_time_dependent_potassium_current__E_K = var_reversal_potentials__E_K;
//	   const Real var_rapid_time_dependent_potassium_current__g_Kr = 0.153;
//	   Real var_rapid_time_dependent_potassium_current__i_Kr = var_rapid_time_dependent_potassium_current__g_Kr * var_rapid_time_dependent_potassium_current__Xr1 * var_rapid_time_dependent_potassium_current__Xr2 * (var_rapid_time_dependent_potassium_current__V - var_rapid_time_dependent_potassium_current__E_K) * sqrt(var_rapid_time_dependent_potassium_current__K_o / 5.4);
//	   Real var_membrane__i_Kr = var_rapid_time_dependent_potassium_current__i_Kr;
//	   const Real var_slow_time_dependent_potassium_current__g_Ks = 0.098;
//	   const Real var_sodium_dynamics__Na_o = 140.0;
//	   Real var_reversal_potentials__Na_o = var_sodium_dynamics__Na_o;
//	   Real var_reversal_potentials__Na_i = var_sodium_dynamics__Na_i;
//	   const Real var_reversal_potentials__P_kna = 0.03;
//	   Real var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
//	   Real var_slow_time_dependent_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
//	   Real var_slow_time_dependent_potassium_current__V = var_membrane__V;
//	   Real var_slow_time_dependent_potassium_current__Xs = var_slow_time_dependent_potassium_current_Xs_gate__Xs;
//	   Real var_slow_time_dependent_potassium_current__i_Ks = var_slow_time_dependent_potassium_current__g_Ks * pow(var_slow_time_dependent_potassium_current__Xs, 2) * (var_slow_time_dependent_potassium_current__V - var_slow_time_dependent_potassium_current__E_Ks);
//	   Real var_membrane__i_Ks = var_slow_time_dependent_potassium_current__i_Ks;
//	   Real var_L_type_Ca_current__Ca_ss = var_calcium_dynamics__Ca_ss;
//	   const Real var_L_type_Ca_current__g_CaL = 3.98e-05;
//	   Real var_L_type_Ca_current__f = var_L_type_Ca_current_f_gate__f;
//	   Real var_L_type_Ca_current__d = var_L_type_Ca_current_d_gate__d;
//	   Real var_L_type_Ca_current__F = var_membrane__F;
//	   Real var_L_type_Ca_current__f2 = var_L_type_Ca_current_f2_gate__f2;
//	   Real var_L_type_Ca_current__fCass = var_L_type_Ca_current_fCass_gate__fCass;
//	   Real var_L_type_Ca_current__V = var_membrane__V;
//	   Real var_L_type_Ca_current__T = var_membrane__T;
//	   const Real var_calcium_dynamics__Ca_o = 2.0;
//	   Real var_L_type_Ca_current__Ca_o = var_calcium_dynamics__Ca_o;
//	   Real var_L_type_Ca_current__R = var_membrane__R;
//	   Real var_L_type_Ca_current__i_CaL = (var_L_type_Ca_current__V != 15.0) ? ((((var_L_type_Ca_current__g_CaL * var_L_type_Ca_current__d * var_L_type_Ca_current__f * var_L_type_Ca_current__f2 * var_L_type_Ca_current__fCass * 4.0 * (var_L_type_Ca_current__V - 15.0) * pow(var_L_type_Ca_current__F, 2)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((0.25 * var_L_type_Ca_current__Ca_ss * exp((2.0 * (var_L_type_Ca_current__V - 15.0) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__Ca_o)) / (exp((2.0 * (var_L_type_Ca_current__V - 15.0) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0)) : ((((var_L_type_Ca_current__g_CaL * var_L_type_Ca_current__d * var_L_type_Ca_current__f * var_L_type_Ca_current__f2 * var_L_type_Ca_current__fCass * 4.0 * (var_L_type_Ca_current__V - 15.0 + 1.0e-06) * pow(var_L_type_Ca_current__F, 2)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((0.25 * var_L_type_Ca_current__Ca_ss * exp((2.0 * (var_L_type_Ca_current__V - 15.0 + 1.0e-06) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__Ca_o)) / (exp((2.0 * (var_L_type_Ca_current__V - 15.0 + 1.0e-06) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0));
//	   Real var_membrane__i_CaL = var_L_type_Ca_current__i_CaL;
//	   Real var_sodium_potassium_pump_current__Na_i = var_sodium_dynamics__Na_i;
//	   Real var_sodium_potassium_pump_current__R = var_membrane__R;
//	   Real var_sodium_potassium_pump_current__T = var_membrane__T;
//	   Real var_sodium_potassium_pump_current__V = var_membrane__V;
//	   const Real var_sodium_potassium_pump_current__K_mk = 1.0;
//	   const Real var_sodium_potassium_pump_current__P_NaK = 2.724;
//	   const Real var_sodium_potassium_pump_current__K_mNa = 40.0;
//	   Real var_sodium_potassium_pump_current__F = var_membrane__F;
//	   Real var_sodium_potassium_pump_current__K_o = var_potassium_dynamics__K_o;
//	   Real var_sodium_potassium_pump_current__i_NaK = ((((var_sodium_potassium_pump_current__P_NaK * var_sodium_potassium_pump_current__K_o) / (var_sodium_potassium_pump_current__K_o + var_sodium_potassium_pump_current__K_mk)) * var_sodium_potassium_pump_current__Na_i) / (var_sodium_potassium_pump_current__Na_i + var_sodium_potassium_pump_current__K_mNa)) / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump_current__V * var_sodium_potassium_pump_current__F) / (var_sodium_potassium_pump_current__R * var_sodium_potassium_pump_current__T))) + (0.0353 * exp(((-var_sodium_potassium_pump_current__V) * var_sodium_potassium_pump_current__F) / (var_sodium_potassium_pump_current__R * var_sodium_potassium_pump_current__T))));
//	   Real var_membrane__i_NaK = var_sodium_potassium_pump_current__i_NaK;
//	   Real var_fast_sodium_current__j = var_fast_sodium_current_j_gate__j;
//	   Real var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
//	   const Real var_fast_sodium_current__g_Na = 14.838;
//	   Real var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
//	   Real var_fast_sodium_current__V = var_membrane__V;
//	   Real var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
//	   Real var_fast_sodium_current__E_Na = var_reversal_potentials__E_Na;
//	   Real var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na);
//	   Real var_membrane__i_Na = var_fast_sodium_current__i_Na;
//	   Real var_sodium_background_current__V = var_membrane__V;
//	   Real var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
//	   const Real var_sodium_background_current__g_bna = 0.00029;
//	   Real var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
//	   Real var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
//	   const Real var_sodium_calcium_exchanger_current__alpha = 2.5;
//	   const Real var_sodium_calcium_exchanger_current__gamma = 0.35;
//	   const Real var_sodium_calcium_exchanger_current__K_sat = 0.1;
//	   const Real var_sodium_calcium_exchanger_current__Km_Ca = 1.38;
//	   const Real var_sodium_calcium_exchanger_current__K_NaCa = 1000.0;
//	   Real var_sodium_calcium_exchanger_current__F = var_membrane__F;
//	   Real var_sodium_calcium_exchanger_current__Ca_i = var_calcium_dynamics__Ca_i;
//	   Real var_sodium_calcium_exchanger_current__Ca_o = var_calcium_dynamics__Ca_o;
//	   Real var_sodium_calcium_exchanger_current__V = var_membrane__V;
//	   Real var_sodium_calcium_exchanger_current__R = var_membrane__R;
//	   const Real var_sodium_calcium_exchanger_current__Km_Nai = 87.5;
//	   Real var_sodium_calcium_exchanger_current__Na_o = var_sodium_dynamics__Na_o;
//	   Real var_sodium_calcium_exchanger_current__Na_i = var_sodium_dynamics__Na_i;
//	   Real var_sodium_calcium_exchanger_current__T = var_membrane__T;
//	   Real var_sodium_calcium_exchanger_current__i_NaCa = (var_sodium_calcium_exchanger_current__K_NaCa * ((exp((var_sodium_calcium_exchanger_current__gamma * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)) * pow(var_sodium_calcium_exchanger_current__Na_i, 3) * var_sodium_calcium_exchanger_current__Ca_o) - (exp(((var_sodium_calcium_exchanger_current__gamma - 1.0) * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)) * pow(var_sodium_calcium_exchanger_current__Na_o, 3) * var_sodium_calcium_exchanger_current__Ca_i * var_sodium_calcium_exchanger_current__alpha))) / ((pow(var_sodium_calcium_exchanger_current__Km_Nai, 3) + pow(var_sodium_calcium_exchanger_current__Na_o, 3)) * (var_sodium_calcium_exchanger_current__Km_Ca + var_sodium_calcium_exchanger_current__Ca_o) * (1.0 + (var_sodium_calcium_exchanger_current__K_sat * exp(((var_sodium_calcium_exchanger_current__gamma - 1.0) * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)))));
//	   Real var_membrane__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa;
//	   Real var_reversal_potentials__Ca_o = var_calcium_dynamics__Ca_o;
//	   Real var_reversal_potentials__Ca_i = var_calcium_dynamics__Ca_i;
//	   Real var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
//	   Real var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
//	   const Real var_calcium_background_current__g_bca = 0.000592;
//	   Real var_calcium_background_current__V = var_membrane__V;
//	   Real var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
//	   Real var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
//	   const Real var_potassium_pump_current__g_pK = 0.0146;
//	   Real var_potassium_pump_current__V = var_membrane__V;
//	   Real var_potassium_pump_current__E_K = var_reversal_potentials__E_K;
//	   Real var_potassium_pump_current__i_p_K = (var_potassium_pump_current__g_pK * (var_potassium_pump_current__V - var_potassium_pump_current__E_K)) / (1.0 + exp((25.0 - var_potassium_pump_current__V) / 5.98));
//	   Real var_membrane__i_p_K = var_potassium_pump_current__i_p_K;
//	   const Real var_calcium_pump_current__K_pCa = 0.0005;
//	   Real var_calcium_pump_current__Ca_i = var_calcium_dynamics__Ca_i;
//	   const Real var_calcium_pump_current__g_pCa = 0.1238;
//	   Real var_calcium_pump_current__i_p_Ca = (var_calcium_pump_current__g_pCa * var_calcium_pump_current__Ca_i) / (var_calcium_pump_current__Ca_i + var_calcium_pump_current__K_pCa);
//	   Real var_membrane__i_p_Ca = var_calcium_pump_current__i_p_Ca;
//	   const Real var_membrane__stim_end = 100000.0;
//	   //const Real var_membrane__stim_amplitude =  -52.0;
//	   const Real var_membrane__stim_amplitude =  0.0;
//	   const Real var_membrane__stim_duration = 1.0;
//	   const Real var_membrane__stim_period = 1000.0;
//	   const Real var_membrane__stim_start = 100.0 + (Real)col;
//	   Real var_membrane__i_Stim = 0.0; //((var_membrane__time >= var_membrane__stim_start) && (var_membrane__time <= var_membrane__stim_end) && (((var_membrane__time - var_membrane__stim_start) - (floor((var_membrane__time - var_membrane__stim_start) / var_membrane__stim_period) * var_membrane__stim_period)) <= var_membrane__stim_duration)) ? var_membrane__stim_amplitude : 0.0;
//	   Real var_rapid_time_dependent_potassium_current_Xr1_gate__V = var_rapid_time_dependent_potassium_current__V;
//	   Real var_rapid_time_dependent_potassium_current_Xr1_gate__xr1_inf = 1.0 / (1.0 + exp(((-26.0) - var_rapid_time_dependent_potassium_current_Xr1_gate__V) / 7.0));
//	   Real var_rapid_time_dependent_potassium_current_Xr1_gate__alpha_xr1 = 450.0 / (1.0 + exp(((-45.0) - var_rapid_time_dependent_potassium_current_Xr1_gate__V) / 10.0));
//	   Real var_rapid_time_dependent_potassium_current_Xr1_gate__beta_xr1 = 6.0 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr1_gate__V + 30.0) / 11.5));
//	   Real var_rapid_time_dependent_potassium_current_Xr1_gate__tau_xr1 = 1.0 * var_rapid_time_dependent_potassium_current_Xr1_gate__alpha_xr1 * var_rapid_time_dependent_potassium_current_Xr1_gate__beta_xr1;
//	   Real var_rapid_time_dependent_potassium_current_Xr2_gate__V = var_rapid_time_dependent_potassium_current__V;
//	   Real var_rapid_time_dependent_potassium_current_Xr2_gate__xr2_inf = 1.0 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr2_gate__V + 88.0) / 24.0));
//	   Real var_rapid_time_dependent_potassium_current_Xr2_gate__alpha_xr2 = 3.0 / (1.0 + exp(((-60.0) - var_rapid_time_dependent_potassium_current_Xr2_gate__V) / 20.0));
//	   Real var_rapid_time_dependent_potassium_current_Xr2_gate__beta_xr2 = 1.12 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr2_gate__V - 60.0) / 20.0));
//	   Real var_rapid_time_dependent_potassium_current_Xr2_gate__tau_xr2 = 1.0 * var_rapid_time_dependent_potassium_current_Xr2_gate__alpha_xr2 * var_rapid_time_dependent_potassium_current_Xr2_gate__beta_xr2;
//	   Real var_slow_time_dependent_potassium_current_Xs_gate__V = var_slow_time_dependent_potassium_current__V;
//	   Real var_slow_time_dependent_potassium_current_Xs_gate__xs_inf = 1.0 / (1.0 + exp(((-5.0) - var_slow_time_dependent_potassium_current_Xs_gate__V) / 14.0));
//	   Real var_slow_time_dependent_potassium_current_Xs_gate__alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - var_slow_time_dependent_potassium_current_Xs_gate__V) / 6.0));
//	   Real var_slow_time_dependent_potassium_current_Xs_gate__beta_xs = 1.0 / (1.0 + exp((var_slow_time_dependent_potassium_current_Xs_gate__V - 35.0) / 15.0));
//	   Real var_slow_time_dependent_potassium_current_Xs_gate__tau_xs = (1.0 * var_slow_time_dependent_potassium_current_Xs_gate__alpha_xs * var_slow_time_dependent_potassium_current_Xs_gate__beta_xs) + 80.0;
//	   Real var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
//	   Real var_fast_sodium_current_m_gate__m_inf = 1.0 / pow(1.0 + exp(((-56.86) - var_fast_sodium_current_m_gate__V) / 9.03), 2);
//	   Real var_fast_sodium_current_m_gate__alpha_m = 1.0 / (1.0 + exp(((-60.0) - var_fast_sodium_current_m_gate__V) / 5.0));
//	   Real var_fast_sodium_current_m_gate__beta_m = (0.1 / (1.0 + exp((var_fast_sodium_current_m_gate__V + 35.0) / 5.0))) + (0.1 / (1.0 + exp((var_fast_sodium_current_m_gate__V - 50.0) / 200.0)));
//	   Real var_fast_sodium_current_m_gate__tau_m = 1.0 * var_fast_sodium_current_m_gate__alpha_m * var_fast_sodium_current_m_gate__beta_m;
//	   Real var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
//	   Real var_fast_sodium_current_h_gate__h_inf = 1.0 / pow(1.0 + exp((var_fast_sodium_current_h_gate__V + 71.55) / 7.43), 2);
//	   Real var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.057 * exp((-(var_fast_sodium_current_h_gate__V + 80.0)) / 6.8)) : 0.0;
//	   Real var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((2.7 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.3485 * var_fast_sodium_current_h_gate__V))) : (0.77 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1)))));
//	   Real var_fast_sodium_current_h_gate__tau_h = 1.0 / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h);
//	   Real var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V;
//	   Real var_fast_sodium_current_j_gate__j_inf = 1.0 / pow(1.0 + exp((var_fast_sodium_current_j_gate__V + 71.55) / 7.43), 2);
//	   Real var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((((((-25428.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (6.948e-06 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / 1.0) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0;
//	   Real var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.02424 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.6 * exp(0.057 * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))));
//	   Real var_fast_sodium_current_j_gate__tau_j = 1.0 / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j);
//	   Real var_L_type_Ca_current_d_gate__V = var_L_type_Ca_current__V;
//	   Real var_L_type_Ca_current_d_gate__d_inf = 1.0 / (1.0 + exp(((-8.0) - var_L_type_Ca_current_d_gate__V) / 7.5));
//	   Real var_L_type_Ca_current_d_gate__alpha_d = (1.4 / (1.0 + exp(((-35.0) - var_L_type_Ca_current_d_gate__V) / 13.0))) + 0.25;
//	   Real var_L_type_Ca_current_d_gate__beta_d = 1.4 / (1.0 + exp((var_L_type_Ca_current_d_gate__V + 5.0) / 5.0));
//	   Real var_L_type_Ca_current_d_gate__gamma_d = 1.0 / (1.0 + exp((50.0 - var_L_type_Ca_current_d_gate__V) / 20.0));
//	   Real var_L_type_Ca_current_d_gate__tau_d = (1.0 * var_L_type_Ca_current_d_gate__alpha_d * var_L_type_Ca_current_d_gate__beta_d) + var_L_type_Ca_current_d_gate__gamma_d;
//	   Real var_L_type_Ca_current_f_gate__V = var_L_type_Ca_current__V;
//	   Real var_L_type_Ca_current_f_gate__f_inf = 1.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 20.0) / 7.0));
//	   Real var_L_type_Ca_current_f_gate__tau_f = (1102.5 * exp((-pow(var_L_type_Ca_current_f_gate__V + 27.0, 2)) / 225.0)) + (200.0 / (1.0 + exp((13.0 - var_L_type_Ca_current_f_gate__V) / 10.0))) + (180.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 30.0) / 10.0))) + 20.0;
//	   Real var_L_type_Ca_current_f2_gate__V = var_L_type_Ca_current__V;
//	   Real var_L_type_Ca_current_f2_gate__f2_inf = (0.67 / (1.0 + exp((var_L_type_Ca_current_f2_gate__V + 35.0) / 7.0))) + 0.33;
//	   Real var_L_type_Ca_current_f2_gate__tau_f2 = (562.0 * exp((-pow(var_L_type_Ca_current_f2_gate__V + 27.0, 2)) / 240.0)) + (31.0 / (1.0 + exp((25.0 - var_L_type_Ca_current_f2_gate__V) / 10.0))) + (80.0 / (1.0 + exp((var_L_type_Ca_current_f2_gate__V + 30.0) / 10.0)));
//	   Real var_L_type_Ca_current_fCass_gate__Ca_ss = var_L_type_Ca_current__Ca_ss;
//	   Real var_L_type_Ca_current_fCass_gate__fCass_inf = (0.6 / (1.0 + pow(var_L_type_Ca_current_fCass_gate__Ca_ss / 0.05, 2))) + 0.4;
//	   Real var_L_type_Ca_current_fCass_gate__tau_fCass = (80.0 / (1.0 + pow(var_L_type_Ca_current_fCass_gate__Ca_ss / 0.05, 2))) + 2.0;
//	   Real var_transient_outward_current_s_gate__V = var_transient_outward_current__V;
//	   Real var_transient_outward_current_s_gate__s_inf = 1.0 / (1.0 + exp((var_transient_outward_current_s_gate__V + 20.0) / 5.0));
//	   Real var_transient_outward_current_s_gate__tau_s = (85.0 * exp((-pow(var_transient_outward_current_s_gate__V + 45.0, 2)) / 320.0)) + (5.0 / (1.0 + exp((var_transient_outward_current_s_gate__V - 20.0) / 5.0))) + 3.0;
//	   Real var_transient_outward_current_r_gate__V = var_transient_outward_current__V;
//	   Real var_transient_outward_current_r_gate__r_inf = 1.0 / (1.0 + exp((20.0 - var_transient_outward_current_r_gate__V) / 6.0));
//	   Real var_transient_outward_current_r_gate__tau_r = (9.5 * exp((-pow(var_transient_outward_current_r_gate__V + 40.0, 2)) / 1800.0)) + 0.8;
//	   const Real var_calcium_dynamics__V_rel = 0.102;
//	   const Real var_calcium_dynamics__k1_prime = 0.15;
//	   const Real var_calcium_dynamics__max_sr = 2.5;
//	   const Real var_calcium_dynamics__EC = 1.5;
//	   const Real var_calcium_dynamics__min_sr = 1.0;
//	   Real var_calcium_dynamics__kcasr = var_calcium_dynamics__max_sr - ((var_calcium_dynamics__max_sr - var_calcium_dynamics__min_sr) / (1.0 + pow(var_calcium_dynamics__EC / var_calcium_dynamics__Ca_SR, 2)));
//	   Real var_calcium_dynamics__k1 = var_calcium_dynamics__k1_prime / var_calcium_dynamics__kcasr;
//	   const Real var_calcium_dynamics__k3 = 0.06;
//	   Real var_calcium_dynamics__O = (var_calcium_dynamics__k1 * pow(var_calcium_dynamics__Ca_ss, 2) * var_calcium_dynamics__R_prime) / (var_calcium_dynamics__k3 + (var_calcium_dynamics__k1 * pow(var_calcium_dynamics__Ca_ss, 2)));
//	   Real var_calcium_dynamics__i_rel = var_calcium_dynamics__V_rel * var_calcium_dynamics__O * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_ss);
//	   const Real var_calcium_dynamics__Vmax_up = 0.006375;
//	   const Real var_calcium_dynamics__K_up = 0.00025;
//	   Real var_calcium_dynamics__i_up = var_calcium_dynamics__Vmax_up / (1.0 + (pow(var_calcium_dynamics__K_up, 2) / pow(var_calcium_dynamics__Ca_i, 2)));
//	   const Real var_calcium_dynamics__V_leak = 0.00036;
//	   Real var_calcium_dynamics__i_leak = var_calcium_dynamics__V_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i);
//	   const Real var_calcium_dynamics__V_xfer = 0.0038;
//	   Real var_calcium_dynamics__i_xfer = var_calcium_dynamics__V_xfer * (var_calcium_dynamics__Ca_ss - var_calcium_dynamics__Ca_i);
//	   const Real var_calcium_dynamics__k2_prime = 0.045;
//	   Real var_calcium_dynamics__k2 = var_calcium_dynamics__k2_prime * var_calcium_dynamics__kcasr;
//	   const Real var_calcium_dynamics__k4 = 0.005;
//	   const Real var_calcium_dynamics__Buf_c = 0.2;
//	   const Real var_calcium_dynamics__K_buf_c = 0.001;
//	   Real var_calcium_dynamics__Ca_i_bufc = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_c * var_calcium_dynamics__K_buf_c) / pow(var_calcium_dynamics__Ca_i + var_calcium_dynamics__K_buf_c, 2)));
//	   const Real var_calcium_dynamics__K_buf_sr = 0.3;
//	   const Real var_calcium_dynamics__Buf_sr = 10.0;
//	   Real var_calcium_dynamics__Ca_sr_bufsr = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_sr * var_calcium_dynamics__K_buf_sr) / pow(var_calcium_dynamics__Ca_SR + var_calcium_dynamics__K_buf_sr, 2)));
//	   const Real var_calcium_dynamics__Buf_ss = 0.4;
//	   const Real var_calcium_dynamics__K_buf_ss = 0.00025;
//	   Real var_calcium_dynamics__Ca_ss_bufss = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_ss * var_calcium_dynamics__K_buf_ss) / pow(var_calcium_dynamics__Ca_ss + var_calcium_dynamics__K_buf_ss, 2)));
//	   const Real var_calcium_dynamics__V_sr = 0.001094;
//	   const Real var_calcium_dynamics__V_ss = 5.468e-05;
//	   Real var_calcium_dynamics__V_c = var_membrane__V_c;
//	   Real var_calcium_dynamics__F = var_membrane__F;
//	   Real var_calcium_dynamics__Cm = var_membrane__Cm;
//	   Real var_calcium_dynamics__i_CaL = var_L_type_Ca_current__i_CaL;
//	   Real var_calcium_dynamics__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa;
//	   Real var_calcium_dynamics__i_p_Ca = var_calcium_pump_current__i_p_Ca;
//	   Real var_calcium_dynamics__i_b_Ca = var_calcium_background_current__i_b_Ca;
//	   Real var_sodium_dynamics__F = var_membrane__F;
//	   Real var_sodium_dynamics__Cm = var_membrane__Cm;
//	   Real var_sodium_dynamics__V_c = var_membrane__V_c;
//	   Real var_sodium_dynamics__i_Na = var_fast_sodium_current__i_Na;
//	   Real var_sodium_dynamics__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa;
//	   Real var_sodium_dynamics__i_NaK = var_sodium_potassium_pump_current__i_NaK;
//	   Real var_sodium_dynamics__i_b_Na = var_sodium_background_current__i_b_Na;
//	   Real var_potassium_dynamics__F = var_membrane__F;
//	   Real var_potassium_dynamics__Cm = var_membrane__Cm;
//	   Real var_potassium_dynamics__V_c = var_membrane__V_c;
//	   Real var_potassium_dynamics__i_K1 = var_inward_rectifier_potassium_current__i_K1;
//	   Real var_potassium_dynamics__i_to = var_transient_outward_current__i_to;
//	   Real var_potassium_dynamics__i_NaK = var_sodium_potassium_pump_current__i_NaK;
//	   Real var_potassium_dynamics__i_Kr = var_rapid_time_dependent_potassium_current__i_Kr;
//	   Real var_potassium_dynamics__i_Ks = var_slow_time_dependent_potassium_current__i_Ks;
//	   Real var_potassium_dynamics__i_p_K = var_potassium_pump_current__i_p_K;
//	   Real var_potassium_dynamics__i_Stim = var_membrane__i_Stim;
//	   Real d_dt_membrane__V = -(var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_CaL + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_NaCa + var_membrane__i_b_Ca + var_membrane__i_p_K + var_membrane__i_p_Ca + var_membrane__i_Stim);
//	   Real d_dt_rapid_time_dependent_potassium_current_Xr1_gate__Xr1 = (var_rapid_time_dependent_potassium_current_Xr1_gate__xr1_inf - var_rapid_time_dependent_potassium_current_Xr1_gate__Xr1) / var_rapid_time_dependent_potassium_current_Xr1_gate__tau_xr1;
//	   Real d_dt_rapid_time_dependent_potassium_current_Xr2_gate__Xr2 = (var_rapid_time_dependent_potassium_current_Xr2_gate__xr2_inf - var_rapid_time_dependent_potassium_current_Xr2_gate__Xr2) / var_rapid_time_dependent_potassium_current_Xr2_gate__tau_xr2;
//	   Real d_dt_slow_time_dependent_potassium_current_Xs_gate__Xs = (var_slow_time_dependent_potassium_current_Xs_gate__xs_inf - var_slow_time_dependent_potassium_current_Xs_gate__Xs) / var_slow_time_dependent_potassium_current_Xs_gate__tau_xs;
//	   Real d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__m_inf - var_fast_sodium_current_m_gate__m) / var_fast_sodium_current_m_gate__tau_m;
//	   Real d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__h_inf - var_fast_sodium_current_h_gate__h) / var_fast_sodium_current_h_gate__tau_h;
//	   Real d_dt_fast_sodium_current_j_gate__j = (var_fast_sodium_current_j_gate__j_inf - var_fast_sodium_current_j_gate__j) / var_fast_sodium_current_j_gate__tau_j;
//	   Real d_dt_L_type_Ca_current_d_gate__d = (var_L_type_Ca_current_d_gate__d_inf - var_L_type_Ca_current_d_gate__d) / var_L_type_Ca_current_d_gate__tau_d;
//	   Real d_dt_L_type_Ca_current_f_gate__f = (var_L_type_Ca_current_f_gate__f_inf - var_L_type_Ca_current_f_gate__f) / var_L_type_Ca_current_f_gate__tau_f;
//	   Real d_dt_L_type_Ca_current_f2_gate__f2 = (var_L_type_Ca_current_f2_gate__f2_inf - var_L_type_Ca_current_f2_gate__f2) / var_L_type_Ca_current_f2_gate__tau_f2;
//	   Real d_dt_L_type_Ca_current_fCass_gate__fCass = (var_L_type_Ca_current_fCass_gate__fCass_inf - var_L_type_Ca_current_fCass_gate__fCass) / var_L_type_Ca_current_fCass_gate__tau_fCass;
//	   Real d_dt_transient_outward_current_s_gate__s = (var_transient_outward_current_s_gate__s_inf - var_transient_outward_current_s_gate__s) / var_transient_outward_current_s_gate__tau_s;
//	   Real d_dt_transient_outward_current_r_gate__r = (var_transient_outward_current_r_gate__r_inf - var_transient_outward_current_r_gate__r) / var_transient_outward_current_r_gate__tau_r;
//	   Real d_dt_calcium_dynamics__R_prime = ((-var_calcium_dynamics__k2) * var_calcium_dynamics__Ca_ss * var_calcium_dynamics__R_prime) + (var_calcium_dynamics__k4 * (1.0 - var_calcium_dynamics__R_prime));
//	   Real d_dt_calcium_dynamics__Ca_i = var_calcium_dynamics__Ca_i_bufc * (((((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_c) + var_calcium_dynamics__i_xfer) - ((((var_calcium_dynamics__i_b_Ca + var_calcium_dynamics__i_p_Ca) - (2.0 * var_calcium_dynamics__i_NaCa)) * var_calcium_dynamics__Cm) / (2.0 * var_calcium_dynamics__V_c * var_calcium_dynamics__F)));
//	   Real d_dt_calcium_dynamics__Ca_SR = var_calcium_dynamics__Ca_sr_bufsr * (var_calcium_dynamics__i_up - (var_calcium_dynamics__i_rel + var_calcium_dynamics__i_leak));
//	   Real d_dt_calcium_dynamics__Ca_ss = var_calcium_dynamics__Ca_ss_bufss * (((((-var_calcium_dynamics__i_CaL) * var_calcium_dynamics__Cm) / (2.0 * var_calcium_dynamics__V_ss * var_calcium_dynamics__F)) + ((var_calcium_dynamics__i_rel * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_ss)) - ((var_calcium_dynamics__i_xfer * var_calcium_dynamics__V_c) / var_calcium_dynamics__V_ss));
//	   Real d_dt_sodium_dynamics__Na_i = ((-(var_sodium_dynamics__i_Na + var_sodium_dynamics__i_b_Na + (3.0 * var_sodium_dynamics__i_NaK) + (3.0 * var_sodium_dynamics__i_NaCa))) / (var_sodium_dynamics__V_c * var_sodium_dynamics__F)) * var_sodium_dynamics__Cm;
//	   Real d_dt_potassium_dynamics__K_i = ((-((var_potassium_dynamics__i_K1 + var_potassium_dynamics__i_to + var_potassium_dynamics__i_Kr + var_potassium_dynamics__i_Ks + var_potassium_dynamics__i_p_K + var_potassium_dynamics__i_Stim) - (2.0 * var_potassium_dynamics__i_NaK))) / (var_potassium_dynamics__V_c * var_potassium_dynamics__F)) * var_potassium_dynamics__Cm;
//
//	   rDY[ 0] = d_dt_membrane__V;
//	   rDY[ 1] = d_dt_rapid_time_dependent_potassium_current_Xr1_gate__Xr1;
//	   rDY[ 2] = d_dt_rapid_time_dependent_potassium_current_Xr2_gate__Xr2;
//	   rDY[ 3] = d_dt_slow_time_dependent_potassium_current_Xs_gate__Xs;
//	   rDY[ 4] = d_dt_fast_sodium_current_m_gate__m;
//	   rDY[ 5] = d_dt_fast_sodium_current_h_gate__h;
//	   rDY[ 6] = d_dt_fast_sodium_current_j_gate__j;
//	   rDY[ 7] = d_dt_L_type_Ca_current_d_gate__d;
//	   rDY[ 8] = d_dt_L_type_Ca_current_f_gate__f;
//	   rDY[ 9] = d_dt_L_type_Ca_current_f2_gate__f2;
//	   rDY[10] = d_dt_L_type_Ca_current_fCass_gate__fCass;
//	   rDY[11] = d_dt_transient_outward_current_s_gate__s;
//	   rDY[12] = d_dt_transient_outward_current_r_gate__r;
//	   rDY[13] = d_dt_calcium_dynamics__Ca_i;
//	   rDY[14] = d_dt_calcium_dynamics__Ca_SR;
//	   rDY[15] = d_dt_calcium_dynamics__Ca_ss;
//	   rDY[16] = d_dt_calcium_dynamics__R_prime;
//	   rDY[17] = d_dt_sodium_dynamics__Na_i;
//	   rDY[18] = d_dt_potassium_dynamics__K_i;

}
