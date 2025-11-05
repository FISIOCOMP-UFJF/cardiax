#include "bondarenko.h"
#include <stdio.h>

void setIC_ode_cpu(Real *sv, int cellID, int NEQ = 40)
{
  // Vector containing the initial conditions of the cell
  Real IC[NEQ];
  
 	IC[0] = -82.4202;	 // V millivolt
	IC[1] = 0.115001;	 // Cai micromolar
	IC[2] = 0.115001;	 // Cass micromolar
	IC[3] = 1299.5;	 // CaJSR micromolar
	IC[4] = 1299.5;	 // CaNSR micromolar
	IC[5] = 0;	 // P_RyR dimensionless
	IC[6] = 11.2684;	 // LTRPN_Ca micromolar
	IC[7] = 125.29;	 // HTRPN_Ca micromolar
	IC[8] = 0.149102e-4;	 // P_O1 dimensionless
	IC[9] = 0.951726e-10;	 // P_O2 dimensionless
	IC[10] = 0.16774e-3;	 // P_C2 dimensionless
	IC[11] = 0.930308e-18;	 // O dimensionless
	IC[12] = 0.124216e-3;	 // C2 dimensionless
	IC[13] = 0.578679e-8;	 // C3 dimensionless
	IC[14] = 0.119816e-12;	 // C4 dimensionless
	IC[15] = 0.497923e-18;	 // I1 dimensionless
	IC[16] = 0.345847e-13;	 // I2 dimensionless
	IC[17] = 0.185106e-13;	 // I3 dimensionless
	IC[18] = 14237.1;	 // Nai micromolar
	IC[19] = 0.020752;	 // C_Na2 dimensionless
	IC[20] = 0.279132e-3;	 // C_Na1 dimensionless
	IC[21] = 0.713483e-6;	 // O_Na dimensionless
	IC[22] = 0.153176e-3;	 // IF_Na dimensionless
	IC[23] = 0.673345e-6;	 // I1_Na dimensionless
	IC[24] = 0.155787e-8;	 // I2_Na dimensionless
	IC[25] = 0.0113879;	 // IC_Na2 dimensionless
	IC[26] = 0.34278;	 // IC_Na3 dimensionless
	IC[27] = 143720;	 // Ki micromolar
	IC[28] = 0.265563e-2;	 // ato_f dimensionless
	IC[29] = 0.999977;	 // ito_f dimensionless
	IC[30] = 0.417069e-3;	 // ato_s dimensionless
	IC[31] = 0.998543;	 // ito_s dimensionless
	IC[32] = 0.262753e-3;	 // nKs dimensionless
	IC[33] = 0.417069e-3;	 // aur dimensionless
	IC[34] = 0.998543;	 // iur dimensionless
	IC[35] = 0.417069e-3;	 // aKss dimensionless
	IC[36] = 1;	 // iKss dimensionless
	IC[37] = 0.641229e-3;	 // C_K2 dimensionless
	IC[38] = 0.992513e-3;	 // C_K1 dimensionless
	IC[39] = 0.175298e-3;	 // O_K dimensionless
	IC[40] = 0.319129e-4;	 // I_K dimensionless

  for(int i = 0; i < NEQ; i++)
    sv[cellID+i] = IC[i];
}

void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, 
		   int NEQ, int col)
{
   solve_Forward_Euler_cpu(time, sv, dt, NEQ, cellID, col);
}

void solve_Forward_Euler_cpu(Real time, Real *sv, Real dt, 
			     int NEQ, int cellID, int col)
{
   Real rY[NEQ], rDY[NEQ];

   for(int i = 0; i < NEQ; i++) rY[i] = sv[cellID+i];

   RHS_Bondarenko_cpu(time, rY, rDY, col);

   for(int i = 0; i < NEQ; i++) sv[cellID+i] = dt*rDY[i] + rY[i];
}

void RHS_Bondarenko_cpu(Real time_, Real *sv_, Real *rDY_, int col)
{

	// State variables
	const Real V_old_ = sv_[0];	 // initial value = -82.4202 millivolt
	const Real Cai_old_ = sv_[1];	 // initial value = 0.115001 micromolar
	const Real Cass_old_ = sv_[2];	 // initial value = 0.115001 micromolar
	const Real CaJSR_old_ = sv_[3];	 // initial value = 1299.5 micromolar
	const Real CaNSR_old_ = sv_[4];	 // initial value = 1299.5 micromolar
	const Real P_RyR_old_ = sv_[5];	 // initial value = 0 dimensionless
	const Real LTRPN_Ca_old_ = sv_[6];	 // initial value = 11.2684 micromolar
	const Real HTRPN_Ca_old_ = sv_[7];	 // initial value = 125.29 micromolar
	const Real P_O1_old_ = sv_[8];	 // initial value = 0.149102e-4 dimensionless
	const Real P_O2_old_ = sv_[9];	 // initial value = 0.951726e-10 dimensionless
	const Real P_C2_old_ = sv_[10];	 // initial value = 0.16774e-3 dimensionless
	const Real O_old_ = sv_[11];	 // initial value = 0.930308e-18 dimensionless
	const Real C2_old_ = sv_[12];	 // initial value = 0.124216e-3 dimensionless
	const Real C3_old_ = sv_[13];	 // initial value = 0.578679e-8 dimensionless
	const Real C4_old_ = sv_[14];	 // initial value = 0.119816e-12 dimensionless
	const Real I1_old_ = sv_[15];	 // initial value = 0.497923e-18 dimensionless
	const Real I2_old_ = sv_[16];	 // initial value = 0.345847e-13 dimensionless
	const Real I3_old_ = sv_[17];	 // initial value = 0.185106e-13 dimensionless
	const Real Nai_old_ = sv_[18];	 // initial value = 14237.1 micromolar
	const Real C_Na2_old_ = sv_[19];	 // initial value = 0.020752 dimensionless
	const Real C_Na1_old_ = sv_[20];	 // initial value = 0.279132e-3 dimensionless
	const Real O_Na_old_ = sv_[21];	 // initial value = 0.713483e-6 dimensionless
	const Real IF_Na_old_ = sv_[22];	 // initial value = 0.153176e-3 dimensionless
	const Real I1_Na_old_ = sv_[23];	 // initial value = 0.673345e-6 dimensionless
	const Real I2_Na_old_ = sv_[24];	 // initial value = 0.155787e-8 dimensionless
	const Real IC_Na2_old_ = sv_[25];	 // initial value = 0.0113879 dimensionless
	const Real IC_Na3_old_ = sv_[26];	 // initial value = 0.34278 dimensionless
	const Real Ki_old_ = sv_[27];	 // initial value = 143720 micromolar
	const Real ato_f_old_ = sv_[28];	 // initial value = 0.265563e-2 dimensionless
	const Real ito_f_old_ = sv_[29];	 // initial value = 0.999977 dimensionless
	const Real ato_s_old_ = sv_[30];	 // initial value = 0.417069e-3 dimensionless
	const Real ito_s_old_ = sv_[31];	 // initial value = 0.998543 dimensionless
	const Real nKs_old_ = sv_[32];	 // initial value = 0.262753e-3 dimensionless
	const Real aur_old_ = sv_[33];	 // initial value = 0.417069e-3 dimensionless
	const Real iur_old_ = sv_[34];	 // initial value = 0.998543 dimensionless
	const Real aKss_old_ = sv_[35];	 // initial value = 0.417069e-3 dimensionless
	const Real iKss_old_ = sv_[36];	 // initial value = 1 dimensionless
	const Real C_K2_old_ = sv_[37];	 // initial value = 0.641229e-3 dimensionless
	const Real C_K1_old_ = sv_[38];	 // initial value = 0.992513e-3 dimensionless
	const Real O_K_old_ = sv_[39];	 // initial value = 0.175298e-3 dimensionless
	const Real I_K_old_ = sv_[40];	 // initial value = 0.319129e-4 dimensionless


	// Parameters
	const Real stim_amplitude = 0.0;	 // picoA_per_picoF
	const Real stim_start = 20;	 // millisecond
	const Real stim_end = 100000;	 // millisecond
	const Real stim_period = 71.43;	 // millisecond
	const Real stim_duration = 0.5;	 // millisecond
	const Real Acap = 1.534e-4;	 // cm2
	const Real Cm = 1;	 // microF_per_cm2
	const Real Vmyo = 25.84e-6;	 // microlitre
	const Real F = 96.5;	 // coulomb_per_millimole
	const Real VJSR = 0.12e-6;	 // microlitre
	const Real Vss = 1.485e-9;	 // microlitre
	const Real VNSR = 2.098e-6;	 // microlitre
	const Real CMDN_tot = 50;	 // micromolar
	const Real Km_CMDN = 0.238;	 // micromolar
	const Real CSQN_tot = 15000;	 // micromolar
	const Real Km_CSQN = 800;	 // micromolar
	const Real v1 = 4.5;	 // per_millisecond
	const Real tau_tr = 20;	 // millisecond
	const Real tau_xfer = 8;	 // millisecond
	const Real v2 = 1.74e-5;	 // per_millisecond
	const Real v3 = 0.45;	 // micromolar_per_millisecond
	const Real Km_up = 0.5;	 // micromolar
	const Real k_plus_htrpn = 0.00237;	 // per_micromolar_millisecond
	const Real HTRPN_tot = 140;	 // micromolar
	const Real k_plus_ltrpn = 0.0327;	 // per_micromolar_millisecond
	const Real LTRPN_tot = 70;	 // micromolar
	const Real k_minus_htrpn = 3.2e-5;	 // per_millisecond
	const Real k_minus_ltrpn = 0.0196;	 // per_millisecond
	const Real i_CaL_max = 7;	 // picoA_per_picoF
	const Real k_plus_a = 0.006075;	 // micromolar4_per_millisecond
	const Real n = 4;	 // dimensionless
	const Real k_minus_b = 0.965;	 // per_millisecond
	const Real k_minus_c = 0.0008;	 // per_millisecond
	const Real k_minus_a = 0.07125;	 // per_millisecond
	const Real k_plus_b = 0.00405;	 // micromolar3_per_millisecond
	const Real m = 3;	 // dimensionless
	const Real k_plus_c = 0.009;	 // per_millisecond
	const Real g_CaL = 0.1729;	 // milliS_per_microF
	const Real E_CaL = 63;	 // millivolt
	const Real Kpcb = 0.0005;	 // per_millisecond
	const Real Kpc_max = 0.23324;	 // per_millisecond
	const Real Kpc_half = 20;	 // micromolar
	const Real i_pCa_max = 1;	 // picoA_per_picoF
	const Real Km_pCa = 0.5;	 // micromolar
	const Real k_NaCa = 292.8;	 // picoA_per_picoF
	const Real K_mNa = 87500;	 // micromolar
	const Real Nao = 140000;	 // micromolar
	const Real K_mCa = 1380;	 // micromolar
	const Real Cao = 1800;	 // micromolar
	const Real k_sat = 0.1;	 // dimensionless
	const Real eta = 0.35;	 // dimensionless
	const Real R = 8.314;	 // joule_per_mole_kelvin
	const Real T = 298;	 // kelvin
	const Real g_Cab = 0.000367;	 // milliS_per_microF
	const Real g_Na = 13;	 // milliS_per_microF
	const Real Ko = 5400;	 // micromolar
	const Real g_Nab = 0.0026;	 // milliS_per_microF
	const Real g_Kto_f = 0.4067;	 // milliS_per_microF
	const Real g_Kto_s = 0;	 // milliS_per_microF
	const Real g_Ks = 0.00575;	 // milliS_per_microF
	const Real g_Kur = 0.16;	 // milliS_per_microF
	const Real g_Kss = 0.05;	 // milliS_per_microF
	const Real g_Kr = 0.078;	 // milliS_per_microF
	const Real kf = 0.023761;	 // per_millisecond
	const Real kb = 0.036778;	 // per_millisecond
	const Real i_NaK_max = 0.88;	 // picoA_per_picoF
	const Real Km_Nai = 21000;	 // micromolar
	const Real Km_Ko = 1500;	 // micromolar
	const Real g_ClCa = 10;	 // milliS_per_microF
	const Real Km_Cl = 10;	 // micromolar
	const Real E_Cl = -40;	 // millivolt

	// Independent Variable
	Real time_new = time_;


	// Algebraic Equations
	Real calc_i_stim = (((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration)))? stim_amplitude: 0.0;	//0
	Real calc_Bi = pow((1.0+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cai_old_),2.0))),(-1.0));	//6
	Real calc_Bss = pow((1.0+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cass_old_),2.0))),(-1.0));	//7
	Real calc_BJSR = pow((1.0+((CSQN_tot*Km_CSQN)/pow((Km_CSQN+CaJSR_old_),2.0))),(-1.0));	//8
	Real calc_J_rel = (v1*(P_O1_old_+P_O2_old_)*(CaJSR_old_-Cass_old_)*P_RyR_old_);	//9
	Real calc_J_tr = ((CaNSR_old_-CaJSR_old_)/tau_tr);	//10
	Real calc_J_xfer = ((Cass_old_-Cai_old_)/tau_xfer);	//11
	Real calc_J_leak = (v2*(CaNSR_old_-Cai_old_));	//12
	Real calc_J_up = ((v3*pow(Cai_old_,2.0))/(pow(Km_up,2.0)+pow(Cai_old_,2.0)));	//13
	Real calc_J_trpn = (((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))+(k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_)))-((k_minus_htrpn*HTRPN_Ca_old_)+(k_minus_ltrpn*LTRPN_Ca_old_)));	//14
	Real calc_P_C1 = (1.0-(P_C2_old_+P_O1_old_+P_O2_old_));	//19
	Real calc_i_CaL = (g_CaL*O_old_*(V_old_-E_CaL));	//22
	Real calc_C1 = (1.0-(O_old_+C2_old_+C3_old_+C4_old_+I1_old_+I2_old_+I3_old_));	//24
	Real calc_alpha = ((0.4*exp(((V_old_+12.0)/10.0))*((1.0+(0.7*exp(((-pow((V_old_+40.0),2.0))/10.0))))-(0.75*exp(((-pow((V_old_+20.0),2.0))/400.0)))))/(1.0+(0.12*exp(((V_old_+12.0)/10.0)))));	//31
	Real calc_beta = (0.05*exp(((-(V_old_+12.0))/13.0)));	//32
	Real calc_gamma = ((Kpc_max*Cass_old_)/(Kpc_half+Cass_old_));	//33
	Real calc_Kpcf = (13.0*(1.0-exp(((-pow((V_old_+14.5),2.0))/100.0))));	//34
	Real calc_i_pCa = ((i_pCa_max*pow(Cai_old_,2.0))/(pow(Km_pCa,2.0)+pow(Cai_old_,2.0)));	//35
	Real calc_i_NaCa = (((((((k_NaCa*1.0)/(pow(K_mNa,3.0)+pow(Nao,3.0)))*1.0)/(K_mCa+Cao))*1.0)/(1.0+(k_sat*exp((((eta-1.0)*V_old_*F)/(R*T))))))*((exp(((eta*V_old_*F)/(R*T)))*pow(Nai_old_,3.0)*Cao)-(exp((((eta-1.0)*V_old_*F)/(R*T)))*pow(Nao,3.0)*Cai_old_)));	//36
	Real calc_E_CaN = (((R*T)/(2.0*F))*log((Cao/Cai_old_)));	//38
	Real calc_E_Na = (((R*T)/F)*log((((0.9*Nao)+(0.1*Ko))/((0.9*Nai_old_)+(0.1*Ki_old_)))));	//41
	Real calc_C_Na3 = (1.0-(O_Na_old_+C_Na1_old_+C_Na2_old_+IF_Na_old_+I1_Na_old_+I2_Na_old_+IC_Na2_old_+IC_Na3_old_));	//42
	Real calc_alpha_Na11 = (3.802/((0.1027*exp(((-(V_old_+2.5))/17.0)))+(0.2*exp(((-(V_old_+2.5))/150.0)))));	//51
	Real calc_alpha_Na12 = (3.802/((0.1027*exp(((-(V_old_+2.5))/15.0)))+(0.23*exp(((-(V_old_+2.5))/150.0)))));	//52
	Real calc_alpha_Na13 = (3.802/((0.1027*exp(((-(V_old_+2.5))/12.0)))+(0.25*exp(((-(V_old_+2.5))/150.0)))));	//53
	Real calc_beta_Na11 = (0.1917*exp(((-(V_old_+2.5))/20.3)));	//54
	Real calc_beta_Na12 = (0.2*exp(((-(V_old_-2.5))/20.3)));	//55
	Real calc_beta_Na13 = (0.22*exp(((-(V_old_-7.5))/20.3)));	//56
	Real calc_alpha_Na3 = (7e-7*exp(((-(V_old_+7.0))/7.7)));	//57
	Real calc_beta_Na3 = (0.00854+(0.00002*V_old_));	//58
	Real calc_alpha_Na2 = (1.0/((0.188495*exp(((-(V_old_+7.0))/16.6)))+0.393956));	//59
	Real calc_E_K = (((R*T)/F)*log((Ko/Ki_old_)));	//68
	Real calc_alpha_a = (0.18064*exp((0.03577*(V_old_+30.0))));	//71
	Real calc_beta_a = (0.3956*exp(((-0.06237)*(V_old_+30.0))));	//72
	Real calc_alpha_i = ((0.000152*exp(((-(V_old_+13.5))/7.0)))/((0.067083*exp(((-(V_old_+33.5))/7.0)))+1.0));	//73
	Real calc_beta_i = ((0.00095*exp(((V_old_+33.5)/7.0)))/((0.051335*exp(((V_old_+33.5)/7.0)))+1.0));	//74
	Real calc_ass = (1.0/(1.0+exp(((-(V_old_+22.5))/7.7))));	//78
	Real calc_iss = (1.0/(1.0+exp(((V_old_+45.2)/5.7))));	//79
	Real calc_tau_ta_s = ((0.493*exp(((-0.0629)*V_old_)))+2.058);	//80
	Real calc_tau_ti_s = (270.0+(1050.0/(1.0+exp(((V_old_+45.2)/5.7)))));	//81
	Real calc_alpha_n = (V_old_ != -26.5)?((0.00000481333*(V_old_+26.5))/(1.0-exp(((-0.128)*(V_old_+26.5))))): 0.000037604;   //85
	Real calc_beta_n = (0.0000953333*exp(((-0.038)*(V_old_+26.5))));	//86
	Real calc_tau_aur = ((0.493*exp(((-0.0629)*V_old_)))+2.058);	//90
	Real calc_tau_iur = (1200.0-(170.0/(1.0+exp(((V_old_+45.2)/5.7)))));	//91
	Real calc_tau_Kss = ((39.3*exp(((-0.0862)*V_old_)))+13.17);	//95
	Real calc_i_Kr = (g_Kr*O_K_old_*(V_old_-(((R*T)/F)*log((((0.98*Ko)+(0.02*Nao))/((0.98*Ki_old_)+(0.02*Nai_old_)))))));	//96
	Real calc_C_K0 = (1.0-(C_K1_old_+C_K2_old_+O_K_old_+I_K_old_));	//97
	Real calc_alpha_a0 = (0.022348*exp((0.01176*V_old_)));	//102
	Real calc_beta_a0 = (0.047002*exp(((-0.0631)*V_old_)));	//103
	Real calc_alpha_a1 = (0.013733*exp((0.038198*V_old_)));	//104
	Real calc_beta_a1 = (0.0000689*exp(((-0.04178)*V_old_)));	//105
	Real calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (0.090821*exp((0.023391*(V_old_+5.0))));	//106
	Real calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (0.006497*exp(((-0.03268)*(V_old_+5.0))));	//107
	Real calc_sigma = ((1.0/7.0)*(exp((Nao/67300.0))-1.0));	//110
	Real calc_O_ClCa = (0.2/(1.0+exp(((-(V_old_-46.7))/7.8))));	//112
	Real calc_beta_Na2 = ((calc_alpha_Na13*calc_alpha_Na2*calc_alpha_Na3)/(calc_beta_Na13*calc_beta_Na3));	//60
	Real calc_alpha_Na4 = (calc_alpha_Na2/1000.0);	//61
	Real calc_beta_Na4 = calc_alpha_Na3;	//62
	Real calc_alpha_Na5 = (calc_alpha_Na2/95000.0);	//63
	Real calc_beta_Na5 = (calc_alpha_Na3/50.0);	//64
	Real calc_i_Nab = (g_Nab*(V_old_-calc_E_Na));	//65
	Real calc_i_Kto_s = (g_Kto_s*ato_s_old_*ito_s_old_*(V_old_-calc_E_K));	//75
	Real calc_i_K1 = ((((0.2938*Ko)/(Ko+210.0))*(V_old_-calc_E_K))/(1.0+exp((0.0896*(V_old_-calc_E_K)))));	//82
	Real calc_i_Ks = (g_Ks*pow(nKs_old_,2.0)*(V_old_-calc_E_K));	//83
	Real calc_i_Kur = (g_Kur*aur_old_*iur_old_*(V_old_-calc_E_K));	//87
	Real calc_i_Kss = (g_Kss*aKss_old_*iKss_old_*(V_old_-calc_E_K));	//92
	Real calc_i_Cab = (g_Cab*(V_old_-calc_E_CaN));	//37
	Real calc_i_Na = (g_Na*O_Na_old_*(V_old_-calc_E_Na));	//40
	Real calc_i_Kto_f = (g_Kto_f*pow(ato_f_old_,3.0)*ito_f_old_*(V_old_-calc_E_K));	//67
	Real calc_f_NaK = (1.0/(1.0+(0.1245*exp((((-0.1)*V_old_*F)/(R*T))))+(0.0365*calc_sigma*exp((((-V_old_)*F)/(R*T))))));	//109
	Real calc_i_ClCa = (((g_ClCa*calc_O_ClCa*Cai_old_)/(Cai_old_+Km_Cl))*(V_old_-E_Cl));	//111
	Real calc_i_NaK = ((((i_NaK_max*calc_f_NaK*1.0)/(1.0+pow((Km_Nai/Nai_old_),1.5)))*Ko)/(Ko+Km_Ko));	//108

	// Differential Equations
	Real d_dt_V = (-(calc_i_CaL+calc_i_pCa+calc_i_NaCa+calc_i_Cab+calc_i_Na+calc_i_Nab+calc_i_NaK+calc_i_Kto_f+calc_i_Kto_s+calc_i_K1+calc_i_Ks+calc_i_Kur+calc_i_Kss+calc_i_Kr+calc_i_ClCa+calc_i_stim));	// 1
	Real d_dt_Cai = (calc_Bi*((calc_J_leak+calc_J_xfer)-(calc_J_up+calc_J_trpn+((((calc_i_Cab+calc_i_pCa)-(2.0*calc_i_NaCa))*Acap*Cm)/(2.0*Vmyo*F)))));	// 2
	Real d_dt_Cass = (calc_Bss*(((calc_J_rel*VJSR)/Vss)-(((calc_J_xfer*Vmyo)/Vss)+((calc_i_CaL*Acap*Cm)/(2.0*Vss*F)))));	// 3
	Real d_dt_CaJSR = (calc_BJSR*(calc_J_tr-calc_J_rel));	// 4
	Real d_dt_CaNSR = ((((calc_J_up-calc_J_leak)*Vmyo)/VNSR)-((calc_J_tr*VJSR)/VNSR));	// 5
	Real d_dt_P_RyR = (((-0.04)*P_RyR_old_)-(((0.1*calc_i_CaL)/i_CaL_max)*exp(((-pow((V_old_-5.0),2.0))/648.0))));	// 15
	Real d_dt_LTRPN_Ca = ((k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_))-(k_minus_ltrpn*LTRPN_Ca_old_));	// 16
	Real d_dt_HTRPN_Ca = ((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))-(k_minus_htrpn*HTRPN_Ca_old_));	// 17
	Real d_dt_P_O1 = (((k_plus_a*pow(Cass_old_,n)*calc_P_C1)+(k_minus_b*P_O2_old_)+(k_minus_c*P_C2_old_))-((k_minus_a*P_O1_old_)+(k_plus_b*pow(Cass_old_,m)*P_O1_old_)+(k_plus_c*P_O1_old_)));	// 18
	Real d_dt_P_O2 = ((k_plus_b*pow(Cass_old_,m)*P_O1_old_)-(k_minus_b*P_O2_old_));	// 20
	Real d_dt_P_C2 = ((k_plus_c*P_O1_old_)-(k_minus_c*P_C2_old_));	// 21
	Real d_dt_O = (((calc_alpha*C4_old_)+(Kpcb*I1_old_)+(0.001*((calc_alpha*I2_old_)-(calc_Kpcf*O_old_))))-((4.0*calc_beta*O_old_)+(calc_gamma*O_old_)));	// 23
	Real d_dt_C2 = (((4.0*calc_alpha*calc_C1)+(2.0*calc_beta*C3_old_))-((calc_beta*C2_old_)+(3.0*calc_alpha*C2_old_)));	// 25
	Real d_dt_C3 = (((3.0*calc_alpha*C2_old_)+(3.0*calc_beta*C4_old_))-((2.0*calc_beta*C3_old_)+(2.0*calc_alpha*C3_old_)));	// 26
	Real d_dt_C4 = (((2.0*calc_alpha*C3_old_)+(4.0*calc_beta*O_old_)+(0.01*((4.0*Kpcb*calc_beta*I1_old_)-(calc_alpha*calc_gamma*C4_old_)))+(0.002*((4.0*calc_beta*I2_old_)-(calc_Kpcf*C4_old_)))+(4.0*calc_beta*Kpcb*I3_old_))-((3.0*calc_beta*C4_old_)+(calc_alpha*C4_old_)+(1.0*calc_gamma*calc_Kpcf*C4_old_)));	// 27
	Real d_dt_I1 = (((calc_gamma*O_old_)+(0.001*((calc_alpha*I3_old_)-(calc_Kpcf*I1_old_)))+(0.01*((calc_alpha*calc_gamma*C4_old_)-(4.0*calc_beta*Kpcb*I1_old_))))-(Kpcb*I1_old_));	// 28
	Real d_dt_I2 = (((0.001*((calc_Kpcf*O_old_)-(calc_alpha*I2_old_)))+(Kpcb*I3_old_)+(0.002*((calc_Kpcf*C4_old_)-(4.0*calc_beta*I2_old_))))-(calc_gamma*I2_old_));	// 29
	Real d_dt_I3 = (((0.001*((calc_Kpcf*I1_old_)-(calc_alpha*I3_old_)))+(calc_gamma*I2_old_)+(1.0*calc_gamma*calc_Kpcf*C4_old_))-((4.0*calc_beta*Kpcb*I3_old_)+(Kpcb*I3_old_)));	// 30
	Real d_dt_Nai = (((-(calc_i_Na+calc_i_Nab+(3.0*calc_i_NaK)+(3.0*calc_i_NaCa)))*Acap*Cm)/(Vmyo*F));	// 39
	Real d_dt_C_Na2 = (((calc_alpha_Na11*calc_C_Na3)+(calc_beta_Na12*C_Na1_old_)+(calc_alpha_Na3*IC_Na2_old_))-((calc_beta_Na11*C_Na2_old_)+(calc_alpha_Na12*C_Na2_old_)+(calc_beta_Na3*C_Na2_old_)));	// 43
	Real d_dt_C_Na1 = (((calc_alpha_Na12*C_Na2_old_)+(calc_beta_Na13*O_Na_old_)+(calc_alpha_Na3*IF_Na_old_))-((calc_beta_Na12*C_Na1_old_)+(calc_alpha_Na13*C_Na1_old_)+(calc_beta_Na3*C_Na1_old_)));	// 44
	Real d_dt_O_Na = (((calc_alpha_Na13*C_Na1_old_)+(calc_beta_Na2*IF_Na_old_))-((calc_beta_Na13*O_Na_old_)+(calc_alpha_Na2*O_Na_old_)));	// 45
	Real d_dt_IF_Na = (((calc_alpha_Na2*O_Na_old_)+(calc_beta_Na3*C_Na1_old_)+(calc_beta_Na4*I1_Na_old_)+(calc_alpha_Na12*IC_Na2_old_))-((calc_beta_Na2*IF_Na_old_)+(calc_alpha_Na3*IF_Na_old_)+(calc_alpha_Na4*IF_Na_old_)+(calc_beta_Na12*IF_Na_old_)));	// 46
	Real d_dt_I1_Na = (((calc_alpha_Na4*IF_Na_old_)+(calc_beta_Na5*I2_Na_old_))-((calc_beta_Na4*I1_Na_old_)+(calc_alpha_Na5*I1_Na_old_)));	// 47
	Real d_dt_I2_Na = ((calc_alpha_Na5*I1_Na_old_)-(calc_beta_Na5*I2_Na_old_));	// 48
	Real d_dt_IC_Na2 = (((calc_alpha_Na11*IC_Na3_old_)+(calc_beta_Na12*IF_Na_old_)+(calc_beta_Na3*C_Na2_old_))-((calc_beta_Na11*IC_Na2_old_)+(calc_alpha_Na12*IC_Na2_old_)+(calc_alpha_Na3*IC_Na2_old_)));	// 49
	Real d_dt_IC_Na3 = (((calc_beta_Na11*IC_Na2_old_)+(calc_beta_Na3*calc_C_Na3))-((calc_alpha_Na11*IC_Na3_old_)+(calc_alpha_Na3*IC_Na3_old_)));	// 50
	Real d_dt_Ki = (((-((calc_i_Kto_f+calc_i_Kto_s+calc_i_K1+calc_i_Ks+calc_i_Kss+calc_i_Kur+calc_i_Kr)-(2.0*calc_i_NaK)))*Acap*Cm)/(Vmyo*F));	// 66
	Real d_dt_ato_f = ((calc_alpha_a*(1.0-ato_f_old_))-(calc_beta_a*ato_f_old_));	// 69
	Real d_dt_ito_f = ((calc_alpha_i*(1.0-ito_f_old_))-(calc_beta_i*ito_f_old_));	// 70
	Real d_dt_ato_s = ((calc_ass-ato_s_old_)/calc_tau_ta_s);	// 76
	Real d_dt_ito_s = ((calc_iss-ito_s_old_)/calc_tau_ti_s);	// 77
	Real d_dt_nKs = ((calc_alpha_n*(1.0-nKs_old_))-(calc_beta_n*nKs_old_));	// 84
	Real d_dt_aur = ((calc_ass-aur_old_)/calc_tau_aur);	// 88
	Real d_dt_iur = ((calc_iss-iur_old_)/calc_tau_iur);	// 89
	Real d_dt_aKss = ((calc_ass-aKss_old_)/calc_tau_Kss);	// 93
	Real d_dt_iKss = 0.0;	// 94
	Real d_dt_C_K2 = (((kf*C_K1_old_)+(calc_beta_a1*O_K_old_))-((kb*C_K2_old_)+(calc_alpha_a1*C_K2_old_)));	// 98
	Real d_dt_C_K1 = (((calc_alpha_a0*calc_C_K0)+(kb*C_K2_old_))-((calc_beta_a0*C_K1_old_)+(kf*C_K1_old_)));	// 99
	Real d_dt_O_K = (((calc_alpha_a1*C_K2_old_)+(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*I_K_old_))-((calc_beta_a1*O_K_old_)+(calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*O_K_old_)));	// 100
	Real d_dt_I_K = ((calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*O_K_old_)-(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*I_K_old_));	// 101

	rDY_[0] = d_dt_V;
	rDY_[1] = d_dt_Cai;
	rDY_[2] = d_dt_Cass;
	rDY_[3] = d_dt_CaJSR;
	rDY_[4] = d_dt_CaNSR;
	rDY_[5] = d_dt_P_RyR;
	rDY_[6] = d_dt_LTRPN_Ca;
	rDY_[7] = d_dt_HTRPN_Ca;
	rDY_[8] = d_dt_P_O1;
	rDY_[9] = d_dt_P_O2;
	rDY_[10] = d_dt_P_C2;
	rDY_[11] = d_dt_O;
	rDY_[12] = d_dt_C2;
	rDY_[13] = d_dt_C3;
	rDY_[14] = d_dt_C4;
	rDY_[15] = d_dt_I1;
	rDY_[16] = d_dt_I2;
	rDY_[17] = d_dt_I3;
	rDY_[18] = d_dt_Nai;
	rDY_[19] = d_dt_C_Na2;
	rDY_[20] = d_dt_C_Na1;
	rDY_[21] = d_dt_O_Na;
	rDY_[22] = d_dt_IF_Na;
	rDY_[23] = d_dt_I1_Na;
	rDY_[24] = d_dt_I2_Na;
	rDY_[25] = d_dt_IC_Na2;
	rDY_[26] = d_dt_IC_Na3;
	rDY_[27] = d_dt_Ki;
	rDY_[28] = d_dt_ato_f;
	rDY_[29] = d_dt_ito_f;
	rDY_[30] = d_dt_ato_s;
	rDY_[31] = d_dt_ito_s;
	rDY_[32] = d_dt_nKs;
	rDY_[33] = d_dt_aur;
	rDY_[34] = d_dt_iur;
	rDY_[35] = d_dt_aKss;
	rDY_[36] = d_dt_iKss;
	rDY_[37] = d_dt_C_K2;
	rDY_[38] = d_dt_C_K1;
	rDY_[39] = d_dt_O_K;
	rDY_[40] = d_dt_I_K;
}
