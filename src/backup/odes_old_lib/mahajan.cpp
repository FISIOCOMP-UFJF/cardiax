#include "mahajan.h"

void setIC_ode_cpu(Real *sv, int cellID, int NEQ = 26)
{
	// Vector containing the initial conditions of the cell
	Real IC[NEQ];

	IC[0] = -87.169816169406;	 // V mV
	IC[1] = 0.990691306716;	 // xh dimensionless
	IC[2] = 0.993888937283;	 // xj dimensionless
	IC[3] = 0.001075453357;	 // xm dimensionless
	IC[4] = 0.000018211252;	 // c1 dimensionless
	IC[5] = 0.979322592773;	 // c2 dimensionless
	IC[6] = 0.001208153482;	 // xi1ca dimensionless
	IC[7] = 0.000033616596;	 // xi1ba dimensionless
	IC[8] = 0.004173008466;	 // xi2ca dimensionless
	IC[9] = 0.015242594688;	 // xi2ba dimensionless
	IC[10] = 0.007074239331;	 // xr dimensionless
	IC[11] = 0.048267587131;	 // xs1 dimensionless
	IC[12] = 0.105468807033;	 // xs2 dimensionless
	IC[13] = 0.00364776906;	 // xtos dimensionless
	IC[14] = 0.174403618112;	 // ytos dimensionless
	IC[15] = 0.003643592594;	 // xtof dimensionless
	IC[16] = 0.993331326442;	 // ytof dimensionless
	IC[17] = 97.505463697266;	 // Ca_JSR uM
	IC[18] = 0.006679257264;	 // xir uM_per_ms
	IC[19] = 11.441712311614;	 // Na_i mM
	IC[20] = 1.716573130685;	 // Ca_dyad uM
	IC[21] = 0.226941113355;	 // Ca_submem uM
	IC[22] = 0.256752008084;	 // Ca_i uM
	IC[23] = 104.450004990523;	 // Ca_NSR uM
	IC[24] = 22.171689894953;	 // tropi uM
	IC[25] = 19.864701949854;	 // trops uM
	for(int i = 0; i < NEQ; i++)
		sv[cellID+i] = IC[i];

}

// Solving the model for each cell in the tissue matrix ni x nj
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

	RHS_mahajan_cpu(time, rY, rDY, col);

	for(int i = 0; i < NEQ; i++)
		sv[cellID+i] = dt*rDY[i] + rY[i];	

}

// SOLVE RHS
void RHS_mahajan_cpu(Real time_, Real *sv_, Real *rDY_, int col)

{

	// State variables
	const Real V_old_ = sv_[0];	 // initial value = -87.169816169406 mV
	const Real xh_old_ = sv_[1];	 // initial value = 0.990691306716 dimensionless
	const Real xj_old_ = sv_[2];	 // initial value = 0.993888937283 dimensionless
	const Real xm_old_ = sv_[3];	 // initial value = 0.001075453357 dimensionless
	const Real c1_old_ = sv_[4];	 // initial value = 0.000018211252 dimensionless
	const Real c2_old_ = sv_[5];	 // initial value = 0.979322592773 dimensionless
	const Real xi1ca_old_ = sv_[6];	 // initial value = 0.001208153482 dimensionless
	const Real xi1ba_old_ = sv_[7];	 // initial value = 0.000033616596 dimensionless
	const Real xi2ca_old_ = sv_[8];	 // initial value = 0.004173008466 dimensionless
	const Real xi2ba_old_ = sv_[9];	 // initial value = 0.015242594688 dimensionless
	const Real xr_old_ = sv_[10];	 // initial value = 0.007074239331 dimensionless
	const Real xs1_old_ = sv_[11];	 // initial value = 0.048267587131 dimensionless
	const Real xs2_old_ = sv_[12];	 // initial value = 0.105468807033 dimensionless
	const Real xtos_old_ = sv_[13];	 // initial value = 0.00364776906 dimensionless
	const Real ytos_old_ = sv_[14];	 // initial value = 0.174403618112 dimensionless
	const Real xtof_old_ = sv_[15];	 // initial value = 0.003643592594 dimensionless
	const Real ytof_old_ = sv_[16];	 // initial value = 0.993331326442 dimensionless
	const Real Ca_JSR_old_ = sv_[17];	 // initial value = 97.505463697266 uM
	const Real xir_old_ = sv_[18];	 // initial value = 0.006679257264 uM_per_ms
	const Real Na_i_old_ = sv_[19];	 // initial value = 11.441712311614 mM
	const Real Ca_dyad_old_ = sv_[20];	 // initial value = 1.716573130685 uM
	const Real Ca_submem_old_ = sv_[21];	 // initial value = 0.226941113355 uM
	const Real Ca_i_old_ = sv_[22];	 // initial value = 0.256752008084 uM
	const Real Ca_NSR_old_ = sv_[23];	 // initial value = 104.450004990523 uM
	const Real tropi_old_ = sv_[24];	 // initial value = 22.171689894953 uM
	const Real trops_old_ = sv_[25];	 // initial value = 19.864701949854 uM


	// Parameters
	const Real F = 96.4853415;	 // coulomb_per_mmole
	const Real R = 8.314472;	 // J_per_moleK
	const Real T = 308;	 // kelvin
	const Real stim_period = 400;	 // ms
	const Real stim_amplitude = 0.0;	 // nA_per_nF
	const Real stim_offset = 0;	 // ms
	const Real stim_duration = 3;	 // ms
	const Real gna = 12;	 // uS_per_nF
	const Real pca = 0.00054;	 // cm_per_s
	const Real Ca_o = 1.8;	 // mM
	const Real vth = 0;	 // mV
	const Real s6 = 8;	 // mV
	const Real taupo = 1;	 // ms
	const Real cat = 3;	 // uM
	const Real r1 = 0.3;	 // per_ms
	const Real r2 = 3;	 // per_ms
	const Real k2 = 1.03615e-4;	 // per_ms
	const Real s1t = 0.00195;	 // per_ms
	const Real k2t = 0.00224;	 // per_ms
	const Real k1t = 0.00413;	 // per_ms
	const Real vx = -40;	 // mV
	const Real sx = 3;	 // mV
	const Real tau3 = 3;	 // ms
	const Real vy = -40;	 // mV
	const Real sy = 4;	 // mV
	const Real tca = 78.0329;	 // ms
	const Real cpt = 6.09365;	 // uM
	const Real vyr = -40;	 // mV
	const Real syr = 11.32;	 // mV
	const Real gca = 182;	 // mmole_per_coulomb_cm
	const Real wca = 8;	 // mV_per_uM
	const Real gkix = 0.3;	 // uS_per_nF
	const Real K_o = 5.4;	 // mM
	const Real gkr = 0.0125;	 // uS_per_nF
	const Real gks = 0.1386;	 // uS_per_nF
	const Real gtos = 0.04;	 // uS_per_nF
	const Real gtof = 0.11;	 // uS_per_nF
	const Real Na_o = 136;	 // mM
	const Real gNaK = 1.5;	 // nA_per_nF
	const Real xkmnai = 12;	 // mM
	const Real xkmko = 1.5;	 // mM
	const Real xkdna = 0.3;	 // uM
	const Real xmcao = 1.3;	 // mM
	const Real xmnao = 87.5;	 // mM
	const Real xmnai = 12.3;	 // mM
	const Real xmcai = 0.0036;	 // mM
	const Real gNaCa = 0.84;	 // uM_per_ms
	const Real taua = 100;	 // ms
	const Real av = 11.3;	 // per_ms
	const Real cstar = 90;	 // uM
	const Real ay = 0.05;	 // per_mV
	const Real gryr = 2.58079;	 // per_ms
	const Real taur = 30;	 // ms
	const Real gbarsr = 26841.8;	 // dimensionless
	const Real ax = 0.3576;	 // per_mV
	const Real gdyad = 9000;	 // mmole_per_coulomb_cm
	const Real vup = 0.4;	 // uM_per_ms
	const Real cup = 0.5;	 // uM
	const Real gleak = 0.00002069;	 // per_ms
	const Real kj = 50;	 // uM
	const Real bcal = 24;	 // uM
	const Real xkcal = 7;	 // uM
	const Real srmax = 47;	 // uM
	const Real srkd = 0.6;	 // uM
	const Real bmem = 15;	 // uM
	const Real kmem = 0.3;	 // uM
	const Real bsar = 42;	 // uM
	const Real ksar = 13;	 // uM
	const Real taud = 4;	 // ms
	const Real xkon = 0.0327;	 // per_uM_per_ms
	const Real btrop = 70;	 // uM
	const Real xkoff = 0.0196;	 // per_ms
	const Real taups = 0.5;	 // ms
	const Real K_i = 140;	 // mM
	const Real prNaK = 0.01833;	 // dimensionless

	// Independent Variable
	Real time_new = time_;


	// Algebraic Equations
	Real calc_FonRT = (F/(R*T));	//0
	Real calc_past = (floor((time_new/stim_period))*stim_period);	//1
	Real calc_am = ((fabs((V_old_+47.13))>0.001))? ((0.32*1.0*(V_old_+47.13))/(1.0-exp(((-0.1)*(V_old_+47.13))))): 3.2;	//5
	Real calc_bm = (0.08*exp(((-V_old_)/11.0)));	//6
	Real calc_ah = ((V_old_<(-40.0)))? (0.135*exp(((80.0+V_old_)/(-6.8)))): 0.0;	//7
	Real calc_bh = ((V_old_<(-40.0)))? ((3.56*exp((0.079*V_old_)))+(310000.0*exp((0.35*V_old_)))): (1.0/(0.13*(1.0+exp(((V_old_+10.66)/(-11.1))))));	//8
	Real calc_aj = ((V_old_<(-40.0)))? (((((-127140.0)*exp((0.2444*V_old_)))-(0.00003474*exp(((-0.04391)*V_old_))))*1.0*(V_old_+37.78))/(1.0+exp((0.311*(V_old_+79.23))))): 0.0;	//9
	Real calc_bj = ((V_old_<(-40.0)))? ((0.1212*exp(((-0.01052)*V_old_)))/(1.0+exp(((-0.1378)*(V_old_+40.14))))): ((0.3*exp(((-0.0000002535)*V_old_)))/(1.0+exp(((-0.1)*(V_old_+32.0)))));	//10
	Real calc_poinf = (1.0/(1.0+exp(((-(V_old_-vth))/s6))));	//17
	Real calc_fca = (1.0/(1.0+pow((cat/Ca_dyad_old_),3.0)));	//20
	Real calc_s2t = ((((s1t*r1)/r2)*k2t)/k1t);	//24
	Real calc_poi = (1.0/(1.0+exp(((-(V_old_-vx))/sx))));	//25
	Real calc_Pr = (1.0-(1.0/(1.0+exp(((-(V_old_-vy))/sy)))));	//28
	Real calc_recov = (10.0+(4954.0*exp((V_old_/15.6))));	//29
	Real calc_tau_ca = ((tca/(1.0+pow((Ca_dyad_old_/cpt),4.0)))+0.1);	//30
	Real calc_Ps = (1.0/(1.0+exp(((-(V_old_-vyr))/syr))));	//33
	Real calc_po = ((((((1.0-xi1ca_old_)-xi2ca_old_)-xi1ba_old_)-xi2ba_old_)-c1_old_)-c2_old_);	//40
	Real calc_xkrv1 = ((fabs((V_old_+7.0))>0.001))? ((0.00138*1.0*(V_old_+7.0))/(1.0-exp(((-0.123)*(V_old_+7.0))))): (0.00138/0.123);	//53
	Real calc_xkrv2 = ((fabs((V_old_+10.0))>0.001))? ((0.00061*1.0*(V_old_+10.0))/(exp((0.145*(V_old_+10.0)))-1.0)): (0.00061/0.145);	//54
	Real calc_xkrinf = (1.0/(1.0+exp(((-(V_old_+50.0))/7.5))));	//56
	Real calc_rg = (1.0/(1.0+exp(((V_old_+33.0)/22.4))));	//57
	Real calc_xs1ss = (1.0/(1.0+exp(((-(V_old_-1.5))/16.7))));	//60
	Real calc_tauxs1 = ((fabs((V_old_+30.0))<(0.001/0.0687)))? (1.0/((0.0000719/0.148)+(0.000131/0.0687))): (1.0/(((0.0000719*(V_old_+30.0))/(1.0-exp(((-0.148)*(V_old_+30.0)))))+((0.000131*(V_old_+30.0))/(exp((0.0687*(V_old_+30.0)))-1.0))));	//62
	Real calc_gksx = (1.0+(0.8/(1.0+pow((0.5/Ca_i_old_),3.0))));	//64
	Real calc_rt1 = ((-(V_old_+3.0))/15.0);	//68
	Real calc_rt2 = ((V_old_+33.5)/10.0);	//69
	Real calc_rt3 = ((V_old_+60.0)/10.0);	//70
	Real calc_rt4 = ((((-V_old_)/30.0)*V_old_)/30.0);	//71
	Real calc_rt5 = ((V_old_+33.5)/10.0);	//72
	Real calc_sigma = ((exp((Na_o/67.3))-1.0)/7.0);	//89
	Real calc_aloss = (1.0/(1.0+pow((xkdna/Ca_submem_old_),3.0)));	//94
	Real calc_yz3 = (xmcai*pow(Na_o,3.0)*(1.0+pow((Na_i_old_/xmnai),3.0)));	//97
	Real calc_bv = (((1.0-av)*cstar)-50.0);	//103
	Real calc_Qr0 = (((Ca_JSR_old_>50.0)&&(Ca_JSR_old_<cstar)))? ((Ca_JSR_old_-50.0)/1.0): (((Ca_JSR_old_>=cstar))? ((av*Ca_JSR_old_)+calc_bv): 0.0);	//104
	Real calc_sparkV = (exp(((-ay)*(V_old_+30.0)))/(1.0+exp(((-ay)*(V_old_+30.0)))));	//106
	Real calc_jup = ((vup*Ca_i_old_*Ca_i_old_)/((Ca_i_old_*Ca_i_old_)+(cup*cup)));	//112
	Real calc_jleak = (((gleak*Ca_NSR_old_*Ca_NSR_old_)/((Ca_NSR_old_*Ca_NSR_old_)+(kj*kj)))*((Ca_NSR_old_*16.667)-Ca_i_old_));	//113
	Real calc_bpxs = ((bcal*xkcal)/((xkcal+Ca_submem_old_)*(xkcal+Ca_submem_old_)));	//115
	Real calc_spxs = ((srmax*srkd)/((srkd+Ca_submem_old_)*(srkd+Ca_submem_old_)));	//116
	Real calc_mempxs = ((bmem*kmem)/((kmem+Ca_submem_old_)*(kmem+Ca_submem_old_)));	//117
	Real calc_sarpxs = ((bsar*ksar)/((ksar+Ca_submem_old_)*(ksar+Ca_submem_old_)));	//118
	Real calc_bpxi = ((bcal*xkcal)/((xkcal+Ca_i_old_)*(xkcal+Ca_i_old_)));	//120
	Real calc_spxi = ((srmax*srkd)/((srkd+Ca_i_old_)*(srkd+Ca_i_old_)));	//121
	Real calc_mempxi = ((bmem*kmem)/((kmem+Ca_i_old_)*(kmem+Ca_i_old_)));	//122
	Real calc_sarpxi = ((bsar*ksar)/((ksar+Ca_i_old_)*(ksar+Ca_i_old_)));	//123
	Real calc_jd = ((Ca_submem_old_-Ca_i_old_)/taud);	//125
	Real calc_xbi = ((xkon*Ca_i_old_*(btrop-tropi_old_))-(xkoff*tropi_old_));	//126
	Real calc_xbs = ((xkon*Ca_submem_old_*(btrop-trops_old_))-(xkoff*trops_old_));	//127
	Real calc_csm = (Ca_submem_old_/1000.0);	//135
	Real calc_i_Stim = ((((time_new-calc_past)>=stim_offset)&&((time_new-calc_past)<=(stim_offset+stim_duration))))? stim_amplitude: 0.0;	//2
	Real calc_za = (V_old_*2.0*calc_FonRT);	//15
	Real calc_alpha = (calc_poinf/taupo);	//18
	Real calc_beta = ((1.0-calc_poinf)/taupo);	//19
	Real calc_s1 = (0.0182688*calc_fca);	//21
	Real calc_k1 = (0.024168*calc_fca);	//22
	Real calc_k3 = ((1.0-calc_poi)/tau3);	//26
	Real calc_tauca = (((calc_recov-calc_tau_ca)*calc_Pr)+calc_tau_ca);	//31
	Real calc_tauba = (((calc_recov-450.0)*calc_Pr)+450.0);	//32
	Real calc_taukr = (1.0/(calc_xkrv1+calc_xkrv2));	//55
	Real calc_xs2ss = calc_xs1ss;	//61
	Real calc_tauxs2 = (4.0*calc_tauxs1);	//63
	Real calc_xtos_inf = (1.0/(1.0+exp(calc_rt1)));	//73
	Real calc_ytos_inf = (1.0/(1.0+exp(calc_rt2)));	//74
	Real calc_rs_inf = (1.0/(1.0+exp(calc_rt2)));	//77
	Real calc_txs = ((9.0/(1.0+exp((-calc_rt1))))+0.5);	//78
	Real calc_tys = ((3000.0/(1.0+exp(calc_rt3)))+30.0);	//79
	Real calc_txf = ((3.5*exp(calc_rt4))+1.5);	//80
	Real calc_tyf = ((20.0/(1.0+exp(calc_rt5)))+20.0);	//81
	Real calc_fNaK = (1.0/(1.0+(0.1245*exp(((-0.1)*V_old_*calc_FonRT)))+(0.0365*calc_sigma*exp(((-V_old_)*calc_FonRT)))));	//90
	Real calc_zw4 = (1.0+(0.2*exp((V_old_*(0.35-1.0)*calc_FonRT))));	//93
	Real calc_Qr = ((Ca_NSR_old_*calc_Qr0)/cstar);	//105
	Real calc_dcsib = (1.0/(1.0+calc_bpxs+calc_spxs+calc_mempxs+calc_sarpxs));	//119
	Real calc_dciib = (1.0/(1.0+calc_bpxi+calc_spxi+calc_mempxi+calc_sarpxi));	//124
	Real calc_dCa_JSR = (((-xir_old_)+calc_jup)-calc_jleak);	//128
	Real calc_ek = ((1.0/calc_FonRT)*log((K_o/K_i)));	//136
	Real calc_eks = ((1.0/calc_FonRT)*log(((K_o+(prNaK*Na_o))/(K_i+(prNaK*Na_i_old_)))));	//137
	Real calc_ena = ((1.0/calc_FonRT)*log((Na_o/Na_i_old_)));	//138
	Real calc_rxa = ((fabs(calc_za)<0.001))? ((4.0*pca*F*calc_FonRT*((calc_csm*exp(calc_za))-(0.341*Ca_o)))/(2.0*calc_FonRT)): ((4.0*pca*V_old_*F*calc_FonRT*((calc_csm*exp(calc_za))-(0.341*Ca_o)))/(exp(calc_za)-1.0));	//16
	Real calc_s2 = ((((calc_s1*r1)/r2)*k2)/calc_k1);	//23
	Real calc_k3t = calc_k3;	//27
	Real calc_k6 = ((calc_fca*calc_Ps)/calc_tauca);	//34
	Real calc_k5 = ((1.0-calc_Ps)/calc_tauca);	//35
	Real calc_k6t = (calc_Ps/calc_tauba);	//36
	Real calc_k5t = ((1.0-calc_Ps)/calc_tauba);	//37
	Real calc_xtof_inf = calc_xtos_inf;	//75
	Real calc_ytof_inf = calc_ytos_inf;	//76
	Real calc_xiNaK = ((((gNaK*calc_fNaK*Na_i_old_)/(Na_i_old_+xkmnai))*K_o)/(K_o+xkmko));	//91
	Real calc_zw3 = ((pow(Na_i_old_,3.0)*Ca_o*exp((V_old_*0.35*calc_FonRT)))-(pow(Na_o,3.0)*calc_csm*exp((V_old_*(0.35-1.0)*calc_FonRT))));	//92
	Real calc_yz1 = ((xmcao*pow(Na_i_old_,3.0))+(pow(xmnao,3.0)*calc_csm));	//95
	Real calc_yz2 = (pow(xmnai,3.0)*Ca_o*(1.0+(calc_csm/xmcai)));	//96
	Real calc_yz4 = ((pow(Na_i_old_,3.0)*Ca_o)+(pow(Na_o,3.0)*calc_csm));	//98
	Real calc_xina = (gna*xh_old_*xj_old_*xm_old_*xm_old_*xm_old_*(V_old_-calc_ena));	//14
	Real calc_k4 = ((((((calc_k3*calc_alpha)/calc_beta)*calc_k1)/k2)*calc_k5)/calc_k6);	//38
	Real calc_k4t = ((((((calc_k3t*calc_alpha)/calc_beta)*k1t)/k2t)*calc_k5t)/calc_k6t);	//39
	Real calc_jca = (gca*calc_po*calc_rxa);	//47
	Real calc_aki = (1.02/(1.0+exp((0.2385*((V_old_-calc_ek)-59.215)))));	//49
	Real calc_bki = (((0.49124*exp((0.08032*((V_old_-calc_ek)+5.476))))+(1.0*exp((0.06175*((V_old_-calc_ek)-594.31)))))/(1.0+exp(((-0.5143)*((V_old_-calc_ek)+4.753)))));	//50
	Real calc_xikr = (gkr*pow((K_o/5.4),1.0/2.0)*xr_old_*calc_rg*(V_old_-calc_ek));	//58
	Real calc_xiks = (gks*calc_gksx*xs1_old_*xs2_old_*(V_old_-calc_eks));	//65
	Real calc_xitos = (gtos*xtos_old_*(ytos_old_+(0.5*calc_rs_inf))*(V_old_-calc_ek));	//86
	Real calc_xitof = (gtof*xtof_old_*ytof_old_*(V_old_-calc_ek));	//87
	Real calc_zw8 = (calc_yz1+calc_yz2+calc_yz3+calc_yz4);	//99
	Real calc_spark_rate = ((gryr/1.0)*calc_po*fabs(calc_rxa)*calc_sparkV);	//107
	Real calc_xirp = ((((calc_po*calc_Qr*fabs(calc_rxa)*gbarsr)/1.0)*exp(((-ax)*(V_old_+30.0))))/(1.0+exp(((-ax)*(V_old_+30.0)))));	//109
	Real calc_xicap = (calc_po*gdyad*fabs(calc_rxa));	//110
	Real calc_xica = (2.0*wca*calc_jca);	//48
	Real calc_xkin = (calc_aki/(calc_aki+calc_bki));	//51
	Real calc_xito = (calc_xitos+calc_xitof);	//88
	Real calc_jNaCa = ((gNaCa*calc_aloss*calc_zw3)/(calc_zw4*calc_zw8));	//100
	Real calc_xiryr = (calc_xirp+calc_xicap);	//111
	Real calc_xik1 = (gkix*pow((K_o/5.4),1.0/2.0)*calc_xkin*(V_old_-calc_ek));	//52
	Real calc_xiNaCa = (wca*calc_jNaCa);	//101
	Real calc_Itotal = (-(calc_xina+calc_xik1+calc_xikr+calc_xiks+calc_xito+calc_xiNaCa+calc_xica+calc_xiNaK+calc_i_Stim));	//3

	// Differential Equations
	Real d_dt_V = calc_Itotal;	// 4
	Real d_dt_xh = ((calc_ah*(1.0-xh_old_))-(calc_bh*xh_old_));	// 11
	Real d_dt_xj = ((calc_aj*(1.0-xj_old_))-(calc_bj*xj_old_));	// 12
	Real d_dt_xm = ((calc_am*(1.0-xm_old_))-(calc_bm*xm_old_));	// 13
	Real d_dt_c1 = (((calc_alpha*c2_old_)+(k2*xi1ca_old_)+(k2t*xi1ba_old_)+(r2*calc_po))-((calc_beta+r1+k1t+calc_k1)*c1_old_));	// 41
	Real d_dt_c2 = (((calc_beta*c1_old_)+(calc_k5*xi2ca_old_)+(calc_k5t*xi2ba_old_))-((calc_k6+calc_k6t+calc_alpha)*c2_old_));	// 42
	Real d_dt_xi1ca = (((calc_k1*c1_old_)+(calc_k4*xi2ca_old_)+(calc_s1*calc_po))-((calc_k3+k2+calc_s2)*xi1ca_old_));	// 43
	Real d_dt_xi1ba = (((k1t*c1_old_)+(calc_k4t*xi2ba_old_)+(s1t*calc_po))-((calc_k3t+k2t+calc_s2t)*xi1ba_old_));	// 44
	Real d_dt_xi2ca = (((calc_k3*xi1ca_old_)+(calc_k6*c2_old_))-((calc_k5+calc_k4)*xi2ca_old_));	// 45
	Real d_dt_xi2ba = (((calc_k3t*xi1ba_old_)+(calc_k6t*c2_old_))-((calc_k5t+calc_k4t)*xi2ba_old_));	// 46
	Real d_dt_xr = ((calc_xkrinf-xr_old_)/calc_taukr);	// 59
	Real d_dt_xs1 = ((calc_xs1ss-xs1_old_)/calc_tauxs1);	// 66
	Real d_dt_xs2 = ((calc_xs2ss-xs2_old_)/calc_tauxs2);	// 67
	Real d_dt_xtos = ((calc_xtos_inf-xtos_old_)/calc_txs);	// 82
	Real d_dt_ytos = ((calc_ytos_inf-ytos_old_)/calc_tys);	// 83
	Real d_dt_xtof = ((calc_xtof_inf-xtof_old_)/calc_txf);	// 84
	Real d_dt_ytof = ((calc_ytof_inf-ytof_old_)/calc_tyf);	// 85
	Real d_dt_Ca_JSR = ((Ca_NSR_old_-Ca_JSR_old_)/taua);	// 102
	Real d_dt_xir = ((calc_spark_rate*calc_Qr)-((xir_old_*(1.0-((taur*calc_dCa_JSR)/Ca_NSR_old_)))/taur));	// 108
	Real d_dt_Na_i = ((-(calc_xina+(3.0*calc_xiNaK)+(3.0*calc_xiNaCa)))/(wca*1000.0));	// 114
	Real d_dt_Ca_dyad = (calc_xiryr-((Ca_dyad_old_-Ca_submem_old_)/taups));	// 129
	Real d_dt_Ca_submem = (calc_dcsib*((50.0*(((xir_old_-calc_jd)-calc_jca)+calc_jNaCa))-calc_xbs));	// 130
	Real d_dt_Ca_i = (calc_dciib*(((calc_jd-calc_jup)+calc_jleak)-calc_xbi));	// 131
	Real d_dt_Ca_NSR = calc_dCa_JSR;	// 132
	Real d_dt_tropi = calc_xbi;	// 133
	Real d_dt_trops = calc_xbs;	// 134

	rDY_[0] = d_dt_V;
	rDY_[1] = d_dt_xh;
	rDY_[2] = d_dt_xj;
	rDY_[3] = d_dt_xm;
	rDY_[4] = d_dt_c1;
	rDY_[5] = d_dt_c2;
	rDY_[6] = d_dt_xi1ca;
	rDY_[7] = d_dt_xi1ba;
	rDY_[8] = d_dt_xi2ca;
	rDY_[9] = d_dt_xi2ba;
	rDY_[10] = d_dt_xr;
	rDY_[11] = d_dt_xs1;
	rDY_[12] = d_dt_xs2;
	rDY_[13] = d_dt_xtos;
	rDY_[14] = d_dt_ytos;
	rDY_[15] = d_dt_xtof;
	rDY_[16] = d_dt_ytof;
	rDY_[17] = d_dt_Ca_JSR;
	rDY_[18] = d_dt_xir;
	rDY_[19] = d_dt_Na_i;
	rDY_[20] = d_dt_Ca_dyad;
	rDY_[21] = d_dt_Ca_submem;
	rDY_[22] = d_dt_Ca_i;
	rDY_[23] = d_dt_Ca_NSR;
	rDY_[24] = d_dt_tropi;
	rDY_[25] = d_dt_trops;
}

