#include "torord_land.hpp"

TorordLand::TorordLand() : CellModel(50) {
  var_names.insert( std::pair<int, std::string>(0, "v") );
  var_names.insert( std::pair<int, std::string>(1, "nai") );
  var_names.insert( std::pair<int, std::string>(2, "nass") );
  var_names.insert( std::pair<int, std::string>(3, "ki") );
  var_names.insert( std::pair<int, std::string>(4, "kss") );
  var_names.insert( std::pair<int, std::string>(5, "cai") );
  var_names.insert( std::pair<int, std::string>(6, "cass") );
  var_names.insert( std::pair<int, std::string>(7, "cansr") );
  var_names.insert( std::pair<int, std::string>(8, "cajsr") );
  var_names.insert( std::pair<int, std::string>(9, "m") );
  var_names.insert( std::pair<int, std::string>(10, "hp") );
  var_names.insert( std::pair<int, std::string>(11, "h") );
  var_names.insert( std::pair<int, std::string>(12, "j") );
  var_names.insert( std::pair<int, std::string>(13, "jp") );
  var_names.insert( std::pair<int, std::string>(14, "mL") );
  var_names.insert( std::pair<int, std::string>(15, "hL") );
  var_names.insert( std::pair<int, std::string>(16, "hLp") );
  var_names.insert( std::pair<int, std::string>(17, "a") );
  var_names.insert( std::pair<int, std::string>(18, "iF") );
  var_names.insert( std::pair<int, std::string>(19, "iS") );
  var_names.insert( std::pair<int, std::string>(20, "ap") );
  var_names.insert( std::pair<int, std::string>(21, "iFp") );
  var_names.insert( std::pair<int, std::string>(22, "iSp") );
  var_names.insert( std::pair<int, std::string>(23, "d") );
  var_names.insert( std::pair<int, std::string>(24, "ff") );
  var_names.insert( std::pair<int, std::string>(25, "fs") );
  var_names.insert( std::pair<int, std::string>(26, "fcaf") );
  var_names.insert( std::pair<int, std::string>(27, "fcas") );
  var_names.insert( std::pair<int, std::string>(28, "jca") );
  var_names.insert( std::pair<int, std::string>(29, "nca") );
  var_names.insert( std::pair<int, std::string>(30, "nca_i") );
  var_names.insert( std::pair<int, std::string>(31, "ffp") );
  var_names.insert( std::pair<int, std::string>(32, "fcafp") );
  var_names.insert( std::pair<int, std::string>(33, "xs1") );
  var_names.insert( std::pair<int, std::string>(34, "xs2") );
  var_names.insert( std::pair<int, std::string>(35, "Jrelnp") );
  var_names.insert( std::pair<int, std::string>(36, "CaMKt") );
  var_names.insert( std::pair<int, std::string>(37, "c0") );
  var_names.insert( std::pair<int, std::string>(38, "c1") );
  var_names.insert( std::pair<int, std::string>(39, "c2") );
  var_names.insert( std::pair<int, std::string>(40, "elta_o") );
  var_names.insert( std::pair<int, std::string>(41, "i") );
  var_names.insert( std::pair<int, std::string>(42, "Jrelp") );
  var_names.insert( std::pair<int, std::string>(43, "XS") );
  var_names.insert( std::pair<int, std::string>(44, "XW") );
  var_names.insert( std::pair<int, std::string>(45, "Ca_TRPN") );
  var_names.insert( std::pair<int, std::string>(46, "TmBlocked") );
  var_names.insert( std::pair<int, std::string>(47, "ZETAS") );
  var_names.insert( std::pair<int, std::string>(48, "ZETAW") );
  var_names.insert( std::pair<int, std::string>(49, "Ta") );

  //monitored.push_back( &active );
  rlvars.insert(49);

}

void TorordLand::init(double * values) const
{
  assert(values != nullptr);
            // Default initial conditions for ENDO cell (from original Matlab script)
            values[0] = -8.863699e+01;
            values[1] = 1.189734e+01;
            values[2] = 1.189766e+01;
            values[3] = 1.412345e+02;
            values[4] = 1.412344e+02;
            values[5] = 7.267473e-05;
            values[6] = 6.337870e-05;
            values[7] = 1.532653e+00;
            values[8] = 1.533946e+00;
            values[9] = 8.280078e-04;
            values[10] = 6.665272e-01;
            values[11] = 8.260208e-01;
            values[12] = 8.260560e-01;
            values[13] = 8.258509e-01;
            values[14] = 1.668686e-04;
            values[15] = 5.228306e-01;
            values[16] = 2.859696e-01;
            values[17] = 9.591370e-04;
            values[18] = 9.996012e-01;
            values[19] = 5.934016e-01;
            values[20] = 4.886961e-04;
            values[21] = 9.996011e-01;
            values[22] = 6.546687e-01;
            values[23] = 9.500075e-32;
            values[24] = 1.000000e+00;
            values[25] = 9.392580e-01;
            values[26] = 1.000000e+00;
            values[27] = 9.998984e-01;
            values[28] = 9.999783e-01;
            values[29] = 4.448162e-04;
            values[30] = 7.550725e-04;
            values[31] = 1.000000e+00;
            values[32] = 1.000000e+00;
            values[33] = 2.424047e-01;
            values[34] = 1.795377e-04;
            values[35] = -6.883086e-25;
            values[36] = 1.117498e-02;
            values[37] = 9.980366e-01;
            values[38] = 8.588018e-04;
            values[39] = 7.097447e-04;
            values[40] = 3.812617e-04;
            values[41] = 1.357116e-05;
            values[42] = 2.302525e-23;
            values[43] = 1.561941e-04;
            values[44] = 2.351289e-04;
            values[45] = 8.077631e-03;
            values[46] = 9.993734e-01;
            values[47] = 0.000000e+00;
            values[48] = 0.000000e+00;
	          values[49] = 0.000000e+00;

            // Default initial conditions for MID cell (from original Matlab script)
            //values[0] = -8.953800e+01;
            //values[1] = 1.492920e+01;
            //values[2] = 1.492967e+01;
            //values[3] = 1.448447e+02;
            //values[4] = 1.448447e+02;
            //values[5] = 7.502288e-05;
            //values[6] = 6.107636e-05;
            //values[7] = 1.790435e+00;
            //values[8] = 1.794842e+00;
            //values[9] = 6.819365e-04;
            //values[10] = 6.953807e-01;
            //values[11] = 8.434888e-01;
            //values[12] = 8.435208e-01;
            //values[13] = 8.432262e-01;
            //values[14] = 1.406211e-04;
            //values[15] = 5.453149e-01;
            //values[16] = 2.924967e-01;
            //values[17] = 9.026127e-04;
            //values[18] = 9.996593e-01;
            //values[19] = 5.631197e-01;
            //values[20] = 4.598833e-04;
            //values[21] = 9.996593e-01;
            //values[22] = 6.236964e-01;
            //values[23] = -1.314189e-33;
            //values[24] = 1.000000e+00;
            //values[25] = 9.204086e-01;
            //values[26] = 1.000000e+00;
            //values[27] = 9.997620e-01;
            //values[28] = 9.999625e-01;
            //values[29] = 3.853595e-04;
            //values[30] = 8.535292e-04;
            //values[31] = 1.000000e+00;
            //values[32] = 1.000000e+00;
            //values[33] = 2.664151e-01;
            //values[34] = 1.623107e-04;
            //values[35] = 1.209762e-24;
            //values[36] = 1.782437e-02;
            //values[37] = 9.979720e-01;
            //values[38] = 8.053991e-04;
            //values[39] = 6.781800e-04;
            //values[40] = 5.265363e-04;
            //values[41] = 1.789565e-05;
            //values[42] = 7.059162e-23;
            //values[43] = 1.670654e-04;
            //values[44] = 2.506794e-04;
            //values[45] = 8.602625e-03;
            //values[46] = 9.993314e-01;
            //values[47] = 0.000000e+00;
            //values[48] = 0.000000e+00;

            // Default initial conditions for EPI cell (from original Matlab script)
            //values[0] = -8.904628e+01;
            //values[1] = 1.272190e+01;
            //values[2] = 1.272220e+01;
            //values[3] = 1.422490e+02;
            //values[4] = 1.422489e+02;
            //values[5] = 6.541058e-05;
            //values[6] = 5.684431e-05;
            //values[7] = 1.809117e+00;
            //values[8] = 1.809702e+00;
            //values[9] = 7.581821e-04;
            //values[10] = 6.798398e-01;
            //values[11] = 8.341502e-01;
            //values[12] = 8.341883e-01;
            //values[13] = 8.340817e-01;
            //values[14] = 1.543877e-04;
            //values[15] = 5.382951e-01;
            //values[16] = 3.027694e-01;
            //values[17] = 9.330351e-04;
            //values[18] = 9.996287e-01;
            //values[19] = 9.996262e-01;
            //values[20] = 4.753907e-04;
            //values[21] = 9.996287e-01;
            //values[22] = 9.996285e-01;
            //values[23] = 1.742134e-37;
            //values[24] = 1.000000e+00;
            //values[25] = 9.479522e-01;
            //values[26] = 1.000000e+00;
            //values[27] = 9.999327e-01;
            //values[28] = 9.999829e-01;
            //values[29] = 2.915447e-04;
            //values[30] = 5.026045e-04;
            //values[31] = 1.000000e+00;
            //values[32] = 1.000000e+00;
            //values[33] = 2.288155e-01;
            //values[34] = 1.714978e-04;
            //values[35] = -1.131190e-26;
            //values[36] = 1.295052e-02;
            //values[37] = 9.981944e-01;
            //values[38] = 8.342321e-04;
            //values[39] = 6.838658e-04;
            //values[40] = 2.778785e-04;
            //values[41] = 9.667759e-06;
            //values[42] = 8.169304e-24;
            //values[43] = 1.259996e-04;
            //values[44] = 1.899522e-04;
            //values[45] = 6.551494e-03;
            //values[46] = 9.994940e-01;
            //values[47] = 0.000000e+00;
            //values[48] = 0.000000e+00;
}

void TorordLand::equation(const double time, const double *rY, double *rDY)
{
    // Current modifiers
    double INa_Multiplier =  1.0;
    double INaL_Multiplier = 1.0;
    double INaCa_Multiplier = 1.0;
    double INaK_Multiplier = 1.0;
    double INab_Multiplier =  1.0;
    double Ito_Multiplier = 1.0;
    double IKr_Multiplier =  1.0;
    double IKs_Multiplier =  1.0;
    double IK1_Multiplier = 1.0;
    double IKb_Multiplier = 1.0;
    double IKCa_Multiplier =  0.0;
    double ICaL_Multiplier = 1.0;
    double ICab_Multiplier = 1.0;
    double IpCa_Multiplier = 1.0;
    double ICaCl_Multiplier =  1.0;
    double IClb_Multiplier = 1.0;
    double Jrel_Multiplier = 1.0;
    double Jup_Multiplier = 1.0;
    double aCaMK_Multiplier = 1.0;
    double taurelp_Multiplier = 1.0;

    // Get the celltype for the current cell
    double celltype = type;

    // Get the stimulus current from the current cell
    double calc_I_stim = i_stim;

    // State variables (same order as the original Matlab script)
    double v = rY[0];
    double nai = rY[1];
    double nass = rY[2];
    double ki = rY[3];
    double kss = rY[4];
    double cai = rY[5];
    double cass = rY[6];
    double cansr = rY[7];
    double cajsr = rY[8];
    double m = rY[9];
    double hp = rY[10];
    double h = rY[11];
    double j = rY[12];
    double jp = rY[13];
    double mL = rY[14];
    double hL = rY[15];
    double hLp = rY[16];
    double a = rY[17];
    double iF = rY[18];
    double iS = rY[19];
    double ap = rY[20];
    double iFp = rY[21];
    double iSp = rY[22];

    // ical
    double d = rY[23];
    double ff = rY[24];
    double fs = rY[25];
    double fcaf = rY[26];
    double fcas = rY[27];
    double jca = rY[28];
    double nca = rY[29];
    double nca_i = rY[30];
    double ffp = rY[31];
    double fcafp = rY[32];

    double xs1 = rY[33];
    double xs2 = rY[34];
    double Jrel_np = rY[35];
    double CaMKt = rY[36];

    // new MM ICaL states
    double ikr_c0 = rY[37];
    double ikr_c1 = rY[38];
    double ikr_c2 = rY[39];
    double ikr_o = rY[40];
    double ikr_i = rY[41];
    double Jrel_p = rY[42];

    const double cli = 24;   // Intracellular Cl  [mM]
    const double clo = 150;  // Extracellular Cl  [mM]
// -----------------------------------------------------
    // Land-Niederer model
    double XS = fmaxf(0,rY[43]);
    double XW = fmaxf(0,rY[44]);
    double Ca_TRPN = fmaxf(0,rY[45]);
    double TmBlocked = rY[46];
    double ZETAS = rY[47];
    double ZETAW = rY[48];
    double TA = rY[49];

// Avoid NaN errors
//if (v >= 100.0) v = 100.0;
//if (v <= -100.0) v = -100.0;

//double aux, aux2, aux3;

// Changeable parameters
double nao = 140.0;
double cao = 1.8;
double ko = 5.0;
// Localization of ICaL and NCX: the fraction in junctional subspace
double ICaL_fractionSS = 0.8; 
double INaCa_fractionSS = 0.35;
// Additional parameters
double ca50 = 0.805;
double trpnmax = 0.07;

// INPUT CODE:
int mode = 0;           // 0 = "intact", 1 = "skinned"
double lambda = 1.0;
double lambda_rate = 0.0;

// EC parameters
double perm50 = 0.35;
double TRPN_n = 2;
double koff = 0.1;
double dr = 0.25;
double wfrac = 0.5;
double TOT_A = 25;
double ktm_unblock = 0.021; 
double beta_1 = -2.4;
double beta_0 = 2.3;
double gamma = 0.0085;
double gamma_wu = 0.615;
double phi = 2.23;

double nperm;
double Tref;
double nu;
double mu;
if (mode == 1) {
    nperm = 2.2;
    //ca50=2.5;
    Tref = 40.5;
    nu = 1;
    mu = 1;
}
else {
    nperm = 2.036; 
    //ca50=0.805;
    //ca50 = 0.5;
    //ca50 = 0.7;
    Tref = 120;
    nu = 7;
    mu = 3;
}

double k_ws = 0.004;
double k_uw = 0.026;

double lambda_min = 0.87;
double lambda_max = 1.2;

//double dydt[6] = {0,0,0,0,0,0};

k_ws = k_ws*mu;
k_uw = k_uw*nu;

double cdw = phi*k_uw*(1-dr)*(1-wfrac)/((1-dr)*wfrac);
double cds = phi*k_ws*(1-dr)*wfrac/dr;
double k_wu = k_uw*(1/wfrac-1)-k_ws;
double k_su = k_ws*(1/dr-1)*wfrac; 
double A = (0.25*TOT_A)/((1-dr)*wfrac+dr)*(dr/0.25);

//double lambda0 = (lambda_max < lambda) ? lambda_max : lambda;
double lambda0 = fminf(lambda_max, lambda);
//double lambda_aux = (lambda_min < lambda0) ? lambda_min : lambda0;
//aux = 1+beta_0*(lambda0+lambda_aux-(1+lambda_min));
//double Lfac = (aux > 0) ? aux : 0;
double Lfac = fmaxf(0,1+beta_0*(lambda0+fminf(lambda_min,lambda0)-(1+lambda_min)));

double XU = (1-TmBlocked)-XW-XS; // unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
double xb_ws = k_ws*XW;
double xb_uw = k_uw*XU;
double xb_wu = k_wu*XW;
double xb_su = k_su*XS;

//aux = ((ZETAS>0)*ZETAS > (ZETAS<-1)*(-ZETAS-1)) ? (ZETAS>0)*ZETAS : (ZETAS<-1)*(-ZETAS-1); 
double gamma_rate=gamma*fmaxf((ZETAS>0)*ZETAS,(ZETAS<-1)*(-ZETAS-1));
double xb_su_gamma=gamma_rate*XS;
double gamma_rate_w=gamma_wu*abs(ZETAW); // weak xbs don't like being strained
double xb_wu_gamma=gamma_rate_w*XW;

//dydt[0] = xb_ws-xb_su-xb_su_gamma;
//dydt[1] = xb_uw-xb_wu-xb_ws-xb_wu_gamma;
double dXS = xb_ws-xb_su-xb_su_gamma;
double dXW = xb_uw-xb_wu-xb_ws-xb_wu_gamma;

//aux = (lambda-1 < 0.2) ? (lambda-1) : 0.2;
//ca50 = ca50+beta_1*aux;
ca50=ca50+beta_1*fminf(0.2,lambda-1);
//dydt[2] = koff*(pow(((cai*1000)/ca50),TRPN_n)*(1-Ca_TRPN)-Ca_TRPN); // untouched
double dCa_TRPN = koff*(pow(((cai*1000)/ca50),TRPN_n)*(1-Ca_TRPN)-Ca_TRPN);

double XSSS = dr*0.5;
double XWSS = (1-dr)*wfrac*0.5;
double ktm_block = ktm_unblock*(pow(perm50,nperm))*0.5/(0.5-XSSS-XWSS);

//aux = pow(Ca_TRPN,-(nperm/2));
//aux2 = pow(Ca_TRPN,(nperm/2));
//aux3 = (100 < aux) ? 100 : aux;
//dydt[3] = ktm_block*aux3*XU-ktm_unblock*aux2*TmBlocked;
//double dTmBlocked = ktm_block*aux3*XU-ktm_unblock*aux2*TmBlocked;
double dTmBlocked = ktm_block*fminf(100, pow(Ca_TRPN,-(nperm/2)))*XU - ktm_unblock*( pow(Ca_TRPN,(nperm/2)))*TmBlocked;

// velocity dependence -- assumes distortion resets on W->S
//dydt[4] = A*lambda_rate-cds*ZETAS;    // - gamma_rate * ZETAS;
//dydt[5] = A*lambda_rate-cdw*ZETAW;    // - gamma_rate_w * ZETAW;
double dZETAS = A*lambda_rate-cds*ZETAS;    // - gamma_rate * ZETAS;
double dZETAW = A*lambda_rate-cdw*ZETAW;    // - gamma_rate_w * ZETAW;

// Active Force
double Ta = Lfac*(Tref/dr)*((ZETAS+1)*XS+(ZETAW)*XW);

active = Ta;


// physical constants
double R=8314.0;
double T=310.0;
double F=96485.0;

// cell geometry
double L=0.01;
double rad=0.0011;
double vcell=1000*3.14*rad*rad*L;
double Ageo=2*3.14*rad*rad+2*3.14*rad*L;
double Acap=2*Ageo;
double vmyo=0.68*vcell;
double vnsr=0.0552*vcell;
double vjsr=0.0048*vcell;
double vss=0.02*vcell;

double fkatp = 0.0;
double gkatp = 4.3195;

// CaMK constants
double KmCaMK=0.15;
double aCaMK=0.05*aCaMK_Multiplier;
double bCaMK=0.00068;
double CaMKo=0.05;
double KmCaM=0.0015;
// update CaMK
double CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
double CaMKa=CaMKb+CaMKt;
double dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;      // Euler

// reversal potentials
double ENa=(R*T/F)*log(nao/nai);
double EK=(R*T/F)*log(ko/ki);
double PKNa=0.01833;
double EKs=(R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));

// convenient shorthand calculations
double vffrt=v*F*F/(R*T);
double vfrt=v*F/(R*T);
double frt = F/(R*T);

double K_o_n = 5.0;
double A_atp = 2.0;
double K_atp = 0.25;
double akik = pow((ko / K_o_n), 0.24);
double bkik = (1.0 / (1.0 + pow((A_atp / K_atp), 2.0)));

double fINap=(1.0/(1.0+KmCaMK/CaMKa));
double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
double fItop=(1.0/(1.0+KmCaMK/CaMKa));
double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));

// INa formulations
// The Grandi implementation updated with INa phosphorylation.
// m gate
double mss = 1 / (pow(1 + exp( -(56.86 + v) / 9.03 ),2));
double taum = 0.1292 * exp(-pow((v+45.79)/15.54,2)) + 0.06487 * exp(-pow((v-4.823)/51.12,2));
double dm = (mss - m) / taum;                     // Rush-Larsen

// h gate
double ah = (v >= -40) ? (0) : (0.057 * exp( -(v + 80) / 6.8 ));
double bh = (v >= -40) ? (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) : ((2.7 * exp( 0.079 * v) + 3.1*pow(10,5) * exp(0.3485 * v)));
double tauh = 1 / (ah + bh);
double hss = 1 / (pow(1 + exp( (v + 71.55)/7.43 ),2));
double dh = (hss - h) / tauh;                     // Rush-Larsen
// j gate
double aj = (v >= -40) ? (0) : (((-2.5428 * pow(10,4)*exp(0.2444*v) - 6.948*pow(10,-6) * exp(-0.04391*v)) * (v + 37.78)) / (1 + exp( 0.311 * (v + 79.23) )));
double bj = (v >= -40) ? ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) : ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )));
double tauj = 1 / (aj + bj);
double jss = 1 / pow((1 + exp( (v + 71.55)/7.43 )),2);
double dj = (jss - j) / tauj;                     // Rush-Larsen

// h phosphorylated
double hssp = 1 / pow((1 + exp( (v + 71.55 + 6)/7.43 )),2);
double dhp = (hssp - hp) / tauh;                  // Rush-Larsen
// j phosphorylated
double taujp = 1.46 * tauj;
double djp = (jss - jp) / taujp;                  // Rush-Larsen
double GNa = 11.7802;
double INa = INa_Multiplier*GNa*(v-ENa)*pow(m,3.0)*((1.0-fINap)*h*j+fINap*hp*jp);

// INaL
// calculate INaL
double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
double tm = 0.1292 * exp(-pow(((v+45.79)/15.54),2)) + 0.06487 * exp(-pow(((v-4.823)/51.12),2)); 
double tmL=tm;
double dmL=(mLss-mL)/tmL;                                         // Rush-Larsen
double hLss=1.0/(1.0+exp((v+87.61)/7.488));
double thL=200.0;
double dhL=(hLss-hL)/thL;                                         // Rush-Larsen
double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
double thLp=3.0*thL;
double dhLp=(hLssp-hLp)/thLp;                                     // Rush-Larsen
double GNaL=0.0279*INaL_Multiplier;
if (celltype==EPI) GNaL=GNaL*0.6;
double INaL = GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

// ITo
// calculate Ito
double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
double da=(ass-a)/ta;                                                         // Rush-Larsen
double iss=1.0/(1.0+exp((v+43.94)/5.711));
double delta_epi = (celltype == EPI) ? 1.0-(0.95/(1.0+exp((v+70.0)/5.0))) : 1.0;
double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
tiF=tiF*delta_epi;
tiS=tiS*delta_epi;
double AiF=1.0/(1.0+exp((v-213.6)/151.2));
double AiS=1.0-AiF;
double diF=(iss-iF)/tiF;                                                      // Rush-Larsen
double diS=(iss-iS)/tiS;                                                      // Rush-Larsen
double i=AiF*iF+AiS*iS;
double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
double dap=(assp-ap)/ta;                                                     // Rush-Larsen
double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
double tiFp=dti_develop*dti_recover*tiF;
double tiSp=dti_develop*dti_recover*tiS;
double diFp=(iss-iFp)/tiFp;                                                   // Rush-Larsen
double diSp=(iss-iSp)/tiSp;                                                   // Rush-Larsen
double ip=AiF*iFp+AiS*iSp;
double Gto=0.16*Ito_Multiplier;
Gto = (celltype == EPI || celltype == MCELL) ? Gto*2.0 : Gto;
double Ito = Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

// ICaL
// a variant updated by jakub, using a changed activation curve
// it computes both ICaL in subspace and myoplasm (_i)

// calculate ICaL, ICaNa, ICaK
double dss=1.0763*exp(-1.0070*exp(-0.0829*(v)));  // magyar
if(v >31.4978) dss = 1; // activation cannot be greater than 1
double td= 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
double dd=(dss-d)/td;                                                                     // Rush-Larsen
double fss=1.0/(1.0+exp((v+19.58)/3.696));
double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
double Aff=0.6;
double Afs=1.0-Aff;
double dff=(fss-ff)/tff;                                                                  // Rush-Larsen
double dfs=(fss-fs)/tfs;                                                                  // Rush-Larsen
double f=Aff*ff+Afs*fs;
double fcass=fss;
double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));

double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
double Afcas=1.0-Afcaf;
double dfcaf=(fcass-fcaf)/tfcaf;                                                          // Rush-Larsen
double dfcas=(fcass-fcas)/tfcas;                                                          // Rush-Larsen
double fca=Afcaf*fcaf+Afcas*fcas;

double tjca = 75;
//double tjca = 72.5;
double jcass = 1.0/(1.0+exp((v+18.08)/(2.7916)));   
double djca=(jcass-jca)/tjca;                                                                  // Rush-Larsen
double tffp=2.5*tff;
double dffp=(fss-ffp)/tffp;                                                                    // Rush-Larsen
double fp=Aff*ffp+Afs*fs;
double tfcafp=2.5*tfcaf;
double dfcafp=(fcass-fcafp)/tfcafp;                                                           // Rush-Larsen
double fcap=Afcaf*fcafp+Afcas*fcas;

// SS nca
double Kmn=0.002;
double k2n=500.0;
double km2n=jca*1;
double anca=1.0/(k2n/km2n+pow((1.0+Kmn/cass),4.0));
double dnca=anca*k2n-nca*km2n;                                                        // Euler

// myoplasmic nca
double anca_i = 1.0/(k2n/km2n+pow((1.0+Kmn/cai),4.0));
double dnca_i = anca_i*k2n-nca_i*km2n;                                                     // Euler

// SS driving force
//double clo = 150;     // Extracellular Cl  [mM]
//double cli = 24.0;    // Intracellular Cl  [mM]
double Io = 0.5*(nao + ko + clo + 4*cao)/1000;         // ionic strength outside. /1000 is for things being in micromolar
double Ii = 0.5*(nass + kss + cli + 4*cass)/1000;     // ionic strength outside. /1000 is for things being in micromolar
//double Ii = 0.5*(nass + kss + clss + 4*cass)/1000;     // (dynCl) ionic strength outside. /1000 is for things being in micromolar
// The ionic strength is too high for basic DebHuc. We'll use Davies
double dielConstant = 74;     // water at 37 degrees.
double temp = 310;            // body temp in kelvins.
double constA = 1.82*pow(10,6)*pow((dielConstant*temp),(-1.5));

double gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
double gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
double gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
double gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
double gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
double gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

double PhiCaL_ss =  4.0*vffrt*(gamma_cai*cass*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
double PhiCaNa_ss =  1.0*vffrt*(gamma_nai*nass*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
double PhiCaK_ss =  1.0*vffrt*(gamma_ki*kss*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);

// Myo driving force
Io = 0.5*(nao + ko + clo + 4*cao)/1000; // ionic strength outside. /1000 is for things being in micromolar
Ii = 0.5*(nai + ki + cli + 4*cai)/1000; // ionic strength outside. /1000 is for things being in micromolar
// The ionic strength is too high for basic DebHuc. We'll use Davies

gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

double gammaCaoMyo = gamma_cao;
double gammaCaiMyo = gamma_cai;

double PhiCaL_i =  4.0*vffrt*(gamma_cai*cai*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
double PhiCaNa_i =  1.0*vffrt*(gamma_nai*nai*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
double PhiCaK_i =  1.0*vffrt*(gamma_ki*ki*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);
// The rest
double PCa=8.3757e-05 * ICaL_Multiplier;
if (celltype==EPI)
    PCa=PCa*1.2;
else if (celltype==MCELL)
    //PCa=PCa*2;
    PCa = PCa*1.8;

double PCap=1.1*PCa;
double PCaNa=0.00125*PCa;
double PCaK=3.574e-4*PCa;
double PCaNap=0.00125*PCap;
double PCaKp=3.574e-4*PCap;

double ICaL_ss=(1.0-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
double ICaNa_ss=(1.0-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
double ICaK_ss=(1.0-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);

double ICaL_i=(1.0-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
double ICaNa_i=(1.0-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
double ICaK_i=(1.0-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

// And we weight ICaL (in ss) and ICaL_i
ICaL_i = ICaL_i * (1-ICaL_fractionSS);
ICaNa_i = ICaNa_i * (1-ICaL_fractionSS);
ICaK_i = ICaK_i * (1-ICaL_fractionSS);
ICaL_ss = ICaL_ss * ICaL_fractionSS;
ICaNa_ss = ICaNa_ss * ICaL_fractionSS;
ICaK_ss = ICaK_ss * ICaL_fractionSS;

double ICaL = ICaL_ss + ICaL_i;
double ICaNa = ICaNa_ss + ICaNa_i;
double ICaK = ICaK_ss + ICaK_i;
double ICaL_tot = ICaL + ICaNa + ICaK;

// Ikr
// Variant based on Lu-Vandenberg
// Extracting state vector
// IMPORTANT: For reason of backward compatibility of naming of an older version of a MM IKr
//            c3 in code is c0 in article diagram, c2 is c1, c1 is c2
double c0 = ikr_c0;
double c1 = ikr_c1;
double c2 = ikr_c2;
//double c0 = ikr_c2;
//double c1 = ikr_c1;
//double c2 = ikr_c0;
double o = ikr_o;
double I = ikr_i;
double b = 0; // no channels blocked in via the mechanism of specific MM states

// transition rates
// from c0 to c1 in l-v model,
double alpha = 0.1161 * exp(0.2990 * vfrt);
// from c1 to c0 in l-v/
double beta =  0.2442 * exp(-1.604 * vfrt);

// from c1 to c2 in l-v/
double alpha1 = 1.25 * 0.1235;
// from c2 to c1 in l-v/
double beta1 =  0.1911;

// from c2 to o/           c1 to o
double alpha2 =0.0578 * exp(0.9710 * vfrt);
// from o to c2/
double beta2 = 0.349e-3* exp(-1.062 * vfrt);

// from o to i
double alphai = 0.2533 * exp(0.5953 * vfrt);
// from i to o
double betai = 1.25* 0.0522 * exp(-0.8209 * vfrt);

// from c2 to i (from c1 in orig)
double alphac2ToI = 0.52e-4 * exp(1.525 * vfrt);
// from i to c2
double betaItoC2 = 0.85e-8 * exp(-1.842 * vfrt);
betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai);
// transitions themselves
// for reason of backward compatibility of naming of an older version of a
// MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

double dc0 = c1 * beta - c0 * alpha;                   // Euler
double dc1 = c0 * alpha + c2*beta1 - c1*(beta+alpha1); // Euler
double dc2 = c1 * alpha1 + o*beta2 + I*betaItoC2 - c2 * (beta1 + alpha2 + alphac2ToI); // Euler
double delta_o = c2 * alpha2 + I*betai - o*(beta2+alphai);    // Euler
double di = c2*alphac2ToI + o*alphai - I*(betaItoC2 + betai); // Euler

double GKr = 0.0321 * sqrt(ko/5) * IKr_Multiplier; // 1st element compensates for change to ko (sqrt(5/5.4)* 0.0362)
if (celltype==EPI)
    GKr=GKr*1.3;
else if (celltype==MCELL)
    GKr=GKr*0.8;

double IKr = GKr * o  * (v-EK);

// calculate IKs
double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
double dxs1=(xs1ss-xs1)/txs1;                              // Rush-Larsen
double xs2ss=xs1ss;
double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
double dxs2=(xs2ss-xs2)/txs2;                               // Rush-Larsen
double KsCa=1.0+0.6/(1.0+pow((3.8e-5/cai),1.4));
double GKs= 0.0011*IKs_Multiplier;
if (celltype==EPI)
    GKs=GKs*1.4;
double IKs = GKs*KsCa*xs1*xs2*(v-EKs);

// IK1
double aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)));
double bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)));
double K1ss = aK1/(aK1+bK1);
double GK1=IK1_Multiplier  * 0.6992; // 0.7266; // * sqrt(5/5.4))
if (celltype==EPI)
    GK1=GK1*1.2;
else if (celltype==MCELL)
    GK1=GK1*1.3;
double IK1=GK1*sqrt(ko/5)*K1ss*(v-EK);

// IKCa
double fIKCass = 0.8;
double kdikca = 6.05e-4;
double ikcan = 3.5;
double GKCa = 0.003 * IKCa_Multiplier;
double IKCa_ss = GKCa * fIKCass * pow(cass,ikcan) / (pow(cass,ikcan) + pow(kdikca,ikcan)) * (v-EK);
double IKCa_i = GKCa * (1.0-fIKCass) * pow(cai,ikcan) / (pow(cai,ikcan) + pow(kdikca,ikcan)) * (v-EK);
double IKCa = IKCa_ss + IKCa_i;

// INaCa
double zca = 2.0;
double kna1=15.0;
double kna2=5.0;
double kna3=88.12;
double kasymm=12.5;
double wna=6.0e4;
double wca=6.0e4;
double wnaca=5.0e3;
double kcaon=1.5e6;
double kcaoff=5.0e3;
double qna=0.5224;
double qca=0.1670;
double hca=exp((qca*v*F)/(R*T));
double hna=exp((qna*v*F)/(R*T));

double h1=1+nai/kna3*(1+hna);
double h2=(nai*hna)/(kna3*h1);
double h3=1.0/h1;
double h4=1.0+nai/kna1*(1+nai/kna2);
double h5=nai*nai/(h4*kna1*kna2);
double h6=1.0/h4;
double h7=1.0+nao/kna3*(1.0+1.0/hna);
double h8=nao/(kna3*hna*h7);
double h9=1.0/h7;
double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
double h11=nao*nao/(h10*kna1*kna2);
double h12=1.0/h10;
double k1=h12*cao*kcaon;
double k2=kcaoff;
double k3p=h9*wca;
double k3pp=h8*wnaca;
double k3=k3p+k3pp;
double k4p=h3*wca/hca;
double k4pp=h2*wnaca;
double k4=k4p+k4pp;
double k5=kcaoff;
double k6=h6*cai*kcaon;
double k7=h5*h2*wna;
double k8=h8*h11*wna;
double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
double E1=x1/(x1+x2+x3+x4);
double E2=x2/(x1+x2+x3+x4);
double E3=x3/(x1+x2+x3+x4);
double E4=x4/(x1+x2+x3+x4);
double KmCaAct=150.0e-6;
double allo=1.0/(1.0+pow((KmCaAct/cai),2.0));
double zna=1.0;
double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
double JncxCa=E2*k2-E1*k1;
double Gncx= 0.0034 * INaCa_Multiplier;
if (celltype==EPI)
    Gncx=Gncx*1.1;
else if (celltype==MCELL)
    Gncx=Gncx*1.4;
double INaCa_i = (1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

// calculate INaCa_ss
h1=1+nass/kna3*(1+hna);
h2=(nass*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+nass/kna1*(1+nass/kna2);
h5=nass*nass/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+nao/kna3*(1.0+1.0/hna);
h8=nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
h11=nao*nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*cass*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150.0e-6 ;
allo=1.0/(1.0+pow((KmCaAct/cass),2.0));
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
double INaCa_ss = INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);

// calculate INaK
double k1p=949.5;
double k1m=182.4;
double k2p=687.2;
double k2m=39.4;
k3p=1899.0;
double k3m=79300.0;
k4p=639.0;
double k4m=40.0;
double Knai0=9.073;
double Knao0=27.78;
double delta=-0.1550;
double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
double Kki=0.5;
double Kko=0.3582;
double MgADP=0.05;
double MgATP=9.8;
double Kmgatp=1.698e-7;
double H=1.0e-7;
double eP=4.2;
double Khp=1.698e-7;
double Knap=224.0;
double Kxkur=292.0;
double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
double a1=(k1p*pow((nai/Knai),3.0))/(pow((1.0+nai/Knai),3.0)+pow((1.0+ki/Kki),2.0)-1.0);
double b1=k1m*MgADP;
double a2=k2p;
double b2=(k2m*pow((nao/Knao),3.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
double a3=(k3p*pow((ko/Kko),2.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
double b4=(k4m*pow((ki/Kki),2.0))/(pow((1.0+nai/Knai),3.0)+pow((1.0+ki/Kki),2.0)-1.0);
x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
double zk=1.0;
double JnakNa=3.0*(E1*a3-E2*b3);
double JnakK=2.0*(E4*b1-E3*a1);
double Pnak= 15.4509;
if (celltype==EPI)
    Pnak=Pnak*0.9;
else if (celltype==MCELL)
    Pnak=Pnak*0.7;

double INaK = Pnak*(zna*JnakNa+zk*JnakK)*INaK_Multiplier;

// Minor/background currents
// calculate IKb
double xkb=1.0/(1.0+exp(-(v-10.8968)/(23.9871)));
double GKb=0.0189*IKb_Multiplier;

if (IKCa_Multiplier > 0.0)
    GKb = GKb*0.9;
if (celltype==EPI)
    GKb=GKb*0.6;
double IKb = GKb*xkb*(v-EK);

// calculate INab
double PNab=1.9239e-09*INab_Multiplier;
double INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

// calculate ICab
double PCab=5.9194e-08*ICab_Multiplier; 
double ICab=PCab*4.0*vffrt*(gammaCaiMyo*cai*exp(2.0*vfrt)-gammaCaoMyo*cao)/(exp(2.0*vfrt)-1.0);

// calculate IpCa
double GpCa=5e-04*IpCa_Multiplier;
double IpCa=GpCa*cai/(0.0005+cai);

// Chloride
// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

double ECl = (R*T/F)*log(cli/clo);            // [mV]
//double EClss = (R*T/F)*log(clss/clo);         // [mV]   dynCl

double Fjunc = 1;   
double Fsl = 1-Fjunc; // fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace

double GClCa = ICaCl_Multiplier * 0.2843;   // [mS/uF]
double GClB = IClb_Multiplier * 1.98e-3;    // [mS/uF] 
double KdClCa = 0.1;                        // [mM]

//double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-EClss);     // dynCl
//double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ECl);            // dynCl
double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-ECl);
double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ECl);

double I_ClCa = I_ClCa_junc+I_ClCa_sl;
double I_Clbk = GClB*(v-ECl);

// Calcium handling
// calculate ryanodione receptor calcium induced calcium release from the jsr
double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));

// Jrel
double jsrMidpoint = 1.7;

double bt=4.75;
double a_rel=0.5*bt;
double Jrel_inf=a_rel*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8.0));
if (celltype==MCELL)
    Jrel_inf=Jrel_inf*1.7;
double tau_rel=bt/(1.0+0.0123/cajsr);

if (tau_rel<0.001)
    tau_rel=0.001;

double dJrelnp=(Jrel_inf-Jrel_np)/tau_rel;                     // Rush-Larsen

double btp=1.25*bt;
double a_relp=0.5*btp;
double Jrel_infp=a_relp*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8.0));
if (celltype==MCELL)
    Jrel_infp=Jrel_infp*1.7;
double tau_relp=btp/(1.0+0.0123/cajsr);

if (tau_relp<0.001)
    tau_relp=0.001;

tau_relp = tau_relp*taurelp_Multiplier;
double dJrelp=(Jrel_infp-Jrel_p)/tau_relp;                     // Rush-Larsen
double Jrel = Jrel_Multiplier * 1.5378 * ((1.0-fJrelp)*Jrel_np+fJrelp*Jrel_p);

double fJupp=(1.0/(1.0+KmCaMK/CaMKa));

// Jup
// calculate serca pump, ca uptake flux
// camkFactor = 2.4;
// gjup = 0.00696;
// Jupnp=Jup_Multiplier * gjup*cai/(cai+0.001);
// Jupp=Jup_Multiplier * camkFactor*gjup*cai/(cai + 8.2500e-04);
// if celltype==1
//     Jupnp=Jupnp*1.3;
//     Jupp=Jupp*1.3;
// end
// 
// 
// Jleak=Jup_Multiplier * 0.00629 * cansr/15.0;
// Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

// calculate serca pump, ca uptake flux
double Jupnp = Jup_Multiplier * 0.005425*cai/(cai+0.00092);
double Jupp = Jup_Multiplier * 2.75*0.005425*cai/(cai+0.00092-0.00017);
if (celltype==EPI) {
    Jupnp=Jupnp*1.3;
    Jupp=Jupp*1.3;
}
double Jleak=Jup_Multiplier* 0.0048825*cansr/15.0;
double Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

//calculate tranlocation flux
double Jtr=(cansr-cajsr)/60;

// I_katp current (fkatp current)
//double I_katp = (fkatp * gkatp * akik * bkik * (v - EK));

// stimulus current
double Istim = calc_I_stim;

//update the membrane voltage
//double dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+Istim);    // Euler
double dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+IKCa+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+Istim);    // Euler
//double dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+I_katp+Istim);    // Euler

// calculate diffusion fluxes
double JdiffNa=(nass-nai)/2.0;
double JdiffK=(kss-ki)/2.0;
double Jdiff=(cass-cai)/0.2;
//double JdiffCl=(clss-cli)/2.0;  // dynCl

// calcium buffer constants 
double cmdnmax= 0.05; 
if (celltype==EPI)
    cmdnmax=cmdnmax*1.3;
double kmcmdn=0.00238; 
//double trpnmax=0.07;
double kmtrpn=0.0005;
double BSRmax=0.047;
double KmBSR = 0.00087;
double BSLmax=1.124;
double KmBSL = 0.0087;
double csqnmax=10.0;
double kmcsqn=0.8;

// update intracellular concentrations, using buffers for cai, cass, cajsr
double dnai=-(ICaNa_i+INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;     // Euler
double dnass=-(ICaNa_ss+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;                                   // Euler

double dki=-(ICaK_i+Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;        // Euler
double dkss=-(ICaK_ss)*Acap/(F*vss)-JdiffK;                                                   // Euler

//double Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0)+trpnmax*kmtrpn/pow((kmtrpn+cai),2.0)); // dynCl
double Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0));
//double dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo); // dynCl
double dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo-dCa_TRPN*trpnmax);

double Bcass=1.0/(1.0+BSRmax*KmBSR/pow((KmBSR+cass),2.0)+BSLmax*KmBSL/pow((KmBSL+cass),2.0));
double dcass=Bcass*(-(ICaL_ss-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

double dcansr=Jup-Jtr*vjsr/vnsr;

double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow((kmcsqn+cajsr),2.0));
double dcajsr=Bcajsr*(Jtr-Jrel);

//double dcli  = - (I_Clbk + I_ClCa_sl)*Acap/(-1*F*vmyo)+JdiffCl*vss/vmyo;        // dynCl
//double dclss = - I_ClCa_junc*Acap/(-1*F*vss)-JdiffCl;                           // dynCl

// Right-hand side
rDY[0]  = dv;
rDY[1]  = dnai;
rDY[2]  = dnass;
rDY[3]  = dki;
rDY[4]  = dkss;
rDY[5]  = dcai;
rDY[6]  = dcass;
rDY[7]  = dcansr;
rDY[8]  = dcajsr;
rDY[9]  = dm;
rDY[10] = dhp;
rDY[11] = dh;
rDY[12] = dj;
rDY[13] = djp;
rDY[14] = dmL;
rDY[15] = dhL;
rDY[16] = dhLp;
rDY[17] = da;
rDY[18] = diF;
rDY[19] = diS;
rDY[20] = dap;
rDY[21] = diFp;
rDY[22] = diSp;
rDY[23] = dd;
rDY[24] = dff;
rDY[25] = dfs;
rDY[26] = dfcaf;
rDY[27] = dfcas;
rDY[28] = djca;
rDY[29] = dnca;
rDY[30] = dnca_i;
rDY[31] = dffp;
rDY[32] = dfcafp;
rDY[33] = dxs1;
rDY[34] = dxs2;
rDY[35] = dJrelnp;
rDY[36] = dCaMKt;
rDY[37] = dc0;
rDY[38] = dc1;
rDY[39] = dc2;
rDY[40] = delta_o;
rDY[41] = di;
rDY[42] = dJrelp;
// -----------------------
// Land-Niederer
rDY[43] = dXS;
rDY[44] = dXW;
rDY[45] = dCa_TRPN;
rDY[46] = dTmBlocked;
rDY[47] = dZETAS;
rDY[48] = dZETAW;
rDY[49] = Ta;

}



