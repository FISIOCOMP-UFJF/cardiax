#include "rice_ord.hpp"

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) < (b) ? (b) : (a))
#define sign(a) ((a) < (0.) ? (-1.) : (1.))
#define heav(a) ((a) < (0.) ? (0.) : (1.))
#define abs(a) ((a) < (0.) ? (-a) : (a))

//#define HF

RiceOHaraRudy::RiceOHaraRudy() : CellModel(52)
{
  // state variable names
  const string vn[] = {
      // O'Hara Rudy vars
      "V", "nai", "nass", "ki", "kss", "cai", "cass", "cansr", "cajsr",
      "m", "hf", "hs", "j", "hsp", "jp", "mL", "hL", "hLp", "a", "iF",
      "is", "ap", "iFp", "iSp", "d", "ff", "fs", "fcaf", "fcas", "jca",
      "nca", "ffp", "fcafp", "xrf", "xrs", "xs1", "xs2", "xk1",
      "Jrelnp", "Jrelp", "CaMKt",
      // Rice vars
      "TRPNCaL", "TRPNCaH", "N_NoXB", "P_NoXB", "N", "XBprer",
      "XBpostr", "xXBprer", "xXBpostr", "SL", "intf"
  };

  for(int i=0; i<51; i++)
    var_names.insert( std::pair<int, std::string>(i,  vn[i]) );

  // RL gating variables  
  for(int i=9;  i<=39; i++) rlvars.insert(i);
  
  // Initialize monitored values
  monitored.push_back( &active_force );  
  active_force = 0.0;
}

void RiceOHaraRudy::init(double * values) const
{
  values[ 0] = -87.5;     //v=-87.5;
  values[ 1] = 7;         //nai=7;
  values[ 2] = values[1]; //nass=nai;
  values[ 3] = 145;       //ki=145;
  values[ 4] = values[3]; //kss=ki;
  values[ 5] = 1.0e-4;    //cai=1.0e-4;
  values[ 6] = values[5]; //cass=cai;
  values[ 7] = 1.2;       //cansr=1.2;
  values[ 8] = values[7]; //cajsr=cansr;  
  values[ 9] = 0;         //m=0;
  values[10] = 1;         //hf=1;
  values[11] = 1;         //hs=1;
  values[12] = 1;         //j=1;
  values[13] = 1;         //hsp=1;
  values[14] = 1;         //jp=1;
  values[15] = 0;         //mL=0;
  values[16] = 1;         //hL=1;
  values[17] = 1;         //hLp=1;
  values[18] = 0;         //a=0;
  values[19] = 1;         //iF=1;
  values[20] = 1;         //iS=1;
  values[21] = 0;         //ap=0;
  values[22] = 1;         //iFp=1;
  values[23] = 1;         //iSp=1;
  values[24] = 0;         //d=0;
  values[25] = 1;         //ff=1;
  values[26] = 1;         //fs=1;
  values[27] = 1;         //fcaf=1;
  values[28] = 1;         //fcas=1;
  values[29] = 1;         //jca=1;
  values[30] = 0;         //nca=0;
  values[31] = 1;         //ffp=1;
  values[32] = 1;         //fcafp=1;
  values[33] = 0;         //xrf=0;
  values[34] = 0;         //xrs=0;
  values[35] = 0;         //xs1=0;
  values[36] = 0;         //xs2=0;
  values[37] = 1;         //xk1=1;
  values[38] = 0;         //Jrelnp=0;
  values[39] = 0;         //Jrelp=0;
  values[40] = 0;         //CaMKt=0;
  
  // Rice variables
  values[41] = 1.311635148577327e-02;  // TRPNCaL - Rice myofilament variables
  values[42] = 0.0001;                 // TRPNCaH
  values[43] = 0.99;                   // N_NoXB
  values[44] = 0.01;                   // P_NoXB
  values[45] = 0.97;                   // N
  values[46] = 2.211689363408625e-04;  // XBprer
  values[47] = 1.079302937362542e-04;  // XBpostr
  values[48] = 1.007048444140799e-18;  // xXBprer
  values[49] = 7.000000000000143e-03;  // xXBpostr
  values[50] = 2.2;                    // SL
  values[51] = 0.0;                    // intf  
}

void RiceOHaraRudy::equation(const double time, const double * statevars, 
                             double * values)
{  
  // ---------------------------------------------------------------------------
  // ORd parameters
  // ---------------------------------------------------------------------------
  
  double const nao=140.0; // extracellular sodium in mM
  double const cao=1.8;   // extracellular calcium in mM
  double const ko=5.4;    // extracellular potassium in mM

  // buffer paramaters
  double const BSRmax=0.047;
  double const KmBSR=0.00087;
  double const BSLmax=1.124;
  double const KmBSL=0.0087;
  double const cmdnmax=0.05;  
  double const kmcmdn=0.00238;
  double const trpnmax=0.07;
  double const kmtrpn=0.0005;
  double const csqnmax=10.0;
  double const kmcsqn=0.8;

  // CaMK paramaters
  double const aCaMK=0.05;
  double const bCaMK=0.00068;
  double const CaMKo=0.05;
  double const KmCaM=0.0015;
  double const KmCaMK=0.15;

  // physical constants
  double const R=8314.0;
  double const T=310.0;
  double const F=96485.0;

  // cell geometry
  double const L=0.01;
  double const rad=0.0011;
  double const vcell=1000*3.14*rad*rad*L;
  double const Ageo=2*3.14*rad*rad+2*3.14*rad*L;
  double const Acap=2*Ageo;
  double const vmyo=0.68*vcell;
  double const vmito=0.26*vcell;
  double const vsr=0.06*vcell;
  double const vnsr=0.0552*vcell;
  double const vjsr=0.0048*vcell;
  double const vss=0.02*vcell;

  // introduce varaibles for reversal potentials, currents, fluxes, and CaMK
  double ENa,EK,EKs;
  double INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist;
  double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
  double CaMKa,CaMKb;

  // GAMBIARRA for RUSH-LARSEN
  double dt=0.05; 
 
  // ---------------------------------------------------------------------------
  // Rice parameters
  // ---------------------------------------------------------------------------
  // Sarcomere Geometry
  const double SLmax   = 2.4;        // (um) maximum sarcomere length
  const double SLmin   = 1.4;        // (um) minimum sarcomere length
  const double len_thin  = 1.2;      // (um) thin filament length
  const double len_thick = 1.65;     // (um) thick filament length
  const double len_hbare = 0.1;      // (um) length of bare portion of thick filament

  // Temperature Dependence
  const double Qkon = 1.5;
  const double Qkoff = 1.3;
  const double Qkn_p = 1.6;
  const double Qkp_n = 1.6;
  const double Qfapp = 6.25;
  const double Qgapp = 2.5;
  const double Qhf = 6.25;
  const double Qhb = 6.25;
  const double Qgxb = 6.25;

  // *** RAT PARAMETERS ***
  // Ca binding to troponin
//  double kon     = 50e-3;      // (1/[ms uM])
//  double koffL   = 250e-3;     // (1/ms)
//  double koffH   = 25e-3;      // (1/ms)
//  double perm50  = 0.5;        // perm variable that controls n to p transition
//  double nperm   = 15;         //   in Hill-like fashion
//  double kn_p    = 500e-3;     // (1/ms)
//  double kp_n    = 50e-3;      // (1/ms)
//  double koffmod = 1.0;        // mod to change species
//
//  double fapp    = 500e-3;     // (1/ms) XB on rate
//  double gapp    = 70e-3;      // (1/ms) XB off rate
//  double gslmod  = 6;          // controls SL effect on gapp
//  double rice_hf = 2000e-3;    // (1/ms) rate between pre-force and force states
//  double hfmdc   = 5;          //
//  double hb      = 400e-3;     // (1/ms) rate between pre-force and force states
//  double hbmdc   = 0;          //
//  double gxb     = 70e-3;      // (1/ms) ATP consuming transition rate
//  double sigmap  = 8;          // distortion dependence of STP using transition gxb
//  double sigman  = 1;          //
//  double xbmodsp = 1.0;        // mouse specific modification for XB cycling rates
  
  // *** ADJUSTED PARAMETERS ***
  // 0.428999 0.00190106 0.697061 0.0779171 0.415683 0.501076 0.0547643
  //
  double kon     = 50e-3;      // (1/[ms uM])
  double koffL   = 250e-3;     // (1/ms)
  double koffH   = 25e-3;      // (1/ms)
  double perm50  = 0.5;        // perm variable that controls n to p transition
  double nperm   = 15;         //   in Hill-like fashion
  double koffmod = 1.0;        // mod to change species
  
  double kn_p    = 0.428999;    // (1/ms)
  double kp_n    = 0.00190106;  // (1/ms)
  double fapp    = 0.697061;    // (1/ms) XB on rate
  double gapp    = 0.0779171;   // (1/ms) XB off rate
  double rice_hf = 0.415683;    // (1/ms) rate between pre-force and force states
  double hb      = 0.501076;    // (1/ms) rate between pre-force and force states
  double gxb     = 0.0547643;   // (1/ms) ATP consuming transition rate
  
  double gslmod  = 6;          // controls SL effect on gapp
  double hfmdc   = 5;          //
  double hbmdc   = 0;          //  
  double sigmap  = 8;          // distortion dependence of STP using transition gxb
  double sigman  = 1;          //
  double xbmodsp = 1.0;        // mouse specific modification for XB cycling rates
  

  // Mean strain of strongly-bound states
  const double x_0     = 0.007;      // (um) strain induced by head rotation
  const double xPsi    = 2;          // scaling factor balancing SL motion and XB cycling

  // Normalized active and passive force
  const double SLrest  = 1.85;       // (um) rest SL length for 0 passive force
  const double PCon_t  = 0.002;      // (norm Force) passive force due to titin
  const double PExp_t  = 10;         //   these apply to trabeculae and single cells only
  const double SL_c    = 2.25;       // (um) resting length for collagen
  const double PCon_c  = 0.02;       // (norm Force) passive force due to collagen
  const double PExp_c  = 70;         //   these apply to trabeculae and single cells only

  // Calculation of complete muscle response
  const double massf   = 0.00005e6;  // ([norm Force ms^2]/um) muscle mass
  const double visc    = 0.003e3;    // ([norm Force ms]/um) muscle viscosity
  const double KSE     = 1;          // (norm Force/um) series elastic element
  const double kxb     = 120;        // (mN/mm^2) maximal force
  const double Trop_conc = 70;       // (uM) troponin concentration

  // RICE Time-Varying Parameters
  const double TempCelsius = 37.0;
  const double contrflag   = 1;  
  const double SLset       = 2.2;
  const double SEon        = 1.0;

  double dXBprer, dXBpostr, dN_NoXB, dP_NoXB, dN, drice_P, dSL;
  double dxXBpostr, dxXBprer, dTRPNCaL, dTRPNCaH, dintf;
  double JTropRice = 0.0;

  // ---------------------------------------------------------------------------
  // calculations
  // ---------------------------------------------------------------------------
  // ORd variables
  double v      = statevars[ 0];  // v
  double nai    = statevars[ 1];  // nai
  double nass   = statevars[ 2];  // nass
  double ki     = statevars[ 3];  // ki
  double kss    = statevars[ 4];  // kss
  double cai    = statevars[ 5];  // cai
  double cass   = statevars[ 6];  // cass
  double cansr  = statevars[ 7];  // cansr
  double cajsr  = statevars[ 8];  // cajsr
  double m      = statevars[ 9];  // m
  double hf     = statevars[10];  // hf
  double hs     = statevars[11];  // hs
  double j      = statevars[12];  // j
  double hsp    = statevars[13];  // hsp;
  double jp     = statevars[14];  // jp
  double mL     = statevars[15];  // mL
  double hL     = statevars[16];  // hL
  double hLp    = statevars[17];  // hLp
  double a      = statevars[18];  // a
  double iF     = statevars[19];  // iF
  double iS     = statevars[20];  // iS
  double ap     = statevars[21];  // ap
  double iFp    = statevars[22];  // iFp
  double iSp    = statevars[23];  // iSp
  double d      = statevars[24];  // d
  double ff     = statevars[25];  // ff
  double fs     = statevars[26];  // fs
  double fcaf   = statevars[27];  // fcaf
  double fcas   = statevars[28];  // fcas
  double jca    = statevars[29];  // jca
  double nca    = statevars[30];  // nca
  double ffp    = statevars[31];  // ffp
  double fcafp  = statevars[32];  // fcafp
  double xrf    = statevars[33];  // xrf
  double xrs    = statevars[34];  // xrs
  double xs1    = statevars[35];  // xs1
  double xs2    = statevars[36];  // xs2
  double xk1    = statevars[37];  // xk1
  double Jrelnp = statevars[38];  // Jrelnp
  double Jrelp  = statevars[39];  // Jrelp
  double CaMKt  = statevars[40];  // CaMKt  
  
  // Rice variables
  double TRPNCaL  = statevars[41];  // TRPNCaL - Rice myofilament variables
  double TRPNCaH  = statevars[42];  // TRPNCaH
  double N_NoXB   = statevars[43];  // N_NoXB
  double P_NoXB   = statevars[44];  // P_NoXB  
  double N        = statevars[45];  // N
  double XBprer   = statevars[46];  // XBprer
  double XBpostr  = statevars[47];  // XBpostr
  double xXBprer  = statevars[48];  // xXBprer
  double xXBpostr = statevars[49];  // xXBpostr
  double SL       = statevars[50];  // SL
  double intf     = statevars[51];  // intf  
    
  //revpots();------------------------------------------------------------------
  ENa = (R*T/F)*log(nao/nai);
  EK  = (R*T/F)*log(ko/ki);
  EKs = (R*T/F)*log((ko+0.01833*nao)/(ki+0.01833*nai));  
  
  // d_dt_variables ------------------------------------------------------------
  double d_dt_v,   d_dt_nai,  d_dt_nass,  d_dt_ki,     d_dt_kss;
  double d_dt_cai, d_dt_cass, d_dt_cansr, d_dt_cajsr,  d_dt_m;
  double d_dt_hf,  d_dt_hs,   d_dt_j,     d_dt_hsp,    d_dt_jp;
  double d_dt_mL,  d_dt_hL,   d_dt_hLp,   d_dt_a,      d_dt_iF;
  double d_dt_iS,  d_dt_ap,   d_dt_iFp,   d_dt_iSp,    d_dt_d;
  double d_dt_ff,  d_dt_fs,   d_dt_fcaf,  d_dt_fcas,   d_dt_jca;
  double d_dt_nca, d_dt_ffp,  d_dt_fcafp, d_dt_xrf,    d_dt_xrs;
  double d_dt_xs1, d_dt_xs2,  d_dt_xk1,   d_dt_Jrelnp, d_dt_Jrelp;
  double d_dt_CaMKt;
  
  //RGC(); ---------------------------------------------------------------------
  CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
  CaMKa=CaMKb+CaMKt;
  double vffrt=v*F*F/(R*T);
  double vfrt=v*F/(R*T);

  double mss=1.0/(1.0+exp((-(v+39.57))/9.871));
  double tm=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
  d_dt_m=mss-(mss-m)*exp(-dt/tm);
  m=mss-(mss-m)*exp(-dt/tm);
  
  double hss=1.0/(1+exp((v+82.90)/6.086));
  double thf=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
  double ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
  double Ahf=0.99;
  double Ahs=1.0-Ahf;
  d_dt_hf=hss-(hss-hf)*exp(-dt/thf);
  d_dt_hs=hss-(hss-hs)*exp(-dt/ths);
  hf=hss-(hss-hf)*exp(-dt/thf);
  hs=hss-(hss-hs)*exp(-dt/ths);
    
  double h=Ahf*hf+Ahs*hs;
  double jss=hss;
  double tj=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
  d_dt_j=jss-(jss-j)*exp(-dt/tj);
  j=jss-(jss-j)*exp(-dt/tj);
  
  double hssp=1.0/(1+exp((v+89.1)/6.086));
  double thsp=3.0*ths;  
  d_dt_hsp=hssp-(hssp-hsp)*exp(-dt/thsp);
  hsp=hssp-(hssp-hsp)*exp(-dt/thsp);
  
  double hp=Ahf*hf+Ahs*hsp;
  double tjp=1.46*tj;
  d_dt_jp=jss-(jss-jp)*exp(-dt/tjp);
  jp=jss-(jss-jp)*exp(-dt/tjp);
  
  double GNa=75;
  double fINap=(1.0/(1.0+KmCaMK/CaMKa));
  INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);

  double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
  double tmL=tm;
  d_dt_mL=mLss-(mLss-mL)*exp(-dt/tmL);
  mL=mLss-(mLss-mL)*exp(-dt/tmL);
  
  double hLss=1.0/(1.0+exp((v+87.61)/7.488));
  double thL=200.0;
  d_dt_hL=hLss-(hLss-hL)*exp(-dt/thL);
  hL=hLss-(hLss-hL)*exp(-dt/thL);
  
  double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
  double thLp=3.0*thL;
  d_dt_hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
  hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
  
  double GNaL=0.0075;
  if (type==EPI)
  {
    GNaL*=0.6;
  }
#ifdef HF
  GNaL *= 6.7;
#endif
  double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
  INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

  double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
  double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
  d_dt_a=ass-(ass-a)*exp(-dt/ta);
  a=ass-(ass-a)*exp(-dt/ta);
  
  double iss=1.0/(1.0+exp((v+43.94)/5.711));
  double delta_epi;
  if (type==EPI)
  {
    delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
  }
  else
  {
    delta_epi=1.0;
  }
  double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
  double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
  tiF*=delta_epi;
  tiS*=delta_epi;
  double AiF=1.0/(1.0+exp((v-213.6)/151.2));
  double AiS=1.0-AiF;
  d_dt_iF=iss-(iss-iF)*exp(-dt/tiF);
  d_dt_iS=iss-(iss-iS)*exp(-dt/tiS);
  iF=iss-(iss-iF)*exp(-dt/tiF);
  iS=iss-(iss-iS)*exp(-dt/tiS);
    
  double i=AiF*iF+AiS*iS;
  double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
  d_dt_ap=assp-(assp-ap)*exp(-dt/ta);
  ap=assp-(assp-ap)*exp(-dt/ta);
  
  double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
  double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
  double tiFp=dti_develop*dti_recover*tiF;
  double tiSp=dti_develop*dti_recover*tiS;
  d_dt_iFp=iss-(iss-iFp)*exp(-dt/tiFp);
  d_dt_iSp=iss-(iss-iSp)*exp(-dt/tiSp);
  iFp=iss-(iss-iFp)*exp(-dt/tiFp);
  iSp=iss-(iss-iSp)*exp(-dt/tiSp);
  
  double ip=AiF*iFp+AiS*iSp;
  double Gto=0.02;
  if (type==EPI)
  {
    Gto*=4.0;
  }
  if (type==MCELL)
  {
    Gto*=4.0;
  }
#ifdef HF
  Gto *= 0.16;
//  if(type==ENDO) Gto *= 1.0;
//  else if(type==MCELL) Gto *= 0.7;
//  else if(type==EPI) Gto *= 0.7;
#endif  
  double fItop=(1.0/(1.0+KmCaMK/CaMKa));
  Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

  double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
  double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
  d_dt_d=dss-(dss-d)*exp(-dt/td);
  d=dss-(dss-d)*exp(-dt/td);
  
  double fss=1.0/(1.0+exp((v+19.58)/3.696));
  double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
  double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
  double Aff=0.6;
  double Afs=1.0-Aff;
  d_dt_ff=fss-(fss-ff)*exp(-dt/tff);
  d_dt_fs=fss-(fss-fs)*exp(-dt/tfs);
  ff=fss-(fss-ff)*exp(-dt/tff);
  fs=fss-(fss-fs)*exp(-dt/tfs);
  
  double f=Aff*ff+Afs*fs;
  double fcass=fss;
  double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
  double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
  double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
  double Afcas=1.0-Afcaf;
  d_dt_fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
  d_dt_fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
  fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
  fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
  
  double fca=Afcaf*fcaf+Afcas*fcas;
  double tjca=75.0;
  d_dt_jca=fcass-(fcass-jca)*exp(-dt/tjca);
  jca=fcass-(fcass-jca)*exp(-dt/tjca);
  
  double tffp=2.5*tff;
  d_dt_ffp=fss-(fss-ffp)*exp(-dt/tffp);
  ffp=fss-(fss-ffp)*exp(-dt/tffp);
  
  double fp=Aff*ffp+Afs*fs;
  double tfcafp=2.5*tfcaf;
  d_dt_fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
  fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
  
  double fcap=Afcaf*fcafp+Afcas*fcas;
  double Kmn=0.002;
  double k2n=1000.0;
  double km2n=jca*1.0;
  double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));  
  //nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
  //d_dt_nca = anca*k2n - nca*km2n;
  d_dt_nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
  nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
  
  double PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
  double PhiCaNa=1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
  double PhiCaK=1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
  double zca=2.0;
  double PCa=0.0001;
  if (type==EPI)
  {
    PCa*=1.2;
  }
  if (type==MCELL)
  {
    PCa*=2.5;
  }
  double PCap=1.1*PCa;
  double PCaNa=0.00125*PCa;
  double PCaK=3.574e-4*PCa;
  double PCaNap=0.00125*PCap;
  double PCaKp=3.574e-4*PCap;
  double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
  ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
  ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
  ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);

  double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
  double txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
  double txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
  double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
  double Axrs=1.0-Axrf;
  d_dt_xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
  d_dt_xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
  xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
  xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
  
  double xr=Axrf*xrf+Axrs*xrs;
  double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
  double GKr=0.046;
  if (type==EPI)
  {
    GKr*=1.3;
  }
  if (type==MCELL)
  {
    GKr*=0.8;
  }
  IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);

  double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
  double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
  d_dt_xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
  xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
  
  double xs2ss=xs1ss;
  double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
  d_dt_xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
  xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
  
  double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
  double GKs=0.0034;
  if (type==EPI)
  {
    GKs*=1.4;
  }
#ifdef HF
  //GKs *= 0.51;
#endif  
  IKs=GKs*KsCa*xs1*xs2*(v-EKs);

  double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
  double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
  d_dt_xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
  xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
  
  double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
  double GK1=0.1908;
  if (type==EPI)
  {
    GK1*=1.2;
  }
  if (type==MCELL)
  {
    GK1*=1.3;
  }
#ifdef HF
  //GK1 *= 0.65;
#endif    
  IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);

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
  double allo=1.0/(1.0+pow(KmCaAct/cai,2.0));
  double zna=1.0;
  double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
  double JncxCa=E2*k2-E1*k1;
  double Gncx=0.0008;
  if (type==EPI)
  {
    Gncx*=1.1;
  }
  if (type==MCELL)
  {
    Gncx*=1.4;
  }
#ifdef HF
  if(type==ENDO)
    Gncx = Gncx * 1.6;
  else if(type==EPI)
    Gncx = Gncx * 2.0;
  else if(type==MCELL)
    Gncx = Gncx * 1.6;
#endif 
  INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

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
  KmCaAct=150.0e-6;
  allo=1.0/(1.0+pow(KmCaAct/cass,2.0));
  JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
  JncxCa=E2*k2-E1*k1;
  INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

  INaCa=INaCa_i+INaCa_ss;

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
  double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
  double b1=k1m*MgADP;
  double a2=k2p;
  double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
  double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
  double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
  double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
  double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
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
  double Pnak=30;
  if (type==EPI)
  {
    Pnak*=0.9;
  }
  if (type==MCELL)
  {
    Pnak*=0.7;
  }
  INaK=Pnak*(zna*JnakNa+zk*JnakK);

  double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
  double GKb=0.003;
  IKb=GKb*xkb*(v-EK);

  double PNab=3.75e-10;
  INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

  double PCab=2.5e-8;
  ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);

  double GpCa=0.0005;
  IpCa=GpCa*cai/(0.0005+cai);  
  
  // RICE calculations ---------------------------------------------------------    
  
  double rice_P = 1.0 - N - XBprer - XBpostr;

  // Compute single overlap fractions
  double sovr_ze   = min(len_thick/2,SL/2);               // z-line end
  double sovr_cle  = max(SL/2-(SL-len_thin),len_hbare/2); // centerline of end
  double len_sovr  = sovr_ze-sovr_cle;                    // single overlap length
  double SOVFThick = len_sovr*2/(len_thick-len_hbare);    // thick filament overlap frac
  double SOVFThin  = len_sovr/len_thin;                   // thin filament overlap frac

  // Compute combined Ca binding to high- (w/XB) and low- (no XB) sites
  double Tropreg = (1-SOVFThin)*TRPNCaL + SOVFThin*TRPNCaH;
  double permtot = pow((1.0 / (1.0 + pow((perm50 / Tropreg), nperm))), 0.5);
  double inprmt  = min(1/permtot, 100);

  // Adjustments for Ca activation, temperature, SL, stress and strain
  double konT    = kon*pow(Qkon,(TempCelsius-37)/10.);
  double koffLT  = koffL*pow(Qkoff,(TempCelsius-37)/10.);
  double koffHT  = koffH*pow(Qkoff,(TempCelsius-37)/10.);
  double kn_pT   = kn_p*permtot*pow(Qkn_p,(TempCelsius-37.)/10.);
  double kp_nT   = kp_n*inprmt*pow(Qkp_n,(TempCelsius-37.)/10.);
  double fappT   = fapp*pow(Qfapp,(TempCelsius-37.)/10.);
  double gapslmd = 1.0 + (1.0-SOVFThick)*gslmod;
  double gappT   = gapp*gapslmd*pow(Qgapp,(TempCelsius-37)/10.);
  double hfmd    = exp(-sign(xXBprer)*hfmdc*((xXBprer/x_0)*(xXBprer/x_0)));
  double hbmd    = exp(sign((xXBpostr-x_0))*hbmdc*(((xXBpostr-x_0)/x_0)*((xXBpostr-x_0)/x_0)));
  double hfT     = rice_hf*hfmd*pow(Qhf,(TempCelsius-37)/10.);
  double hbT     = hb*hbmd*pow(Qhb,(TempCelsius-37)/10.);
  double gxbmd   = heav(x_0-xXBpostr)*exp(sigmap*((x_0-xXBpostr)/x_0)*((x_0-xXBpostr)/x_0))
             + (1.-heav(x_0-xXBpostr)*exp(sigman*(((xXBpostr-x_0)/x_0)*(xXBpostr-x_0)/x_0)));
  double gxbT    = gxb*gxbmd*pow(Qgxb,(TempCelsius-37.)/10.);

  // Regulation and crossbridge cycling state derivatives
  double CAI = 1000.0 * cai; // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  dTRPNCaL  = konT*CAI*(1.0-TRPNCaL) - koffLT*TRPNCaL;
  dTRPNCaH  = konT*CAI*(1.0-TRPNCaH) - koffHT*TRPNCaH;
  dN_NoXB   = -kn_pT*N_NoXB + kp_nT*P_NoXB;
  dP_NoXB   = -kp_nT*P_NoXB + kn_pT*N_NoXB;
  dN        = -kn_pT*N + kp_nT*rice_P;
  // dP      = -kp_nT*rice_P + kn_pT*N - fappT*rice_P + gappT*XBprer + gxbT*XBpostr;
  dXBprer   = fappT*rice_P - gappT*XBprer - hfT*XBprer + hbT*XBpostr;
  dXBpostr  = hfT*XBprer - hbT*XBpostr - gxbT*XBpostr;
  drice_P   = -(dN+dXBprer+dXBpostr);

  // steady-state fractions in XBprer and XBpostr using King-Altman rule
  double SSXBprer = (hb*fapp+gxb*fapp)/(gxb*rice_hf+fapp*rice_hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
  double SSXBpostr = fapp*rice_hf/(gxb*rice_hf+fapp*rice_hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);

  // normalization for scaling active and passive force (maximal force)
  double Fnordv = kxb*x_0*SSXBpostr;

  // Calculate Forces (active, passive, preload, afterload)
  double force = kxb*SOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer);
  active_force = force/Fnordv;
  
  double ppforce_t = sign(SL-SLrest)*PCon_t*(exp(PExp_t*abs(SL-SLrest))-1);
  double ppforce_c = heav(SL-SL_c)*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1);
  // ppforce_c = 0;
  double ppforce = ppforce_t + ppforce_c;
  double preload = sign(SLset-SLrest)*PCon_t*(exp(PExp_t*abs(SLset-SLrest))-1);
  double afterload = 0;  // either constant or due to series elastic element
  if (SEon == 1) {
    afterload = KSE*(SLset-SL);
  }

  dintf = (-ppforce+preload-active_force+afterload);
  dSL = contrflag*((intf+(SLset-SL)*visc)/massf)*heav(SL-SLmin)*heav(SLmax-SL);

  // Mean strain of strongly-bound states due to SL motion and XB cycling
  double dutyprer  = (hbT*fappT+gxbT*fappT)/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
  double dutypostr = fappT*hfT/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
  dxXBprer  = dSL/2.0 + xPsi*(1.0/dutyprer)*(-xXBprer*fappT+(xXBpostr-x_0-xXBprer)*hbT);
  dxXBpostr = dSL/2.0 + xPsi*(1.0/dutypostr)*(x_0+xXBprer-xXBpostr)*hfT;

  // Ca buffering by low-affinity troponin C (LTRPNCa)
  double FrSBXB  = (XBpostr+XBprer)/(SSXBpostr + SSXBprer);
  double dFrSBXB = (dXBpostr+dXBprer)/(SSXBpostr + SSXBprer);

  double dsovr_ze   = -(dSL/2.0)*heav(len_thick-SL);
  double dsovr_cle  = -(dSL/2.0)*heav((2.0*len_thin-SL)-len_hbare);
  double dlen_sovr  = dsovr_ze-dsovr_cle;
  double dSOVFThin  = dlen_sovr/len_thin;
  double dSOVFThick = 2.0*dlen_sovr/(len_thick-len_hbare);

  double TropTot  = Trop_conc*((1-SOVFThin)*TRPNCaL + SOVFThin*(FrSBXB*TRPNCaH+(1-FrSBXB)*TRPNCaL));
  double dTropTot = Trop_conc*(-dSOVFThin*TRPNCaL+(1-SOVFThin)*dTRPNCaL +
                               dSOVFThin*(FrSBXB*TRPNCaH+(1-FrSBXB)*TRPNCaL) +
                               SOVFThin*(dFrSBXB*TRPNCaH+FrSBXB*dTRPNCaH-dFrSBXB*TRPNCaL +
                               (1-FrSBXB)*dTRPNCaL));

  JTropRice = dTropTot;

  double dforce = kxb*dSOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer) +
                  kxb*SOVFThick*(dxXBpostr*XBpostr+xXBpostr*dXBpostr +
                  dxXBprer*XBprer+xXBprer*dXBprer);

  double dactive = 0.5*dforce/Fnordv;
  double dppforce_t = sign(SL-SLrest)*PCon_t*PExp_t*dSL*exp(PExp_t*abs(SL-SLrest));
  double dppforce_c = heav(SL-SL_c)*PCon_c*PExp_c*dSL*exp(PExp_c*abs(SL-SL_c));
  // ppforce_c = 0;
  double dppforce = dppforce_t + dppforce_c;
  double dsfib = dppforce+dactive;
  
  //stimulus();-----------------------------------------------------------------
  //const double CL=1000;      // pacing cycle length
  //const double amp=-80;      // stimulus amplitude in uA/uF
  //const double start=0;      // start time of the stimulus, relative to each beat
  //const double duration=0.5; // duration of teh stimulus in ms
  //const double end=1000 * CL;
  
  //if(((time>=start)&&(time<=end)&&(((time-start)-(floor(((time-start)/CL))*CL))<=duration)))
  //  Ist = amp;
  //else
  //  Ist = 0;  
  
  Ist = i_stim;
      
  //voltage();------------------------------------------------------------------
  //v+=-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist);
  d_dt_v = -(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist);  
  
  //FBC();----------------------------------------------------------------------
  double CaMKb_new=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
  
  CaMKa=CaMKb_new+CaMKt;
  d_dt_CaMKt = (aCaMK*CaMKb_new*(CaMKb_new+CaMKt)-bCaMK*CaMKt);
  CaMKt+=dt*(aCaMK*CaMKb_new*(CaMKb_new+CaMKt)-bCaMK*CaMKt);
  
  JdiffNa=(nass-nai)/2.0;
  JdiffK=(kss-ki)/2.0;
  Jdiff=(cass-cai)/0.2;

  double bt=4.75;
  double a_rel=0.5*bt;
  double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
  if (type==MCELL)
  {
    Jrel_inf*=1.7;
  }
  double tau_rel=bt/(1.0+0.0123/cajsr);
  if (tau_rel<0.005)
  {
    tau_rel=0.005;
  }
  d_dt_Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
  Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
  
  double btp=1.25*bt;
  double a_relp=0.5*btp;
  double Jrel_infp=a_relp*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
  if (type==MCELL)
  {
    Jrel_infp*=1.7;
  }
  double tau_relp=btp/(1.0+0.0123/cajsr);
  if (tau_relp<0.005)
  {
    tau_relp=0.005;
  }
  d_dt_Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
  Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
  
  double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
  Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;

  double Jupnp=0.004375*cai/(cai+0.00092);
  double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
  if (type==EPI)
  {
    Jupnp*=1.3;
    Jupp*=1.3;
  }
#ifdef HF

#endif  
  double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
  Jleak=0.0039375*cansr/15.0;
  Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
  Jtr=(cansr-cajsr)/100.0;

  d_dt_nai  = (-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
  d_dt_nass = (-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

  d_dt_ki  = (-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
  d_dt_kss = (-(ICaK)*Acap/(F*vss)-JdiffK);

  double Bcai;
  if (type==EPI)
  {
    //Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0));
  }
  else
  {
    //Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0));
  }  
  //cai+=dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo - (JTropRice/1000.0) ));  
  //d_dt_cai=(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));
  d_dt_cai=(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo- (JTropRice/1000.0) ));

  double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
  d_dt_cass=(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));
  d_dt_cansr=(Jup-Jtr*vjsr/vnsr);

  double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0));
  d_dt_cajsr=(Bcajsr*(Jtr-Jrel));
  
  // copy new values -----------------------------------------------------------
  
  values[ 0] = d_dt_v;
  values[ 1] = d_dt_nai;
  values[ 2] = d_dt_nass;
  values[ 3] = d_dt_ki;
  values[ 4] = d_dt_kss;
  values[ 5] = d_dt_cai;
  values[ 6] = d_dt_cass;
  values[ 7] = d_dt_cansr;
  values[ 8] = d_dt_cajsr;
  values[ 9] = d_dt_m;
  values[10] = d_dt_hf;
  values[11] = d_dt_hs;
  values[12] = d_dt_j;
  values[13] = d_dt_hsp;
  values[14] = d_dt_jp;
  values[15] = d_dt_mL;
  values[16] = d_dt_hL;
  values[17] = d_dt_hLp;
  values[18] = d_dt_a;
  values[19] = d_dt_iF;
  values[20] = d_dt_iS;
  values[21] = d_dt_ap;
  values[22] = d_dt_iFp;
  values[23] = d_dt_iSp;
  values[24] = d_dt_d;
  values[25] = d_dt_ff;
  values[26] = d_dt_fs;
  values[27] = d_dt_fcaf;
  values[28] = d_dt_fcas;
  values[29] = d_dt_jca;
  values[30] = d_dt_nca;
  values[31] = d_dt_ffp;
  values[32] = d_dt_fcafp;
  values[33] = d_dt_xrf;
  values[34] = d_dt_xrs;
  values[35] = d_dt_xs1;
  values[36] = d_dt_xs2;
  values[37] = d_dt_xk1;  
  values[38] = d_dt_Jrelnp;
  values[39] = d_dt_Jrelp;
  values[40] = d_dt_CaMKt;
    
  // Rice update
  values[41] = dTRPNCaL;
  values[42] = dTRPNCaH;
  values[43] = dN_NoXB;
  values[44] = dP_NoXB;
  values[45] = dN;
  values[46] = dXBprer;
  values[47] = dXBpostr;
  values[48] = dxXBprer;
  values[49] = dxXBpostr;
  values[50] = dSL;
  values[51] = dintf;
}
