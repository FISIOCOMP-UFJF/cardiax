#include "iyer.h"

// Setting initial conditions to the Iyer 2004 model
void setIC_ode_cpu(Real *sv, int cellID, int NEQ = 67)
{
   // Vector containing the initial conditions of the cell
   Real IC[NEQ];

   IC[ 0] = -9.065755929e+01;     // V       
   IC[ 1] = 1.2555894320e+02;     // Ki      
   IC[ 2] = 8.6011920159e-05;     // Cai     
   IC[ 3] = 1.4202152450e-04;     // CaSS    
   IC[ 4] = 2.8522394460e-01;     // CaJSR   
   IC[ 5] = 2.8552949150e-01;     // CaNSR   
   IC[ 6] = 9.7983041621e+00;     // Nai     
   IC[ 7] = 8.0465849730e-02;     // LTRPNCa 
   IC[ 8] = 9.7721524810e-01;     // HTRPNCa 
   IC[ 9] = 4.4480389460e-01;     // C1_RyR  
   IC[10] = 3.9739200600e-09;     // O2_RyR  
   IC[11] = 5.5453592400e-01;     // C2_RyR  
   IC[12] = 6.6017832870e-04;     // O1_RyR  
   IC[13] = 8.6219340540e-01;     // C0      
   IC[14] = 1.1413982110e-02;     // C1      
   IC[15] = 5.6663252250e-05;     // C2      
   IC[16] = 1.2502135700e-07;     // C3      
   IC[17] = 1.0344267480e-10;     // C4      
   IC[18] = 7.7579815630e-12;     // Open    
   IC[19] = 1.1985470810e-01;     // CCa0    
   IC[20] = 6.3467943020e-03;     // CCa1    
   IC[21] = 1.2603264880e-04;     // CCa2    
   IC[22] = 1.1123152380e-06;     // CCa3    
   IC[23] = 3.6813158920e-09;     // CCa4    
   IC[24] = 9.9971570740e-01;     // yCa     
   IC[25] = 9.5137213510e-01;     // C0Kv43  
   IC[26] = 2.6682880890e-02;     // C1Kv43  
   IC[27] = 2.8063803580e-04;     // C2Kv43  
   IC[28] = 1.3118375790e-06;     // C3Kv43  
   IC[29] = 2.2995562400e-09;     // OKv43   
   IC[30] = 1.5130252040e-02;     // CI0Kv43 
   IC[31] = 5.4429646010e-03;     // CI1Kv43 
   IC[32] = 9.9183733591e-04;     // CI2Kv43 
   IC[33] = 9.5143860571e-05;     // CI3Kv43 
   IC[34] = 2.7426773820e-06;     // OIKv43  
   IC[35] = 5.9770997650e-01;     // C0Kv14  
   IC[36] = 1.7309905280e-01;     // C1Kv14  
   IC[37] = 1.8810723860e-02;     // C2Kv14  
   IC[38] = 9.1607013500e-04;     // C3Kv14  
   IC[39] = 1.9755413570e-05;     // OKv14   
   IC[40] = 3.5390843460e-02;     // CI0Kv14 
   IC[41] = 5.4288243530e-03;     // CI1Kv14 
   IC[42] = 2.2878588690e-02;     // CI2Kv14 
   IC[43] = 3.2338000030e-02;     // CI3Kv14 
   IC[44] = 1.1340820580e-01;     // OIKv14  
   IC[45] = 1.4375756490e-01;     // na1     
   IC[46] = 4.1777620800e-02;     // na2     
   IC[47] = 4.5528983640e-03;     // na3     
   IC[48] = 2.2052074300e-04;     // na4     
   IC[49] = 4.0052664840e-06;     // na5     
   IC[50] = 1.5744274900e-08;     // na6     
   IC[51] = 2.8566071790e-09;     // na7     
   IC[52] = 4.7499500080e-01;     // na8     
   IC[53] = 2.7072140970e-01;     // na9     
   IC[54] = 5.7861200570e-02;     // na10    
   IC[55] = 5.4962922790e-03;     // na11    
   IC[56] = 1.9578745190e-04;     // na12    
   IC[57] = 4.1767939600e-04;     // na13    
   IC[58] = 9.9669733800e-01;     // C1Herg  
   IC[59] = 4.3408796480e-04;     // C2Herg  
   IC[60] = 7.6340997551e-05;     // C3Herg  
   IC[61] = 9.5117891130e-06;     // OHerg   
   IC[62] = 1.5333470070e-06;     // IHerg   
   IC[63] = 9.6456062950e-01;     // C0ks    
   IC[64] = 3.5426135680e-02;     // C1ks    
   IC[65] = 2.4917106960e-07;     // O1ks    
   IC[66] = 1.2985478220e-05;     // O2ks    
                                       
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
                                       
   RHS_Iyer_cpu(time, rY, rDY, col);
                                       
   for(int i = 0; i < NEQ; i++)        
      sv[cellID+i] = dt*rDY[i] + rY[i];
}                                      
                                       
// The Iyer 2004 model                 
void RHS_Iyer_cpu(Real time, Real *rY, Real *rDY, int col)
{
   Real V       = rY[ 0];
   Real Ki      = rY[ 1];
   Real Cai     = rY[ 2];
   Real CaSS    = rY[ 3];
   Real CaJSR   = rY[ 4];
   Real CaNSR   = rY[ 5];
   Real Nai     = rY[ 6];
   Real LTRPNCa = rY[ 7];
   Real HTRPNCa = rY[ 8];
   Real C1_RyR  = rY[ 9];
   Real O2_RyR  = rY[10];
   Real C2_RyR  = rY[11];
   Real O1_RyR  = rY[12];
   Real C0      = rY[13];
   Real C1      = rY[14];
   Real C2      = rY[15];
   Real C3      = rY[16];
   Real C4      = rY[17];
   Real Open    = rY[18];
   Real CCa0    = rY[19];
   Real CCa1    = rY[20];
   Real CCa2    = rY[21];
   Real CCa3    = rY[22];
   Real CCa4    = rY[23];
   Real yCa     = rY[24];
   Real C0Kv43  = rY[25];
   Real C1Kv43  = rY[26];
   Real C2Kv43  = rY[27];
   Real C3Kv43  = rY[28];
   Real OKv43   = rY[29];
   Real CI0Kv43 = rY[30];
   Real CI1Kv43 = rY[31];
   Real CI2Kv43 = rY[32];
   Real CI3Kv43 = rY[33];
   Real OIKv43  = rY[34];
   Real C0Kv14  = rY[35];
   Real C1Kv14  = rY[36];
   Real C2Kv14  = rY[37];
   Real C3Kv14  = rY[38];
   Real OKv14   = rY[39];
   Real CI0Kv14 = rY[40];
   Real CI1Kv14 = rY[41];
   Real CI2Kv14 = rY[42];
   Real CI3Kv14 = rY[43];
   Real OIKv14  = rY[44];
   Real na1     = rY[45];
   Real na2     = rY[46];
   Real na3     = rY[47];
   Real na4     = rY[48];
   Real na5     = rY[49];
   Real na6     = rY[50];
   Real na7     = rY[51];
   Real na8     = rY[52];
   Real na9     = rY[53];
   Real na10    = rY[54];
   Real na11    = rY[55];
   Real na12    = rY[56];
   Real na13    = rY[57];
   Real C1Herg  = rY[58];
   Real C2Herg  = rY[59];
   Real C3Herg  = rY[60];
   Real OHerg   = rY[61];
   Real IHerg   = rY[62];
   Real C0ks    = rY[63];
   Real C1ks    = rY[64];
   Real O1ks    = rY[65];
   Real O2ks    = rY[66];

   // PARAMETERS
   Real Rgas;    // joule_per_mole_kelvin
   Real Temp;    // kelvin
   Real Faraday;    // coulomb_per_millimole
   Real Acap;    // cm2
   Real Vmyo;    // uL
   Real VSS;     // uL
   Real VJSR;    // uL
   Real VNSR;    // uL
   Real stim_start;     // ms
   Real stim_duration;  // ms
   Real stim_period;    // ms
   Real stim_end;       // ms
   Real stim_amplitude; // uA_per_uF
   Real Kfb;     // mM
   Real Nfb;     // dimensionless
   Real Krb;     // mM
   Real Nrb;     // dimensionless
   Real KSR;     // mM
   Real vmaxf;   // per_ms
   Real vmaxr;   // per_ms
   Real v1;   // per_ms
   Real tautr;   // ms
   Real tauxfer;    // ms
   Real kltrpn_plus;   // per_mM_per_ms
   Real kltrpn_minus;     // per_ms
   Real khtrpn_plus;   // per_mM_per_ms
   Real khtrpn_minus;     // per_ms
   Real LTRPNtot;   // mM
   Real HTRPNtot;   // mM
   Real CMDNtot;    // mM
   Real KmCMDN;     // mM
   Real EGTAtot;    // mM
   Real KmEGTA;     // mM
   Real CSQNtot;    // mM
   Real KmCSQN;     // mM
   Real kaplus;     // per_uM4_per_ms
   int ncoop;   // dimensionless
   Real kaminus;    // per_ms
   Real kbplus;     // per_uM3_per_ms
   int mcoop;   // dimensionless
   Real kbminus;    // per_ms
   Real kcplus;     // per_ms
   Real kcminus;    // per_ms
   Real aL;   // dimensionless
   Real bL;   // dimensionless
   Real fL;   // per_ms
   Real gL;   // per_ms
   Real alphaa0Kv43;   // per_ms
   Real aaKv43;     // per_mV
   Real betaa0Kv43;    // per_ms
   Real baKv43;     // per_mV
   Real alphai0Kv43;   // per_ms
   Real aiKv43;     // per_mV
   Real betai0Kv43;    // per_ms
   Real biKv43;     // per_mV
   Real b1Kv43;     // dimensionless
   Real b2Kv43;     // dimensionless
   Real b3Kv43;     // dimensionless
   Real b4Kv43;     // dimensionless
   Real f1Kv43;     // dimensionless
   Real f2Kv43;     // dimensionless
   Real f3Kv43;     // dimensionless
   Real f4Kv43;     // dimensionless
   Real alphaa0Kv14;   // per_ms
   Real aaKv14;     // per_mV
   Real betaa0Kv14;    // per_ms
   Real baKv14;     // per_mV
   Real alphai0Kv14;   // per_ms
   Real betai0Kv14;    // per_ms
   Real b1Kv14;     // dimensionless
   Real b2Kv14;     // dimensionless
   Real b3Kv14;     // dimensionless
   Real b4Kv14;     // dimensionless
   Real f1Kv14;     // dimensionless
   Real f2Kv14;     // dimensionless
   Real f3Kv14;     // dimensionless
   Real f4Kv14;     // dimensionless
   Real Nao;     // mM
   Real Ko;   // mM
   Real Cao;     // mM
   Real Kv43Frac;   // dimensionless
   Real KvScale;    // dimensionless
   Real GNa;     // mS_per_uF
   Real GKr;     // mS_per_uF
   Real GKs;     // mS_per_uF
   Real GK1;     // mS_per_uF
   Real GNab;    // mS_per_uF
   Real GCab;    // mS_per_uF
   Real IpCamax;    // uA_per_uF
   Real KmpCa;   // mM
   Real KmKo;    // mM
   Real KmNai;   // mM
   Real INaKmax;    // uA_per_uF
   Real eta;     // dimensionless
   Real ksat;    // dimensionless
   Real KmCa;    // mM
   Real KmNa;    // mM
   Real kNaCa;   // uA_per_uF
   Real Pscale;     // dimensionless
   Real ICahalf;    // uA_per_uF
   Real RToverF;    // mV
   Real TNa;     // kelvin
   Real T_Const_HERG;     // dimensionless
   Real A0_HERG;    // per_ms
   Real B0_HERG;    // per_mV
   Real A1_HERG;    // per_ms
   Real B1_HERG;    // per_mV
   Real A2_HERG;    // per_ms
   Real B2_HERG;    // per_mV
   Real A3_HERG;    // per_ms
   Real B3_HERG;    // per_mV
   Real A4_HERG;    // per_ms
   Real B4_HERG;    // per_mV
   Real A5_HERG;    // per_ms
   Real B5_HERG;    // per_mV
   Real A6_HERG;    // per_ms
   Real B6_HERG;    // per_mV
   Real calc_RT_over_F;   // mV
   Real calc_C;     // mF
   Real calc_a1;    // mM
   Real calc_a2;    // dimensionless
   Real calc_i_tot;    // uA_per_uF
   Real calc_i_Stim;   // uA_per_uF
   Real calc_fb;    // dimensionless
   Real calc_rb;    // dimensionless
   Real calc_Jup;   // mM_per_ms
   Real calc_Jrel;     // mM_per_ms
   Real calc_Jtr;   // mM_per_ms
   Real calc_Jxfer;    // mM_per_ms
   Real calc_dLTRPNCa;    // per_ms
   Real calc_dHTRPNCa;    // per_ms
   Real calc_Jtrpn;    // mM_per_ms
   Real calc_beta_SS;     // dimensionless
   Real calc_beta_JSR;    // dimensionless
   Real calc_beta_i;   // dimensionless
   Real calc_dC1_RyR;     // per_ms
   Real calc_dO2_RyR;     // per_ms
   Real calc_dC2_RyR;     // per_ms
   Real calc_dO1_RyR;     // per_ms
   Real calc_alpha;    // per_ms
   Real calc_beta;     // per_ms
   Real calc_alpha_prime;    // per_ms
   Real calc_beta_prime;     // per_ms
   Real calc_gamma;    // per_ms
   Real calc_omega;    // per_ms
   Real calc_C0_to_C1;    // per_ms
   Real calc_C1_to_C2;    // per_ms
   Real calc_C2_to_C3;    // per_ms
   Real calc_C3_to_C4;    // per_ms
   Real calc_CCa0_to_CCa1;   // per_ms
   Real calc_CCa1_to_CCa2;   // per_ms
   Real calc_CCa2_to_CCa3;   // per_ms
   Real calc_CCa3_to_CCa4;   // per_ms
   Real calc_C1_to_C0;    // per_ms
   Real calc_C2_to_C1;    // per_ms
   Real calc_C3_to_C2;    // per_ms
   Real calc_C4_to_C3;    // per_ms
   Real calc_CCa1_to_CCa0;   // per_ms
   Real calc_CCa2_to_CCa1;   // per_ms
   Real calc_CCa3_to_CCa2;   // per_ms
   Real calc_CCa4_to_CCa3;   // per_ms
   Real calc_C0_to_CCa0;     // per_ms
   Real calc_C1_to_CCa1;     // per_ms
   Real calc_C2_to_CCa2;     // per_ms
   Real calc_C3_to_CCa3;     // per_ms
   Real calc_C4_to_CCa4;     // per_ms
   Real calc_CCa0_to_C0;     // per_ms
   Real calc_CCa1_to_C1;     // per_ms
   Real calc_CCa2_to_C2;     // per_ms
   Real calc_CCa3_to_C3;     // per_ms
   Real calc_CCa4_to_C4;     // per_ms
   Real calc_a1_C0;    // per_ms
   Real calc_a2_C0;    // per_ms
   Real calc_a1_C1;    // per_ms
   Real calc_a2_C1;    // per_ms
   Real calc_a1_C2;    // per_ms
   Real calc_a2_C2;    // per_ms
   Real calc_a1_C3;    // per_ms
   Real calc_a2_C3;    // per_ms
   Real calc_a1_C4;    // per_ms
   Real calc_a2_C4;    // per_ms
   Real calc_a1_Ca0;   // per_ms
   Real calc_a2_Ca0;   // per_ms
   Real calc_a1_Ca1;   // per_ms
   Real calc_a2_Ca1;   // per_ms
   Real calc_a1_Ca2;   // per_ms
   Real calc_a2_Ca2;   // per_ms
   Real calc_a1_Ca3;   // per_ms
   Real calc_a2_Ca3;   // per_ms
   Real calc_a1_Ca4;   // per_ms
   Real calc_a2_Ca4;   // per_ms
   Real calc_a1_Cainf;    // dimensionless
   Real calc_yCa_inf;     // dimensionless
   Real calc_tau_yCa;     // ms
   Real calc_alpha_act43;    // per_ms
   Real calc_beta_act43;     // per_ms
   Real calc_alpha_inact43;     // per_ms
   Real calc_beta_inact43;   // per_ms
   Real calc_C0Kv43_to_C1Kv43;     // per_ms
   Real calc_C1Kv43_to_C2Kv43;     // per_ms
   Real calc_C2Kv43_to_C3Kv43;     // per_ms
   Real calc_C3Kv43_to_OKv43;   // per_ms
   Real calc_CI0Kv43_to_CI1Kv43;   // per_ms
   Real calc_CI1Kv43_to_CI2Kv43;   // per_ms
   Real calc_CI2Kv43_to_CI3Kv43;   // per_ms
   Real calc_CI3Kv43_to_OIKv43;    // per_ms
   Real calc_C1Kv43_to_C0Kv43;     // per_ms
   Real calc_C2Kv43_to_C1Kv43;     // per_ms
   Real calc_C3Kv43_to_C2Kv43;     // per_ms
   Real calc_OKv43_to_C3Kv43;   // per_ms
   Real calc_CI1Kv43_to_CI0Kv43;   // per_ms
   Real calc_CI2Kv43_to_CI1Kv43;   // per_ms
   Real calc_CI3Kv43_to_CI2Kv43;   // per_ms
   Real calc_OIKv43_to_CI3Kv43;    // per_ms
   Real calc_C0Kv43_to_CI0Kv43;    // per_ms
   Real calc_C1Kv43_to_CI1Kv43;    // per_ms
   Real calc_C2Kv43_to_CI2Kv43;    // per_ms
   Real calc_C3Kv43_to_CI3Kv43;    // per_ms
   Real calc_OKv43_to_OIKv43;   // per_ms
   Real calc_CI0Kv43_to_C0Kv43;    // per_ms
   Real calc_CI1Kv43_to_C1Kv43;    // per_ms
   Real calc_CI2Kv43_to_C2Kv43;    // per_ms
   Real calc_CI3Kv43_to_C3Kv43;    // per_ms
   Real calc_OIKv43_to_OKv43;   // per_ms
   Real calc_a1_C043;     // per_ms
   Real calc_a2_C043;     // per_ms
   Real calc_a1_C143;     // per_ms
   Real calc_a2_C143;     // per_ms
   Real calc_a1_C243;     // per_ms
   Real calc_a2_C243;     // per_ms
   Real calc_a1_C343;     // per_ms
   Real calc_a2_C343;     // per_ms
   Real calc_a1_O43;   // per_ms
   Real calc_a2_O43;   // per_ms
   Real calc_a1_I043;     // per_ms
   Real calc_a2_I043;     // per_ms
   Real calc_a1_I143;     // per_ms
   Real calc_a2_I143;     // per_ms
   Real calc_a1_I243;     // per_ms
   Real calc_a2_I243;     // per_ms
   Real calc_a1_I343;     // per_ms
   Real calc_a2_I343;     // per_ms
   Real calc_a1_OI43;     // per_ms
   Real calc_a2_OI43;     // per_ms
   Real calc_alpha_act14;    // per_ms
   Real calc_beta_act14;     // per_ms
   Real calc_alpha_inact14;     // per_ms
   Real calc_beta_inact14;   // per_ms
   Real calc_C0Kv14_to_C1Kv14;     // per_ms
   Real calc_C1Kv14_to_C2Kv14;     // per_ms
   Real calc_C2Kv14_to_C3Kv14;     // per_ms
   Real calc_C3Kv14_to_OKv14;   // per_ms
   Real calc_CI0Kv14_to_CI1Kv14;   // per_ms
   Real calc_CI1Kv14_to_CI2Kv14;   // per_ms
   Real calc_CI2Kv14_to_CI3Kv14;   // per_ms
   Real calc_CI3Kv14_to_OIKv14;    // per_ms
   Real calc_C1Kv14_to_C0Kv14;     // per_ms
   Real calc_C2Kv14_to_C1Kv14;     // per_ms
   Real calc_C3Kv14_to_C2Kv14;     // per_ms
   Real calc_OKv14_to_C3Kv14;   // per_ms
   Real calc_CI1Kv14_to_CI0Kv14;   // per_ms
   Real calc_CI2Kv14_to_CI1Kv14;   // per_ms
   Real calc_CI3Kv14_to_CI2Kv14;   // per_ms
   Real calc_OIKv14_to_CI3Kv14;    // per_ms
   Real calc_C0Kv14_to_CI0Kv14;    // per_ms
   Real calc_C1Kv14_to_CI1Kv14;    // per_ms
   Real calc_C2Kv14_to_CI2Kv14;    // per_ms
   Real calc_C3Kv14_to_CI3Kv14;    // per_ms
   Real calc_OKv14_to_OIKv14;   // per_ms
   Real calc_CI0Kv14_to_C0Kv14;    // per_ms
   Real calc_CI1Kv14_to_C1Kv14;    // per_ms
   Real calc_CI2Kv14_to_C2Kv14;    // per_ms
   Real calc_CI3Kv14_to_C3Kv14;    // per_ms
   Real calc_OIKv14_to_OKv14;   // per_ms
   Real calc_a1_O;     // per_ms
   Real calc_a2_O;     // per_ms
   Real calc_a1_CI0;   // per_ms
   Real calc_a2_CI0;   // per_ms
   Real calc_a1_CI1;   // per_ms
   Real calc_a2_CI1;   // per_ms
   Real calc_a1_CI2;   // per_ms
   Real calc_a2_CI2;   // per_ms
   Real calc_a1_CI3;   // per_ms
   Real calc_a2_CI3;   // per_ms
   Real calc_a1_OI;    // per_ms
   Real calc_a2_OI;    // per_ms
   Real calc_ENa;   // mV
   Real calc_EK;    // mV
   Real calc_EKs;   // mV
   Real calc_ECa;   // mV
   Real calc_PKv14;    // litre_per_farad_second
   Real calc_GKv43;    // mS_per_uF
   Real calc_INa;   // uA_per_uF
   Real calc_fKo;   // dimensionless
   Real calc_IKr;   // uA_per_uF
   Real calc_IKs;   // uA_per_uF
   Real calc_IKv43;    // uA_per_uF
   Real calc_VF_over_RT;     // dimensionless
   Real calc_VFsq_over_RT;   // coulomb_per_millimole
   Real calc_a1_K;     // mM
   Real calc_IKv14_K;     // uA_per_uF
   Real calc_a1_Na;    // dimensionless
   Real calc_IKv14_Na;    // uA_per_uF
   Real calc_IKv14;    // uA_per_uF
   Real calc_Ito1;     // uA_per_uF
   Real calc_K1_inf;   // dimensionless
   Real calc_IK1;   // uA_per_uF
   Real calc_INab;     // uA_per_uF
   Real calc_ICab;     // uA_per_uF
   Real calc_IpCa;     // uA_per_uF
   Real calc_sigma;    // dimensionless
   Real calc_a2_Na;    // dimensionless
   Real calc_fNaK;     // dimensionless
   Real calc_a2_K;     // dimensionless
   Real calc_INaK;     // uA_per_uF
   Real calc_a1_ncx;   // mM4
   Real calc_a2_ncx;   // mM4
   Real calc_a3_ncx;   // dimensionless
   Real calc_a4_ncx;   // mM
   Real calc_a5_ncx;   // mM3
   Real calc_INaCa;    // uA_per_uF
   Real calc_PK;    // litre_per_farad_second
   Real calc_PCa;   // litre_per_farad_second
   Real calc_a1_Ca;    // mM
   Real calc_a2_Ca;    // dimensionless
   Real calc_ICamax;   // uA_per_uF
   Real calc_ICa;   // uA_per_uF
   Real calc_Icabar;   // uA_per_uF
   Real calc_PKprime;     // litre_per_farad_second
   Real calc_ICaK;     // uA_per_uF
   Real calc_FoverRT;     // per_mV
   Real calc_KToverH;     // per_ms
   Real calc_RTNa;     // joule_per_mole
   Real calc_RTNaF;    // mV
   Real calc_Temp_Scale;     // dimensionless
   Real calc_alpha1;   // per_ms
   Real calc_beta1;    // per_ms
   Real calc_gamma1;   // per_ms
   Real calc_Delta1;   // per_ms
   Real calc_On;    // per_ms
   Real calc_Of;    // per_ms
   Real calc_GammaGamma;     // per_ms
   Real calc_DeltaDelta;     // per_ms
   Real calc_epsilon;     // per_ms
   Real calc_omega_na;    // per_ms
   Real calc_rho;   // per_ms
   Real calc_mu;    // per_ms
   Real calc_Cn;    // per_ms
   Real calc_Cf;    // per_ms
   Real calc_parameter_a;    // dimensionless
   Real calc_k12;   // per_ms
   Real calc_k23;   // per_ms
   Real calc_k34;   // per_ms
   Real calc_k45;   // per_ms
   Real calc_k56;   // per_ms
   Real calc_k67;   // per_ms
   Real calc_k89;   // per_ms
   Real calc_k910;     // per_ms
   Real calc_k1011;    // per_ms
   Real calc_k1112;    // per_ms
   Real calc_k1213;    // per_ms
   Real calc_k57;   // per_ms
   Real calc_k21;   // per_ms
   Real calc_k32;   // per_ms
   Real calc_k43;   // per_ms
   Real calc_k54;   // per_ms
   Real calc_k65;   // per_ms
   Real calc_k76;   // per_ms
   Real calc_k98;   // per_ms
   Real calc_k109;     // per_ms
   Real calc_k1110;    // per_ms
   Real calc_k1211;    // per_ms
   Real calc_k1312;    // per_ms
   Real calc_k75;   // per_ms
   Real calc_k81;   // per_ms
   Real calc_k92;   // per_ms
   Real calc_k103;     // per_ms
   Real calc_k114;     // per_ms
   Real calc_k125;     // per_ms
   Real calc_k136;     // per_ms
   Real calc_k18;   // per_ms
   Real calc_k29;   // per_ms
   Real calc_k310;     // per_ms
   Real calc_k411;     // per_ms
   Real calc_k512;     // per_ms
   Real calc_k613;     // per_ms
   Real calc_C2H_to_C3H;     // per_ms
   Real calc_C3H_to_C2H;     // per_ms
   Real calc_C1H_to_C2H;     // per_ms
   Real calc_C2H_to_C1H;     // per_ms
   Real calc_C3H_to_OH;   // per_ms
   Real calc_OH_to_C3H;   // per_ms
   Real calc_OH_to_IH;    // per_ms
   Real calc_IH_to_OH;    // per_ms
   Real calc_C3H_to_IH;   // per_ms
   Real calc_IH_to_C3H;   // per_ms
   Real calc_a1_I;     // per_ms
   Real calc_a2_I;     // per_ms
   Real calc_C0ks_C1ks;   // per_ms
   Real calc_C1ks_O1ks;   // per_ms
   Real calc_O1ks_O2ks;   // per_ms
   Real calc_O1ks_C1ks;   // per_ms
   Real calc_O2ks_O1ks;   // per_ms
   Real calc_C1ks_C0ks;   // per_ms
   Real calc_a1_C0_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a2_C0_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a1_C1_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a2_C1_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a1_C2_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a2_C2_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a1_C3_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a2_C3_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES;   // (null)
   Real calc_a1_duplicated_COMPUTE_REVERSAL_POTENTIALS;    // (null)
   Real calc_a2_duplicated_COMPUTE_REVERSAL_POTENTIALS;    // (null)
   Real calc_a2_duplicated_COMPUTE_INa_IKr_IKs_Ito1_IK1_INab_IKp;   // (null)
   Real calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa;    // (null)
   Real calc_a1_Na_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa;   // (null)
   Real calc_a1_K_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa;    // (null)
   Real calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK;    // (null)
   Real calc_VFsq_over_RT_duplicated_COMPUTE_ICa_ICaK;     // (null)
   Real calc_a1_K_duplicated_COMPUTE_ICa_ICaK;    // (null)
   Real calc_a2_K_duplicated_COMPUTE_ICa_ICaK;    // (null)
   Real calc_a1_C2_duplicated_IKr;    // (null)
   Real calc_a2_C2_duplicated_IKr;    // (null)
   Real calc_a1_C3_duplicated_IKr;    // (null)
   Real calc_a2_C3_duplicated_IKr;    // (null)
   Real calc_a1_O_duplicated_IKr;     // (null)
   Real calc_a2_O_duplicated_IKr;     // (null)

   // Stimulus
   stim_amplitude =  -52.0;
//    stim_amplitude =  0.0;
   stim_start = 100.0 + (Real)col;
   stim_duration = 1.0;
   stim_period = 1000.0;
   stim_end = 100000.0;
   
   Rgas = 8.3150e+00;
   Temp = 3.10e+02;
   Faraday = 9.650e+01;
   Acap = 1.5340e-04;
   Vmyo = 2.5840e-05;
   VSS = 1.20e-09;
   VJSR = 1.60e-07;
   VNSR = 2.10e-06;
   stim_period = 1.0e+03;
   stim_amplitude = -1.50e+01;
   stim_duration = 3.0e+00;
   Kfb = 1.680e-04;
   Nfb = 1.20e+00;
   Krb = 3.290e+00;
   Nrb = 1.0e+00;
   KSR = 1.20e+00;
   vmaxf = 7.480e-05;
   vmaxr = 3.180e-04;
   v1 = 1.80e+00;
   tautr = 5.7470e-01;
   tauxfer = 2.670e+01;
   kltrpn_plus = 4.0e+01;
   kltrpn_minus = 4.0e-02;
   khtrpn_plus = 2.0e+01;
   khtrpn_minus = 6.60e-05;
   LTRPNtot = 7.00e-02;
   HTRPNtot = 1.40e-01;
   CMDNtot = 5.0e-02;
   KmCMDN = 2.380e-03;
   EGTAtot = 0.0e+00;
   KmEGTA = 1.50e-04;
   CSQNtot = 1.50e+01;
   KmCSQN = 8.0e-01;
   kaplus = 1.2150e-02;
   ncoop = 4;
   kaminus = 5.760e-01;
   kbplus = 4.050e-03;
   mcoop = 3.0e+00;
   kbminus = 1.930e+00;
   kcplus = 1.0e-01;
   kcminus = 8.0e-04;
   aL = 2.0e+00;
   bL = 2.0e+00;
   fL = 3.0e-01;
   gL = 4.0e+00;
   alphaa0Kv43 = 5.437080e-01;
   aaKv43 = 2.898300e-02;
   betaa0Kv43 = 8.01850e-02;
   baKv43 = 4.684370e-02;
   alphai0Kv43 = 4.984240e-02;
   aiKv43 = 3.730160e-04;
   betai0Kv43 = 8.194820e-04;
   biKv43 = 5.3740e-08;
   b1Kv43 = 6.773480e+00;
   b2Kv43 = 1.562127051520e+01;
   b3Kv43 = 2.875326033130e+01;
   b4Kv43 = 5.245762066790001e+02;
   f1Kv43 = 1.893600e+00;
   f2Kv43 = 1.422464745600e+01;
   f3Kv43 = 1.585743783890e+02;
   f4Kv43 = 1.429366453510e+02;
   alphaa0Kv14 = 1.840024145540e+00;
   aaKv14 = 7.685480310e-03;
   betaa0Kv14 = 1.0817483400e-02;
   baKv14 = 7.7933781740e-02;
   alphai0Kv14 = 3.0576791600e-03;
   betai0Kv14 = 2.449360e-06;
   b1Kv14 = 1.009478471050e+00;
   b2Kv14 = 1.171005405670e+00;
   b3Kv14 = 6.390276875799999e-01;
   b4Kv14 = 2.120353790950e+00;
   f1Kv14 = 5.246507399600e-01;
   f2Kv14 = 1.7518854086390e+01;
   f3Kv14 = 9.385876453455600e+02;
   f4Kv14 = 5.474919473332601e+04;
   Nao = 1.380e+02;
   Ko = 4.0e+00;
   Cao = 2.0e+00;
   Kv43Frac = 8.890e-01;
   KvScale = 8.720e-01;
   GNa = 5.6320e+01;
   GKr = 1.860e-02;
   GKs = 3.50e-03;
   GK1 = 1.253051261188080e-01;
   GNab = 1.0e-03;
   GCab = 7.6840e-05;
   IpCamax = 5.0e-02;
   KmpCa = 5.0e-04;
   KmKo = 1.50e+00;
   KmNai = 2.0e+01;
   INaKmax = 2.3870e+00;
   eta = 3.50e-01;
   ksat = 2.0e-01;
   KmCa = 1.380e+00;
   KmNa = 8.750e+01;
   kNaCa = 4.40e-01;
   Pscale = 7.0e+00;
   ICahalf = -2.650e-01;
   RToverF = 0.0e+00;
   TNa = 2.941600e+02;
   T_Const_HERG = 5.320e+00;
   A0_HERG = 1.714764173308600e-02;
   B0_HERG = 3.3046080388350e-02;
   A1_HERG = 3.9693283811410e-02;
   B1_HERG = -4.306054163980e-02;
   A2_HERG = 2.0574486059770e-02;
   B2_HERG = 2.6174127151180e-02;
   A3_HERG = 1.343666044230e-03;
   B3_HERG = -2.6913854983990e-02;
   A4_HERG = 1.066631649128800e-01;
   B4_HERG = 5.689088597170e-03;
   A5_HERG = 6.463939100490e-03;
   B5_HERG = -4.5366429595430e-02;
   A6_HERG = 8.039374403000001e-05;
   B6_HERG = 6.9808924000e-07;
          
   // Mathematics                      
   calc_RT_over_F = ((Rgas*Temp)/Faraday);   //0
   calc_C = (Acap*1.0e-03); //1
   calc_fb = pow((Cai/Kfb),Nfb);  //14
   calc_rb = pow((CaNSR/Krb),Nrb);   //15
   calc_Jrel = (v1*(O1_RyR+O2_RyR)*(CaJSR-CaSS));   //17
   calc_Jtr = ((CaNSR-CaJSR)/tautr);  //18
   calc_Jxfer = ((CaSS-Cai)/tauxfer); //19
   calc_dLTRPNCa = ((kltrpn_plus*Cai*(1.0e+00-LTRPNCa))-(kltrpn_minus*LTRPNCa));  //20
   calc_dHTRPNCa = ((khtrpn_plus*Cai*(1.0e+00-HTRPNCa))-(khtrpn_minus*HTRPNCa));  //21
   calc_beta_SS = (1.0e+00/(1.0e+00+((CMDNtot*KmCMDN)/pow((CaSS+KmCMDN),2))+((EGTAtot*KmEGTA)/pow((CaSS+KmEGTA),2))));  //25
   calc_beta_JSR = (1.0e+00/(1.0e+00+((CSQNtot*KmCSQN)/pow((CaJSR+KmCSQN),2))));  //26
   calc_beta_i = (1.0e+00/(1.0e+00+((CMDNtot*KmCMDN)/pow((Cai+KmCMDN),2))+((EGTAtot*KmEGTA)/pow((Cai+KmEGTA),2))));  //27
   calc_dC1_RyR = (((-kaplus)*pow((CaSS*1.0e+03),ncoop)*C1_RyR)+(kaminus*O1_RyR));   //28
   calc_dO2_RyR = ((kbplus*pow((CaSS*1.0e+03),mcoop)*O1_RyR)-(kbminus*O2_RyR));   //29
   calc_dC2_RyR = ((kcplus*O1_RyR)-(kcminus*C2_RyR));   //30
   calc_alpha = (4.0e+00*1.20e+00*4.160e-01*exp((1.20e-02*(V-3.50e+01))));   //36
   calc_beta = (4.0e+00*4.50e-01*4.90e-02*exp(((-6.50e-02)*(V-2.20e+01)))); //37
   calc_gamma = (6.0e-01*9.2330e-02*CaSS);   //40
   calc_omega = (2.50e-01*1.0e-02); //41
   calc_a1_Cainf = 8.20e-01; //99
   calc_tau_yCa = (1.0e+00/((3.363362094520e-03/(5.0e-01+exp((V/(-5.538998740360550e+00)))))+(7.790465707370e-03*exp((V/(-4.951039631160386e+01)))))); //101
   calc_alpha_act43 = (alphaa0Kv43*exp((aaKv43*V)));   //103
   calc_beta_act43 = (betaa0Kv43*exp(((-baKv43)*V)));  //104
   calc_alpha_inact43 = (alphai0Kv43*exp(((-aiKv43)*V))); //105
   calc_beta_inact43 = (betai0Kv43*exp((biKv43*V)));   //106
   calc_alpha_act14 = (alphaa0Kv14*exp((aaKv14*V)));   //163
   calc_beta_act14 = (betaa0Kv14*exp(((-baKv14)*V)));  //164
   calc_alpha_inact14 = alphai0Kv14;   //165
   calc_beta_inact14 = betai0Kv14;  //166
   calc_a1_duplicated_COMPUTE_REVERSAL_POTENTIALS = (Ko+(1.8330e-02*Nao));  //225
   calc_a2_duplicated_COMPUTE_REVERSAL_POTENTIALS = (Ki+(1.8330e-02*Nai)); //226
   calc_PKv14 = ((1.0e+00-Kv43Frac)*KvScale*4.298600e-07);   //229
   calc_GKv43 = (Kv43Frac*KvScale*1.0e-01);   //230
   calc_fKo = pow((Ko/4.0e+00),1.0/2.0);   //232
   calc_IpCa = ((IpCamax*Cai)/(KmpCa+Cai)); //249
   calc_sigma = ((exp((Nao/6.730e+01))-1.0e+00)/7.0e+00);  //251
   calc_a1_K_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa = (Ko/(Ko+KmKo));  //255
   calc_a2_K = (1.0e+00+sqrt((KmNai/Nai))); //256
   calc_a4_ncx = (KmCa+Cao);  //261
   calc_a5_ncx = ((pow(KmNa,3)+pow(Nao,3))/5.0e+03);  //262
   calc_PK = (Pscale*4.5740e-07); //264
   calc_PCa = (Pscale*2.4690e-04);   //265
   calc_FoverRT = (1.0e+00/RToverF); //277
   calc_KToverH = ((1.3810e-23*TNa)/6.6260e-31);  //278
   calc_RTNa = (Rgas*TNa); //279
   calc_RTNaF = ((Rgas*TNa)/Faraday);  //280
   calc_Temp_Scale = 1.388622912528710e+00;  //281
   calc_parameter_a = 1.400426254774010e+00; //296
   calc_C2H_to_C3H = (T_Const_HERG*2.6083620433370e-02);  //346
   calc_C3H_to_C2H = (T_Const_HERG*1.483297813214500e-01);  //347
   calc_C1H_to_C2H = (T_Const_HERG*A0_HERG*exp((B0_HERG*V)));   //348
   calc_C2H_to_C1H = (T_Const_HERG*A1_HERG*exp((B1_HERG*V)));   //349
   calc_C3H_to_OH = (T_Const_HERG*A2_HERG*exp((B2_HERG*V))); //350
   calc_OH_to_C3H = (T_Const_HERG*A3_HERG*exp((B3_HERG*V))); //351
   calc_OH_to_IH = (T_Const_HERG*A4_HERG*exp((B4_HERG*V)));  //352
   calc_IH_to_OH = (T_Const_HERG*A5_HERG*exp((B5_HERG*V)));  //353
   calc_C3H_to_IH = (T_Const_HERG*A6_HERG*exp((B6_HERG*V))); //354
   calc_C0ks_C1ks = 7.956007980040e-03;   //369
   calc_C1ks_O1ks = 3.9667206760710e-02;   //370
   calc_O1ks_O2ks = (7.672543630630e-03*exp((8.6629459146550e-02*V))); //371
   calc_O1ks_C1ks = (7.008066289290e-03*exp(((-1.499975470028500e-01)*V))); //372
   calc_O2ks_O1ks = (3.797379983680e-03*exp(((-1.4256681268810e-02)*V))); //373
   calc_C1ks_C0ks = (2.1625575895850e-01*exp(((-1.8891230210e-05)*V))); //374
   calc_a1 = (calc_C/(Vmyo*Faraday));  //2
   calc_a2 = (calc_C/(2.0e+00*VSS*Faraday));  //3
   calc_i_Stim = ((time >= stim_start) && (time <= stim_end) && (((time - stim_start) - (floor((time - stim_start) / stim_period) * stim_period)) <= stim_duration)) ? stim_amplitude : 0.0;
   calc_Jup = ((KSR*((vmaxf*calc_fb)-(vmaxr*calc_rb)))/(1.0e+00+calc_fb+calc_rb));   //16
   calc_Jtrpn = ((LTRPNtot*calc_dLTRPNCa)+(HTRPNtot*calc_dHTRPNCa)); //24
   calc_dO1_RyR = (-(calc_dC1_RyR+calc_dO2_RyR+calc_dC2_RyR)); //31
   calc_alpha_prime = (aL*calc_alpha); //38
   calc_beta_prime = (calc_beta/bL);   //39
   calc_C0_to_C1 = (4.0e+00*calc_alpha);   //42
   calc_C1_to_C2 = (3.0e+00*calc_alpha);   //43
   calc_C2_to_C3 = (2.0e+00*calc_alpha);   //44
   calc_C3_to_C4 = calc_alpha;   //45
   calc_C1_to_C0 = calc_beta; //50
   calc_C2_to_C1 = (2.0e+00*calc_beta); //51
   calc_C3_to_C2 = (3.0e+00*calc_beta); //52
   calc_C4_to_C3 = (4.0e+00*calc_beta); //53
   calc_C0_to_CCa0 = calc_gamma; //58
   calc_CCa0_to_C0 = calc_omega; //63
   calc_yCa_inf = (((calc_a1_Cainf/(1.0e+00+exp(((V+2.850e+01)/7.80e+00))))+1.0e+00)-calc_a1_Cainf); //100
   calc_C0Kv43_to_C1Kv43 = (4.0e+00*calc_alpha_act43); //107
   calc_C1Kv43_to_C2Kv43 = (3.0e+00*calc_alpha_act43); //108
   calc_C2Kv43_to_C3Kv43 = (2.0e+00*calc_alpha_act43); //109
   calc_C3Kv43_to_OKv43 = calc_alpha_act43;  //110
   calc_CI0Kv43_to_CI1Kv43 = (4.0e+00*b1Kv43*calc_alpha_act43); //111
   calc_CI1Kv43_to_CI2Kv43 = ((3.0e+00*b2Kv43*calc_alpha_act43)/b1Kv43); //112
   calc_CI2Kv43_to_CI3Kv43 = ((2.0e+00*b3Kv43*calc_alpha_act43)/b2Kv43); //113
   calc_CI3Kv43_to_OIKv43 = ((b4Kv43*calc_alpha_act43)/b3Kv43);   //114
   calc_C1Kv43_to_C0Kv43 = calc_beta_act43;  //115
   calc_C2Kv43_to_C1Kv43 = (2.0e+00*calc_beta_act43);  //116
   calc_C3Kv43_to_C2Kv43 = (3.0e+00*calc_beta_act43);  //117
   calc_OKv43_to_C3Kv43 = (4.0e+00*calc_beta_act43);   //118
   calc_CI1Kv43_to_CI0Kv43 = (calc_beta_act43/f1Kv43);   //119
   calc_CI2Kv43_to_CI1Kv43 = ((2.0e+00*f1Kv43*calc_beta_act43)/f2Kv43);  //120
   calc_CI3Kv43_to_CI2Kv43 = ((3.0e+00*f2Kv43*calc_beta_act43)/f3Kv43);  //121
   calc_OIKv43_to_CI3Kv43 = ((4.0e+00*f3Kv43*calc_beta_act43)/f4Kv43);   //122
   calc_C0Kv43_to_CI0Kv43 = calc_beta_inact43;  //123
   calc_C1Kv43_to_CI1Kv43 = (f1Kv43*calc_beta_inact43);  //124
   calc_C2Kv43_to_CI2Kv43 = (f2Kv43*calc_beta_inact43);  //125
   calc_C3Kv43_to_CI3Kv43 = (f3Kv43*calc_beta_inact43);  //126
   calc_OKv43_to_OIKv43 = (f4Kv43*calc_beta_inact43); //127
   calc_CI0Kv43_to_C0Kv43 = calc_alpha_inact43; //128
   calc_CI1Kv43_to_C1Kv43 = (calc_alpha_inact43/b1Kv43); //129
   calc_CI2Kv43_to_C2Kv43 = (calc_alpha_inact43/b2Kv43); //130
   calc_CI3Kv43_to_C3Kv43 = (calc_alpha_inact43/b3Kv43); //131
   calc_OIKv43_to_OKv43 = (calc_alpha_inact43/b4Kv43);   //132
   calc_C0Kv14_to_C1Kv14 = (4.0e+00*calc_alpha_act14); //167
   calc_C1Kv14_to_C2Kv14 = (3.0e+00*calc_alpha_act14); //168
   calc_C2Kv14_to_C3Kv14 = (2.0e+00*calc_alpha_act14); //169
   calc_C3Kv14_to_OKv14 = calc_alpha_act14;  //170
   calc_CI0Kv14_to_CI1Kv14 = (4.0e+00*b1Kv14*calc_alpha_act14); //171
   calc_CI1Kv14_to_CI2Kv14 = ((3.0e+00*b2Kv14*calc_alpha_act14)/b1Kv14); //172
   calc_CI2Kv14_to_CI3Kv14 = ((2.0e+00*b3Kv14*calc_alpha_act14)/b2Kv14); //173
   calc_CI3Kv14_to_OIKv14 = ((b4Kv14*calc_alpha_act14)/b3Kv14);   //174
   calc_C1Kv14_to_C0Kv14 = calc_beta_act14;  //175
   calc_C2Kv14_to_C1Kv14 = (2.0e+00*calc_beta_act14);  //176
   calc_C3Kv14_to_C2Kv14 = (3.0e+00*calc_beta_act14);  //177
   calc_OKv14_to_C3Kv14 = (4.0e+00*calc_beta_act14);   //178
   calc_CI1Kv14_to_CI0Kv14 = (calc_beta_act14/f1Kv14);   //179
   calc_CI2Kv14_to_CI1Kv14 = ((2.0e+00*f1Kv14*calc_beta_act14)/f2Kv14);  //180
   calc_CI3Kv14_to_CI2Kv14 = ((3.0e+00*f2Kv14*calc_beta_act14)/f3Kv14);  //181
   calc_OIKv14_to_CI3Kv14 = ((4.0e+00*f3Kv14*calc_beta_act14)/f4Kv14);   //182
   calc_C0Kv14_to_CI0Kv14 = calc_beta_inact14;  //183
   calc_C1Kv14_to_CI1Kv14 = (f1Kv14*calc_beta_inact14);  //184
   calc_C2Kv14_to_CI2Kv14 = (f2Kv14*calc_beta_inact14);  //185
   calc_C3Kv14_to_CI3Kv14 = (f3Kv14*calc_beta_inact14);  //186
   calc_OKv14_to_OIKv14 = (f4Kv14*calc_beta_inact14); //187
   calc_CI0Kv14_to_C0Kv14 = calc_alpha_inact14; //188
   calc_CI1Kv14_to_C1Kv14 = (calc_alpha_inact14/b1Kv14); //189
   calc_CI2Kv14_to_C2Kv14 = (calc_alpha_inact14/b2Kv14); //190
   calc_CI3Kv14_to_C3Kv14 = (calc_alpha_inact14/b3Kv14); //191
   calc_OIKv14_to_OKv14 = (calc_alpha_inact14/b4Kv14);   //192
   calc_ENa = (calc_RT_over_F*log((Nao/Nai)));   //223
   calc_EK = (calc_RT_over_F*log((Ko/Ki)));   //224
   calc_EKs = (calc_RT_over_F*log((calc_a1_duplicated_COMPUTE_REVERSAL_POTENTIALS/calc_a2_duplicated_COMPUTE_REVERSAL_POTENTIALS)));   //227
   calc_ECa = (5.0e-01*calc_RT_over_F*log((Cao/Cai)));  //228
   calc_VF_over_RT = (V/calc_RT_over_F);   //236
   calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa = (V/calc_RT_over_F);  //250
   calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK = (V/calc_RT_over_F);  //266
   calc_alpha1 = (calc_Temp_Scale*calc_KToverH*exp((((-1.140074627002320e+05)/calc_RTNa)+(2.241140e+02/Rgas)+((2.863742685962350e-01*V)/calc_RTNaF))));   //282
   calc_beta1 = (calc_Temp_Scale*calc_KToverH*exp((((-2.724702734896810e+05)/calc_RTNa)+(7.081460e+02/Rgas)+(((-2.285284175864240e+00)*V)/calc_RTNaF)))); //283
   calc_gamma1 = (calc_Temp_Scale*calc_KToverH*exp((((-1.963365757359230e+05)/calc_RTNa)+(5.299520e+02/Rgas)+((2.780849185960450e+00*V)/calc_RTNaF))));   //284
   calc_Delta1 = (calc_Temp_Scale*calc_KToverH*exp((((-1.3368993040910e+05)/calc_RTNa)+(2.292050e+02/Rgas)+(((-1.558042145538830e+00)*V)/calc_RTNaF))));   //285
   calc_On = (calc_Temp_Scale*calc_KToverH*exp((((-6.212307843804810e+04)/calc_RTNa)+(3.929500e+01/Rgas)+((2.888160427432320e-01*V)/calc_RTNaF)))); //286
   calc_Of = (calc_Temp_Scale*calc_KToverH*exp((((-9.765784971370149e+04)/calc_RTNa)+(1.510e+00/Rgas)+((6.848619931006850e-02*V)/calc_RTNaF)))); //287
   calc_GammaGamma = (calc_Temp_Scale*calc_KToverH*exp(((1.164311421423480e+05/calc_RTNa)+((-5.783170e+02)/Rgas)+((7.641260117457070e-01*V)/calc_RTNaF))));  //288
   calc_DeltaDelta = (calc_Temp_Scale*calc_KToverH*exp((((-5.570066246583070e+04)/calc_RTNa)+((-1.306390e+02)/Rgas)+(((-3.649816729270780e+00)*V)/calc_RTNaF))));  //289
   calc_epsilon = (calc_Temp_Scale*calc_KToverH*exp((((-8.580036755783259e+04)/calc_RTNa)+(7.007800e+01/Rgas)))); //290
   calc_omega_na = (calc_Temp_Scale*calc_KToverH*exp((((-1.219551661548640e+05)/calc_RTNa)+(2.251750e+02/Rgas))));   //291
   calc_rho = (calc_Temp_Scale*calc_KToverH*exp((((-1.478139900050350e+05)/calc_RTNa)+(3.389150e+02/Rgas)+((2.136004370212600e+00*V)/calc_RTNaF))));   //292
   calc_mu = (calc_Temp_Scale*calc_KToverH*exp((((-1.213221432752420e+05)/calc_RTNa)+(1.932650e+02/Rgas)+(((-1.742902670209030e+00)*V)/calc_RTNaF)))); //293
   calc_Cn = (calc_Temp_Scale*calc_KToverH*exp((((-2.879134465309530e+05)/calc_RTNa)+(7.862170e+02/Rgas))));   //294
   calc_Cf = (calc_Temp_Scale*calc_KToverH*exp((((-5.956522362845840e+04)/calc_RTNa)+(7.110e-03/Rgas))));   //295
   calc_IH_to_C3H = ((calc_OH_to_C3H*calc_IH_to_OH*calc_C3H_to_IH)/(calc_C3H_to_OH*calc_OH_to_IH));   //355
   calc_a1_C2_duplicated_IKr = ((calc_C1H_to_C2H*C1Herg)+(calc_C3H_to_C2H*C3Herg));   //357
   calc_a2_C2_duplicated_IKr = ((calc_C2H_to_C1H+calc_C2H_to_C3H)*C2Herg);  //358
   calc_a2_C3_duplicated_IKr = ((calc_C3H_to_IH+calc_C3H_to_OH+calc_C3H_to_C2H)*C3Herg);   //361
   calc_a1_O_duplicated_IKr = ((calc_C3H_to_OH*C3Herg)+(calc_IH_to_OH*IHerg));  //363
   calc_a2_O_duplicated_IKr = ((calc_OH_to_C3H+calc_OH_to_IH)*OHerg); //364
   calc_a1_I = ((calc_C3H_to_IH*C3Herg)+(calc_OH_to_IH*OHerg));  //366
   calc_CCa0_to_CCa1 = (4.0e+00*calc_alpha_prime);  //46
   calc_CCa1_to_CCa2 = (3.0e+00*calc_alpha_prime);  //47
   calc_CCa2_to_CCa3 = (2.0e+00*calc_alpha_prime);  //48
   calc_CCa3_to_CCa4 = calc_alpha_prime;  //49
   calc_CCa1_to_CCa0 = calc_beta_prime;   //54
   calc_CCa2_to_CCa1 = (2.0e+00*calc_beta_prime);   //55
   calc_CCa3_to_CCa2 = (3.0e+00*calc_beta_prime);   //56
   calc_CCa4_to_CCa3 = (4.0e+00*calc_beta_prime);   //57
   calc_C1_to_CCa1 = (aL*calc_C0_to_CCa0);   //59
   calc_CCa1_to_C1 = (calc_CCa0_to_C0/bL);   //64
   calc_a1_C0 = ((calc_C0_to_C1+calc_C0_to_CCa0)*C0);  //68
   calc_a2_C0 = ((calc_C1_to_C0*C1)+(calc_CCa0_to_C0*CCa0));  //69
   calc_a1_C043 = ((calc_C0Kv43_to_C1Kv43+calc_C0Kv43_to_CI0Kv43)*C0Kv43);  //133
   calc_a2_C043 = ((calc_C1Kv43_to_C0Kv43*C1Kv43)+(calc_CI0Kv43_to_C0Kv43*CI0Kv43));  //134
   calc_a1_C143 = ((calc_C1Kv43_to_C2Kv43+calc_C1Kv43_to_C0Kv43+calc_C1Kv43_to_CI1Kv43)*C1Kv43); //136
   calc_a2_C143 = ((calc_C2Kv43_to_C1Kv43*C2Kv43)+(calc_CI1Kv43_to_C1Kv43*CI1Kv43)+(calc_C0Kv43_to_C1Kv43*C0Kv43));  //137
   calc_a1_C243 = ((calc_C2Kv43_to_C3Kv43+calc_C2Kv43_to_C1Kv43+calc_C2Kv43_to_CI2Kv43)*C2Kv43); //139
   calc_a2_C243 = ((calc_C3Kv43_to_C2Kv43*C3Kv43)+(calc_CI2Kv43_to_C2Kv43*CI2Kv43)+(calc_C1Kv43_to_C2Kv43*C1Kv43));  //140
   calc_a1_C343 = ((calc_C3Kv43_to_OKv43+calc_C3Kv43_to_C2Kv43+calc_C3Kv43_to_CI3Kv43)*C3Kv43);  //142
   calc_a2_C343 = ((calc_OKv43_to_C3Kv43*OKv43)+(calc_CI3Kv43_to_C3Kv43*CI3Kv43)+(calc_C2Kv43_to_C3Kv43*C2Kv43)); //143
   calc_a1_O43 = ((calc_OKv43_to_C3Kv43+calc_OKv43_to_OIKv43)*OKv43); //145
   calc_a2_O43 = ((calc_C3Kv43_to_OKv43*C3Kv43)+(calc_OIKv43_to_OKv43*OIKv43)); //146
   calc_a1_I043 = ((calc_CI0Kv43_to_C0Kv43+calc_CI0Kv43_to_CI1Kv43)*CI0Kv43);  //148
   calc_a2_I043 = ((calc_C0Kv43_to_CI0Kv43*C0Kv43)+(calc_CI1Kv43_to_CI0Kv43*CI1Kv43));   //149
   calc_a1_I143 = ((calc_CI1Kv43_to_CI2Kv43+calc_CI1Kv43_to_C1Kv43+calc_CI1Kv43_to_CI0Kv43)*CI1Kv43);  //151
   calc_a2_I143 = ((calc_CI2Kv43_to_CI1Kv43*CI2Kv43)+(calc_C1Kv43_to_CI1Kv43*C1Kv43)+(calc_CI0Kv43_to_CI1Kv43*CI0Kv43));   //152
   calc_a1_I243 = ((calc_CI2Kv43_to_CI3Kv43+calc_CI2Kv43_to_C2Kv43+calc_CI2Kv43_to_CI1Kv43)*CI2Kv43);  //154
   calc_a2_I243 = ((calc_CI3Kv43_to_CI2Kv43*CI3Kv43)+(calc_C2Kv43_to_CI2Kv43*C2Kv43)+(calc_CI1Kv43_to_CI2Kv43*CI1Kv43));   //155
   calc_a1_I343 = ((calc_CI3Kv43_to_OIKv43+calc_CI3Kv43_to_C3Kv43+calc_CI3Kv43_to_CI2Kv43)*CI3Kv43);   //157
   calc_a2_I343 = ((calc_OIKv43_to_CI3Kv43*OIKv43)+(calc_C3Kv43_to_CI3Kv43*C3Kv43)+(calc_CI2Kv43_to_CI3Kv43*CI2Kv43));  //158
   calc_a1_OI43 = ((calc_OIKv43_to_OKv43+calc_OIKv43_to_CI3Kv43)*OIKv43);   //160
   calc_a2_OI43 = ((calc_OKv43_to_OIKv43*OKv43)+(calc_CI3Kv43_to_OIKv43*CI3Kv43)); //161
   calc_a1_C0_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C0Kv14_to_C1Kv14+calc_C0Kv14_to_CI0Kv14)*C0Kv14);   //193
   calc_a2_C0_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C1Kv14_to_C0Kv14*C1Kv14)+(calc_CI0Kv14_to_C0Kv14*CI0Kv14));   //194
   calc_a1_C1_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C1Kv14_to_C2Kv14+calc_C1Kv14_to_C0Kv14+calc_C1Kv14_to_CI1Kv14)*C1Kv14);  //196
   calc_a2_C1_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C2Kv14_to_C1Kv14*C2Kv14)+(calc_CI1Kv14_to_C1Kv14*CI1Kv14)+(calc_C0Kv14_to_C1Kv14*C0Kv14));   //197
   calc_a1_C2_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C2Kv14_to_C3Kv14+calc_C2Kv14_to_C1Kv14+calc_C2Kv14_to_CI2Kv14)*C2Kv14);  //199
   calc_a2_C2_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C3Kv14_to_C2Kv14*C3Kv14)+(calc_CI2Kv14_to_C2Kv14*CI2Kv14)+(calc_C1Kv14_to_C2Kv14*C1Kv14));   //200
   calc_a1_C3_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_C3Kv14_to_OKv14+calc_C3Kv14_to_C2Kv14+calc_C3Kv14_to_CI3Kv14)*C3Kv14);   //202
   calc_a2_C3_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES = ((calc_OKv14_to_C3Kv14*OKv14)+(calc_CI3Kv14_to_C3Kv14*CI3Kv14)+(calc_C2Kv14_to_C3Kv14*C2Kv14));  //203
   calc_a1_O = ((calc_OKv14_to_C3Kv14+calc_OKv14_to_OIKv14)*OKv14);   //205
   calc_a2_O = ((calc_C3Kv14_to_OKv14*C3Kv14)+(calc_OIKv14_to_OKv14*OIKv14));   //206
   calc_a1_CI0 = ((calc_CI0Kv14_to_C0Kv14+calc_CI0Kv14_to_CI1Kv14)*CI0Kv14);   //208
   calc_a2_CI0 = ((calc_C0Kv14_to_CI0Kv14*C0Kv14)+(calc_CI1Kv14_to_CI0Kv14*CI1Kv14)); //209
   calc_a1_CI1 = ((calc_CI1Kv14_to_CI2Kv14+calc_CI1Kv14_to_C1Kv14+calc_CI1Kv14_to_CI0Kv14)*CI1Kv14);   //211
   calc_a2_CI1 = ((calc_CI2Kv14_to_CI1Kv14*CI2Kv14)+(calc_C1Kv14_to_CI1Kv14*C1Kv14)+(calc_CI0Kv14_to_CI1Kv14*CI0Kv14)); //212
   calc_a1_CI2 = ((calc_CI2Kv14_to_CI3Kv14+calc_CI2Kv14_to_C2Kv14+calc_CI2Kv14_to_CI1Kv14)*CI2Kv14);   //214
   calc_a2_CI2 = ((calc_CI3Kv14_to_CI2Kv14*CI3Kv14)+(calc_C2Kv14_to_CI2Kv14*C2Kv14)+(calc_CI1Kv14_to_CI2Kv14*CI1Kv14)); //215
   calc_a1_CI3 = ((calc_CI3Kv14_to_OIKv14+calc_CI3Kv14_to_C3Kv14+calc_CI3Kv14_to_CI2Kv14)*CI3Kv14); //217
   calc_a2_CI3 = ((calc_OIKv14_to_CI3Kv14*OIKv14)+(calc_C3Kv14_to_CI3Kv14*C3Kv14)+(calc_CI2Kv14_to_CI3Kv14*CI2Kv14));   //218
   calc_a1_OI = ((calc_OIKv14_to_OKv14+calc_OIKv14_to_CI3Kv14)*OIKv14);  //220
   calc_a2_OI = ((calc_OKv14_to_OIKv14*OKv14)+(calc_CI3Kv14_to_OIKv14*CI3Kv14));   //221
   calc_INa = (GNa*(na6+na7)*(V-calc_ENa));  //231
   calc_IKr = (GKr*calc_fKo*OHerg*(V-calc_EK));   //233
   calc_IKs = (GKs*(O1ks+O2ks)*(V-calc_EK)); //234
   calc_IKv43 = (calc_GKv43*OKv43*(V-calc_EK));   //235
   calc_VFsq_over_RT = (1.0e+03*Faraday*calc_VF_over_RT); //237
   calc_a1_K = ((Ki*exp(calc_VF_over_RT))-Ko);   //238
   calc_a2_duplicated_COMPUTE_INa_IKr_IKs_Ito1_IK1_INab_IKp = (exp(calc_VF_over_RT)-1.0e+00); //239
   calc_a1_Na = ((Nai*exp(calc_VF_over_RT))-Nao);   //241
   calc_K1_inf = (1.0e+00/(9.399999999999999e-01+exp(((1.260e+00/calc_RT_over_F)*(V-calc_EK))))); //245
   calc_INab = (GNab*(V-calc_ENa));  //247
   calc_ICab = (GCab*(V-calc_ECa));  //248
   calc_a1_Na_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa = (1.0e+00+(1.2450e-01*exp(((-1.0e-01)*calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa))));   //252
   calc_a2_Na = (3.650e-02*calc_sigma*exp(((-1.330e+00)*calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa)));  //253
   calc_a1_ncx = (exp((eta*calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa))*pow(Nai,3)*Cao);   //258
   calc_a2_ncx = (exp(((eta-1.0e+00)*calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa))*pow(Nao,3)*Cai);   //259
   calc_a3_ncx = (1.0e+00+(ksat*exp(((eta-1.0e+00)*calc_VF_over_RT_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa))));  //260
   calc_VFsq_over_RT_duplicated_COMPUTE_ICa_ICaK = (1.0e+03*Faraday*calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK);  //267
   calc_a1_Ca = ((1.0e-03*exp((2.0e+00*calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK)))-(Cao*3.410e-01)); //268
   calc_a2_Ca = (exp((2.0e+00*calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK))-1.0e+00); //269
   calc_a1_K_duplicated_COMPUTE_ICa_ICaK = ((Ki*exp(calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK))-Ko); //274
   calc_a2_K_duplicated_COMPUTE_ICa_ICaK = (exp(calc_VF_over_RT_duplicated_COMPUTE_ICa_ICaK)-1.0e+00); //275
   calc_k12 = (4.0e+00*calc_alpha1); //297
   calc_k23 = (3.0e+00*calc_alpha1); //298
   calc_k34 = (2.0e+00*calc_alpha1); //299
   calc_k45 = calc_alpha1; //300
   calc_k56 = calc_gamma1; //301
   calc_k67 = calc_epsilon;   //302
   calc_k1213 = calc_GammaGamma; //307
   calc_k57 = calc_rho; //308
   calc_k21 = calc_beta1;  //309
   calc_k32 = (2.0e+00*calc_beta1);  //310
   calc_k43 = (3.0e+00*calc_beta1);  //311
   calc_k54 = (4.0e+00*calc_beta1);  //312
   calc_k65 = calc_Delta1; //313
   calc_k76 = calc_omega_na;  //314
   calc_k1312 = calc_DeltaDelta; //319
   calc_k75 = calc_mu;  //320
   calc_k81 = calc_Cf;  //321
   calc_k136 = calc_Of; //326
   calc_k18 = calc_Cn;  //327
   calc_k613 = calc_On; //332
   calc_a1_C3_duplicated_IKr = ((calc_C2H_to_C3H*C2Herg)+(calc_OH_to_C3H*OHerg)+(calc_IH_to_C3H*IHerg)); //360
   calc_a2_I = ((calc_IH_to_C3H+calc_IH_to_OH)*IHerg); //367
   calc_C2_to_CCa2 = (aL*calc_C1_to_CCa1);   //60
   calc_CCa2_to_C2 = (calc_CCa1_to_C1/bL);   //65
   calc_a1_C1 = ((calc_C1_to_C0+calc_C1_to_C2+calc_C1_to_CCa1)*C1);   //71
   calc_a2_C1 = ((calc_C0_to_C1*C0)+(calc_C2_to_C1*C2)+(calc_CCa1_to_C1*CCa1));  //72
   calc_a1_Ca0 = ((calc_CCa0_to_CCa1+calc_CCa0_to_C0)*CCa0); //84
   calc_a2_Ca0 = ((calc_CCa1_to_CCa0*CCa1)+(calc_C0_to_CCa0*C0));   //85
   calc_a1_Ca1 = ((calc_CCa1_to_CCa0+calc_CCa1_to_CCa2+calc_CCa1_to_C1)*CCa1); //87
   calc_a2_Ca1 = ((calc_CCa0_to_CCa1*CCa0)+(calc_CCa2_to_CCa1*CCa2)+(calc_C1_to_CCa1*C1));   //88
   calc_IKv14_K = ((calc_PKv14*OKv14*calc_VFsq_over_RT*calc_a1_K)/calc_a2_duplicated_COMPUTE_INa_IKr_IKs_Ito1_IK1_INab_IKp);   //240
   calc_IKv14_Na = ((2.0e-02*calc_PKv14*OKv14*calc_VFsq_over_RT*calc_a1_Na)/calc_a2_duplicated_COMPUTE_INa_IKr_IKs_Ito1_IK1_INab_IKp);   //242
   calc_IK1 = (GK1*pow((Ko/1.0e+00),1.0/2.0)*calc_K1_inf*(V-calc_EK)); //246
   calc_fNaK = (1.0e+00/(calc_a1_Na_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa+calc_a2_Na));  //254
   calc_INaCa = ((kNaCa*(calc_a1_ncx-calc_a2_ncx))/(calc_a4_ncx*calc_a3_ncx*calc_a5_ncx));   //263
   calc_ICamax = ((calc_PCa*4.0e+00*calc_VFsq_over_RT_duplicated_COMPUTE_ICa_ICaK*calc_a1_Ca)/calc_a2_Ca);   //270
   calc_k89 = (calc_k12*calc_parameter_a);   //303
   calc_k910 = (calc_k23*calc_parameter_a);  //304
   calc_k1011 = (calc_k34*calc_parameter_a); //305
   calc_k1112 = (calc_k45*calc_parameter_a); //306
   calc_k98 = (calc_k21/calc_parameter_a);   //315
   calc_k109 = (calc_k32/calc_parameter_a);  //316
   calc_k1110 = (calc_k43/calc_parameter_a); //317
   calc_k1211 = (calc_k54/calc_parameter_a); //318
   calc_k92 = (calc_k81/calc_parameter_a);   //322
   calc_k29 = (calc_k18*calc_parameter_a);   //328
   calc_C3_to_CCa3 = (aL*calc_C2_to_CCa2);   //61
   calc_CCa3_to_C3 = (calc_CCa2_to_C2/bL);   //66
   calc_a1_C2 = ((calc_C2_to_C1+calc_C2_to_C3+calc_C2_to_CCa2)*C2);   //74
   calc_a2_C2 = ((calc_C1_to_C2*C1)+(calc_C3_to_C2*C3)+(calc_CCa2_to_C2*CCa2));  //75
   calc_a1_Ca2 = ((calc_CCa2_to_CCa1+calc_CCa2_to_CCa3+calc_CCa2_to_C2)*CCa2); //90
   calc_a2_Ca2 = ((calc_CCa1_to_CCa2*CCa1)+(calc_CCa3_to_CCa2*CCa3)+(calc_C2_to_CCa2*C2));   //91
   calc_IKv14 = (calc_IKv14_K+calc_IKv14_Na);   //243
   calc_INaK = ((INaKmax*calc_fNaK*calc_a1_K_duplicated_COMPUTE_INaK_INaCa_ICab_IpCa)/calc_a2_K);  //257
   calc_ICa = (calc_ICamax*yCa*Open); //271
   calc_Icabar = (calc_ICamax >= 0.0) ? 0.0 : calc_ICamax;
   calc_k103 = (calc_k92/calc_parameter_a);  //323
   calc_k310 = (calc_k29*calc_parameter_a);  //329
   calc_C4_to_CCa4 = (aL*calc_C3_to_CCa3);   //62
   calc_CCa4_to_C4 = (calc_CCa3_to_C3/bL);   //67
   calc_a1_C3 = ((calc_C3_to_C2+calc_C3_to_C4+calc_C3_to_CCa3)*C3);   //77
   calc_a2_C3 = ((calc_C2_to_C3*C2)+(calc_C4_to_C3*C4)+(calc_CCa3_to_C3*CCa3));  //78
   calc_a1_Ca3 = ((calc_CCa3_to_CCa2+calc_CCa3_to_CCa4+calc_CCa3_to_C3)*CCa3); //93
   calc_a2_Ca3 = ((calc_CCa2_to_CCa3*CCa2)+(calc_CCa4_to_CCa3*CCa4)+(calc_C3_to_CCa3*C3));   //94
   calc_Ito1 = (calc_IKv43+calc_IKv14);   //244
   calc_PKprime = (calc_PK/(1.0e+00+(calc_Icabar/ICahalf))); //273
   calc_k114 = (calc_k103/calc_parameter_a); //324
   calc_k411 = (calc_k310*calc_parameter_a); //330
   calc_a1_C4 = ((calc_C4_to_C3+fL+calc_C4_to_CCa4)*C4);  //80
   calc_a2_C4 = ((calc_C3_to_C4*C3)+(gL*Open)+(calc_CCa4_to_C4*CCa4));  //81
   calc_a1_Ca4 = ((calc_CCa4_to_CCa3+calc_CCa4_to_C4)*CCa4); //96
   calc_a2_Ca4 = ((calc_CCa3_to_CCa4*CCa3)+(calc_C4_to_CCa4*C4));   //97
   calc_ICaK = ((calc_PKprime*Open*yCa*calc_VFsq_over_RT_duplicated_COMPUTE_ICa_ICaK*calc_a1_K_duplicated_COMPUTE_ICa_ICaK)/calc_a2_K_duplicated_COMPUTE_ICa_ICaK); //276
   calc_k125 = (calc_k114/calc_parameter_a); //325
   calc_k512 = (calc_k411*calc_parameter_a); //331
   calc_i_tot = (calc_INa+calc_ICa+calc_ICaK+calc_IKr+calc_IKs+calc_IK1+calc_INaCa+calc_INaK+calc_Ito1+calc_IpCa+calc_ICab+calc_INab+calc_i_Stim); //10

   rDY[ 0] = (-calc_i_tot);   // 11
   rDY[ 1] = ((-(((calc_IKr+calc_IKs+calc_IK1+calc_ICaK+calc_i_Stim)-(2.0e+00*calc_INaK))+calc_IKv43+calc_IKv14_K))*calc_a1);  // 5
   rDY[ 2] = (calc_beta_i*(((calc_Jxfer-calc_Jup)-calc_Jtrpn)-(((calc_ICab-(2.0e+00*calc_INaCa))+calc_IpCa)*5.0e-01*calc_a1))); // 6
   rDY[ 3] = (calc_beta_SS*((((calc_Jrel*VJSR)/VSS)-((calc_Jxfer*Vmyo)/VSS))-(calc_ICa*calc_a2))); // 7
   rDY[ 4] = (calc_beta_JSR*(calc_Jtr-calc_Jrel)); // 8
   rDY[ 5] = (((calc_Jup*Vmyo)/VNSR)-((calc_Jtr*VJSR)/VNSR));  // 9
   rDY[ 6] = ((-(calc_INa+calc_INab+(3.0e+00*(calc_INaCa+calc_INaK))+calc_IKv14_Na))*calc_a1);   // 4
   rDY[ 7] = calc_dLTRPNCa;   // 22
   rDY[ 8] = calc_dHTRPNCa;   // 23
   rDY[ 9] = calc_dC1_RyR; // 32
   rDY[10] = calc_dO2_RyR;   // 33
   rDY[11] = calc_dC2_RyR;   // 34
   rDY[12] = calc_dO1_RyR;   // 35
   rDY[13] = (calc_a2_C0-calc_a1_C0); // 70
   rDY[14] = (calc_a2_C1-calc_a1_C1); // 73
   rDY[15] = (calc_a2_C2-calc_a1_C2); // 76
   rDY[16] = (calc_a2_C3-calc_a1_C3); // 79
   rDY[17] = (calc_a2_C4-calc_a1_C4); // 82
   rDY[18] = ((fL*C4)-(gL*Open)); // 83
   rDY[19] = (calc_a2_Ca0-calc_a1_Ca0);  // 86
   rDY[20] = (calc_a2_Ca1-calc_a1_Ca1);  // 89
   rDY[21] = (calc_a2_Ca2-calc_a1_Ca2);  // 92
   rDY[22] = (calc_a2_Ca3-calc_a1_Ca3);  // 95
   rDY[23] = (calc_a2_Ca4-calc_a1_Ca4);  // 98
   rDY[24] = ((calc_yCa_inf-yCa)/calc_tau_yCa); // 102
   rDY[25] = (calc_a2_C043-calc_a1_C043);   // 135
   rDY[26] = (calc_a2_C143-calc_a1_C143);   // 138
   rDY[27] = (calc_a2_C243-calc_a1_C243);   // 141
   rDY[28] = (calc_a2_C343-calc_a1_C343);   // 144
   rDY[29] = (calc_a2_O43-calc_a1_O43);  // 147
   rDY[30] = (calc_a2_I043-calc_a1_I043);   // 150
   rDY[31] = (calc_a2_I143-calc_a1_I143);   // 153
   rDY[32] = (calc_a2_I243-calc_a1_I243);   // 156
   rDY[33] = (calc_a2_I343-calc_a1_I343);   // 159
   rDY[34] = (calc_a2_OI43-calc_a1_OI43);   // 162
   rDY[35] = (calc_a2_C0_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES-calc_a1_C0_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES);  // 195
   rDY[36] = (calc_a2_C1_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES-calc_a1_C1_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES);  // 198
   rDY[37] = (calc_a2_C2_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES-calc_a1_C2_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES);  // 201
   rDY[38] = (calc_a2_C3_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES-calc_a1_C3_duplicated_COMPUTE_DERIVATIVES_OF_Kv1_4_CHANNEL_STATES);  // 204
   rDY[39] = (calc_a2_O-calc_a1_O);   // 207
   rDY[40] = (calc_a2_CI0-calc_a1_CI0);  // 210
   rDY[41] = (calc_a2_CI1-calc_a1_CI1);  // 213
   rDY[42] = (calc_a2_CI2-calc_a1_CI2);  // 216
   rDY[43] = (calc_a2_CI3-calc_a1_CI3);  // 219
   rDY[44] = (calc_a2_OI-calc_a1_OI); // 222
   rDY[45] = (((-(calc_k18+calc_k12))*na1)+(calc_k21*na2)+(calc_k81*na8));   // 333
   rDY[46] = (((calc_k12*na1)-((calc_k21+calc_k23+calc_k29)*na2))+(calc_k32*na3)+(calc_k92*na9));  // 334
   rDY[47] = (((calc_k23*na2)-((calc_k32+calc_k34+calc_k310)*na3))+(calc_k43*na4)+(calc_k103*na10));  // 335
   rDY[48] = (((calc_k34*na3)-((calc_k43+calc_k45+calc_k411)*na4))+(calc_k54*na5)+(calc_k114*na11));  // 336
   rDY[49] = (((calc_k45*na4)-((calc_k54+calc_k56+calc_k57+calc_k512)*na5))+(calc_k65*na6)+(calc_k75*na7)+(calc_k125*na12));   // 337
   rDY[50] = (((calc_k56*na5)-((calc_k65+calc_k67+calc_k613)*na6))+(calc_k76*na7)+(calc_k136*na13));  // 338
   rDY[51] = (((calc_k57*na5)+(calc_k67*na6))-((calc_k75+calc_k76)*na7)); // 339
   rDY[52] = (((calc_k18*na1)-((calc_k81+calc_k89)*na8))+(calc_k98*na9)); // 340
   rDY[53] = ((((calc_k29*na2)+(calc_k89*na8))-((calc_k98+calc_k92+calc_k910)*na9))+(calc_k109*na10));   // 341
   rDY[54] = ((((calc_k310*na3)+(calc_k910*na9))-((calc_k1011+calc_k103+calc_k109)*na10))+(calc_k1110*na11));  // 342
   rDY[55] = ((((calc_k411*na4)+(calc_k1011*na10))-((calc_k1110+calc_k114+calc_k1112)*na11))+(calc_k1211*na12));  // 343
   rDY[56] = ((((calc_k512*na5)+(calc_k1112*na11))-((calc_k1211+calc_k125+calc_k1213)*na12))+(calc_k1312*na13));  // 344
   rDY[57] = (((calc_k613*na6)+(calc_k1213*na12))-((calc_k1312+calc_k136)*na13));  // 345
   rDY[58] = ((calc_C2H_to_C1H*C2Herg)-(calc_C1H_to_C2H*C1Herg));  // 356
   rDY[59] = (calc_a1_C2_duplicated_IKr-calc_a2_C2_duplicated_IKr); // 359
   rDY[60] = (calc_a1_C3_duplicated_IKr-calc_a2_C3_duplicated_IKr); // 362
   rDY[61] = (calc_a1_O_duplicated_IKr-calc_a2_O_duplicated_IKr);   // 365
   rDY[62] = (calc_a1_I-calc_a2_I);   // 368
   rDY[63] = (((-calc_C0ks_C1ks)*C0ks)+(calc_C1ks_C0ks*C1ks));  // 375
   rDY[64] = (((calc_C0ks_C1ks*C0ks)-((calc_C1ks_C0ks+calc_C1ks_O1ks)*C1ks))+(calc_O1ks_C1ks*O1ks)); // 376
   rDY[65] = (((calc_C1ks_O1ks*C1ks)-((calc_O1ks_C1ks+calc_O1ks_O2ks)*O1ks))+(calc_O2ks_O1ks*O2ks)); // 377
   rDY[66] = ((calc_O1ks_O2ks*O1ks)-(calc_O2ks_O1ks*O2ks));  // 378
}
