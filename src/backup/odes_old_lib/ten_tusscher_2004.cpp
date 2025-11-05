/*
 * ten_Tusscher_2004.cpp
 *
 *  Created on: 27/10/2009
 *
 */

#include "ten_Tusscher_2004.h"

int cell_type = 5;
// Setting initial conditions to the ten Tusscher 2004 model
void setIC_ode_cpu(Real *sv, int cellID, int NEQ = NUMBER_EQUATIONS_CELL_MODEL)
{
   // Vector containing the initial conditions of the cell (EPI)
   Real IC[NEQ];

   // ENDO
   if(cell_type == 3) {
      // svolt; Units: millivolt; Initial value: -86.042839
      IC[ 0] = -86.042839;
      // sm; Units: dimensionless; Initial value: 0.766497
      IC[ 1] =  0.766497;
      // sh; Units: dimensionless; Initial value: 0.766376
      IC[ 2] =  0.766376;
      // sj; Units: dimensionless; Initial value: 0.766306
      IC[ 3] =  0.766306;
      // sxr1; Units: dimensionless; Initial value: 0.000188
      IC[ 4] =  0.000188;
      // sxr2; Units: dimensionless; Initial value: 0.479623
      IC[ 5] =  0.479623;
      // sxs; Units: dimensionless; Initial value: 0.003052
      IC[ 6] =  0.003052;
      // ss; Units: dimensionless; Initial value: 0.663683
      IC[ 7] =  0.663683;
      // sr; Units: dimensionless; Initial value: 0.000000
      IC[ 8] =  0.000000;
      // sd; Units: dimensionless; Initial value: 0.000020
      IC[ 9] =  0.000020;
      // sf; Units: dimensionless; Initial value: 0.999682
      IC[10] =  0.999682;
      // sfca; Units: dimensionless; Initial value: 1.006651
      IC[11] =  1.006651;
      // sg; Units: dimensionless; Initial value: 0.999979
      IC[12] =  0.999979;
      // Cai; Units: millimolar; Initial value: 0.000055
      IC[13] =  0.000055;
      // CaSR; Units: millimolar; Initial value: 0.345311
      IC[14] =  0.345311;
      // Nai; Units: millimolar; Initial value: 11.404968
      IC[15] =  11.404968;
      // Ki; Units: millimolar; Initial value: 138.300003
      IC[16] =  138.300003;
   }
   // M
   else if(cell_type == 4) {
      // svolt; Units: millivolt; Initial value: -85.990288
      IC[ 0] = -85.990288;
      // sm; Units: dimensionless; Initial value: 0.001459
      IC[ 1] =  0.001459;
      // sh; Units: dimensionless; Initial value: 0.765144
      IC[ 2] =  0.765144;
      // sj; Units: dimensionless; Initial value: 0.764952
      IC[ 3] =  0.764952;
      // sxr1; Units: dimensionless; Initial value: 0.000190
      IC[ 4] =  0.000190;
      // sxr2; Units: dimensionless; Initial value: 0.479077
      IC[ 5] =  0.479077;
      // sxs; Units: dimensionless; Initial value: 0.003064
      IC[ 6] =  0.003064;
      // ss; Units: dimensionless; Initial value: 0.999973
      IC[ 7] =  0.999973;
      // sr; Units: dimensionless; Initial value: 0.000000
      IC[ 8] =  0.000000;
      // sd; Units: dimensionless; Initial value: 0.000020
      IC[ 9] =  0.000020;
      // sf; Units: dimensionless; Initial value: 0.999681
      IC[10] =  0.999681;
      // sfca; Units: dimensionless; Initial value: 1.006487
      IC[11] =  1.006487;
      // sg; Units: dimensionless; Initial value: 0.999971
      IC[12] =  0.999971;
      // Cai; Units: millimolar; Initial value: 0.000059
      IC[13] =  0.000059;
      // CaSR; Units: millimolar; Initial value: 0.396961
      IC[14] =  0.396961;
      // Nai; Units: millimolar; Initial value: 11.322280
      IC[15] =  11.322280;
      // Ki; Units: millimolar; Initial value: 138.300003
      IC[16] =  138.300003;
   }
   // EPIC (default)
   else {
      // svolt; Units: millivolt; Initial value: -86.038116
      IC[ 0] = -86.038116;
      // sm; Units: dimensionless; Initial value: 0.001445
      IC[ 1] =  0.001445;
      // sh; Units: dimensionless; Initial value: 0.766376
      IC[ 2] =  0.766376;
      // sj; Units: dimensionless; Initial value: 0.766184
      IC[ 3] =  0.766184;
      // sxr1; Units: dimensionless; Initial value: 0.000189
      IC[ 4] =  0.000189;
      // sxr2; Units: dimensionless; Initial value: 0.479574
      IC[ 5] =  0.479574;
      // sxs; Units: dimensionless; Initial value: 0.003053
      IC[ 6] =  0.003053;
      // ss; Units: dimensionless; Initial value: 0.999973
      IC[ 7] =  0.999973;
      // sr; Units: dimensionless; Initial value: 0.000000
      IC[ 8] =  0.000000;
      // sd; Units: dimensionless; Initial value: 0.000020
      IC[ 9] =  0.000020;
      // sf; Units: dimensionless; Initial value: 0.999682
      IC[10] =  0.999682;
      // sfca; Units: dimensionless; Initial value: 1.006624
      IC[11] =  1.006624;
      // sg; Units: dimensionless; Initial value: 0.999978
      IC[12] =  0.999978;
      // Cai; Units: millimolar; Initial value: 0.000056
      IC[13] =  0.000056;
      // CaSR; Units: millimolar; Initial value: 0.353775
      IC[14] =  0.353775;
      // Nai; Units: millimolar; Initial value: 11.413716
      IC[15] =  11.413716;
      // Ki; Units: millimolar; Initial value: 138.300003
      IC[16] =  138.300003;
   }


   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = IC[i];
}

// Solving the ten_Tusscher model for each cell in the tissue matrix ni x nj
void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, int NEQ, int col)
{
   solve_Rush_Larsen_cpu(time, sv, dt, NEQ, cellID, col);
}

// Solving the ten Tusscher 2004 model using the Rush-Larsen method
void solve_Rush_Larsen_cpu(Real time, Real *sv, Real dt, int NEQ, int cellID, int col)
{
   // sv contains the membrane state variables (NEQ) of all cells in the 'tissue' (ni*nj)
   // cellID defines where the state variable starts in sv

   Real rY[NEQ], rDY[NEQ];

   for(int i = 0; i < NEQ; i++)
      rY[i] = sv[cellID+i];

   Real svolt  = rY[ 0];
   Real sm     = rY[ 1];
   Real sh     = rY[ 2];
   Real sj     = rY[ 3];
   Real sxr1   = rY[ 4];
   Real sxr2   = rY[ 5];
   Real sxs    = rY[ 6];
   Real ss     = rY[ 7];
   Real sr     = rY[ 8];
   Real sd     = rY[ 9];
   Real sf     = rY[10];
   Real sfca   = rY[11];
   Real sg     = rY[12];
   Real Cai    = rY[13];
   Real CaSR   = rY[14];
   Real Nai    = rY[15];
   Real Ki     = rY[16];

   // Stimulus current
   //const real stim_amplitude = -52.0;
   const Real stim_amplitude =  0.0;
   const Real stim_start     =  100.0 + (Real)col;
   const Real stim_duration  =  1.0;
   const Real stim_period    =  1000.0;
   const Real stim_end       =  100000.0;
   Real Istim = ((time >= stim_start) && (time <= stim_end) && (((time - stim_start) - (floor((time - stim_start) / stim_period) * stim_period)) <= stim_duration)) ? stim_amplitude : 0.0;

   //External concentrations
   Real Ko  = 5.4;
   Real Cao = 2.0;
   Real Nao = 140.0;

   //Intracellular volumes
   Real Vc  = 0.016404;
   Real Vsr = 0.001094;

   //Calcium dynamics
   Real Bufc   = 0.15;
   Real Kbufc  = 0.001;
   Real Bufsr  = 10.;
   Real Kbufsr = 0.3;
   Real taufca = 2.;
   Real taug   = 2.;
   Real Vmaxup = 0.000425;
   Real Kup    = 0.00025;

   //Constants
   Real R = 8314.472;
   Real F = 96485.3415;
   Real T = 310.0;
   Real RTONF;

   //Cellular capacitance
   Real CAPACITANCE = 0.185;

   //Parameters for currents

   //Parameters for IKr
   Real Gkr  = 0.096;

   //Parameters for Iks
   Real pKNa = 0.03;
   Real Gks;

   //Parameters for Ik1
   Real GK1 = 0.5*5.405;

   //Parameters for Ito
   Real Gto;

   //Parameters for INa
   Real GNa = 14.838;

   //Parameters for IbNa
   Real GbNa = 0.00029;

   //Parameters for INaK
   Real KmK  = 1.0;
   Real KmNa = 40.0;
   Real knak = 1.362;

   //Parameters for ICaL
   Real GCaL = 0.000175;

   //Parameters for IbCa
   Real GbCa = 0.000592;

   //Parameters for INaCa
   Real knaca = 1000;
   Real KmNai = 87.5;
   Real KmCa  = 1.38;
   Real ksat  = 0.1;
   Real n     = 0.35;

   //Parameters for IpCa
   Real GpCa = 0.825;
   Real KpCa = 0.0005;

   //Parameters for IpK;
   Real GpK = 0.0146;

   Real IKr;
   Real IKs;
   Real IK1;
   Real Ito;
   Real INa;
   Real IbNa;
   Real ICaL;
   Real IbCa;
   Real INaCa;
   Real IpCa;
   Real IpK;
   Real INaK;
   Real Irel;
   Real Ileak;
   Real sItot;

   Real dNai;
   Real dKi;
   Real dCai;
   Real dCaSR;

   Real A;
   Real SERCA;
   Real Caisquare;
   Real CaSRsquare;
   Real CaCurrent;
   Real CaSRCurrent;

   Real fcaold;
   Real gold;
   Real Ek;
   Real Ena;
   Real Eks;
   Real Eca;
   Real CaCSQN;
   Real bjsr;
   Real cjsr;
   Real CaBuf;
   Real bc;
   Real cc;
   Real Ak1;
   Real Bk1;
   Real rec_iK1;
   Real rec_ipK;
   Real rec_iNaK;
   Real AM;
   Real BM;
   Real AH_1;
   Real BH_1;
   Real AH_2;
   Real BH_2;
   Real AJ_1;
   Real BJ_1;
   Real AJ_2;
   Real BJ_2;
   Real M_INF;
   Real H_INF;
   Real J_INF;
   Real TAU_M;
   Real TAU_H;
   Real TAU_J;
   Real axr1;
   Real bxr1;
   Real axr2;
   Real bxr2;
   Real Xr1_INF;
   Real Xr2_INF;
   Real TAU_Xr1;
   Real TAU_Xr2;
   Real Axs;
   Real Bxs;
   Real Xs_INF;
   Real TAU_Xs;
   Real R_INF;
   Real TAU_R;
   Real S_INF;
   Real TAU_S;
   Real Ad;
   Real Bd;
   Real Cd;
   Real TAU_D;
   Real D_INF;
   Real TAU_F;
   Real F_INF;
   Real FCa_INF;
   Real G_INF;

   Real inverseVcF2;
   Real inverseVcF;
   Real Kupsquare;
   Real exptaufca;
   Real exptaug;

   // if endo cell
   if(cell_type == 3) {
      Gks = 0.245;
      Gto = 0.073;
   }
   // if M cell
   else if(cell_type == 4) {
      Gks = 0.062;
      Gto = 0.294;
   }
   // if epi cell
   else {
      Gks = 0.245;
      Gto = 0.294;
   }

   RTONF = ((R*T)/F);
   inverseVcF2  = 1./(2.*Vc*F);
   inverseVcF   = 1./(Vc*F);
   Kupsquare    = Kup*Kup;
   exptaufca    = exp(-dt/taufca);
   exptaug      = exp(-dt/taug);

   //Needed to compute currents
   Ek  = RTONF*(log((Ko/Ki)));
   Ena = RTONF*(log((Nao/Nai)));
   Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
   Eca = 0.5*RTONF*(log((Cao/Cai)));
   Ak1 = 0.1/(1.+exp(0.06*(svolt-Ek-200.)));
   Bk1 = (3.*exp(0.0002*(svolt-Ek+100.))+exp(0.1*(svolt-Ek-10.)))/(1.+exp(-0.5*(svolt-Ek)));
   rec_iK1  = Ak1/(Ak1+Bk1);
   rec_iNaK = (1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
   rec_ipK  = 1./(1.+exp((25-svolt)/5.98));

   //Compute currents
   INa   = GNa*sm*sm*sm*sh*sj*(svolt-Ena);
   ICaL  = (svolt != 0.0) ? (GCaL*sd*sf*sfca*4.*svolt*(F*F/(R*T))*(exp(2.*svolt*F/(R*T))*Cai-0.341*Cao)/(exp(2.*svolt*F/(R*T))-1.)) : (GCaL*sd*sf*sfca*4.*(svolt+1.0e-06)*(F*F/(R*T))*(exp(2.*(svolt+1.0e-06)*F/(R*T))*Cai-0.341*Cao)/(exp(2.*(svolt+1.0e-06)*F/(R*T))-1.));
   Ito   = Gto*sr*ss*(svolt-Ek);
   IKr   = Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
   IKs   = Gks*sxs*sxs*(svolt-Eks);
   IK1   = GK1*rec_iK1*(svolt-Ek);
   INaCa = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*(1./(1.+ksat*exp((n-1)*svolt*F/(R*T))))*(exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-exp((n-1.)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
   INaK  = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
   IpCa  = GpCa*Cai/(KpCa+Cai);
   IpK   = GpK*rec_ipK*(svolt-Ek);
   IbNa  = GbNa*(svolt-Ena);
   IbCa  = GbCa*(svolt-Eca);

   //Determine total current
   (sItot) = IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + Istim;

   //update concentrations
   Caisquare   =  Cai*Cai;
   CaSRsquare  =  CaSR*CaSR;
   CaCurrent   =- (ICaL+IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE;
   A           = 0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
   Irel        = A*sd*sg;
   Ileak       = 0.00008*(CaSR-Cai);
   SERCA       = Vmaxup/(1.+(Kupsquare/Caisquare));
   CaSRCurrent = SERCA-Irel-Ileak;
   CaCSQN      = Bufsr*CaSR/(CaSR+Kbufsr);
   dCaSR       = dt*(Vc/Vsr)*CaSRCurrent;
   bjsr        = Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
   cjsr        = Kbufsr*(CaCSQN+dCaSR+CaSR);
   CaSR        = (sqrt(bjsr*bjsr+4.*cjsr)-bjsr)/2.;
   CaBuf       = Bufc*Cai/(Cai+Kbufc);
   dCai        = dt*(CaCurrent-CaSRCurrent);
   bc          = Bufc-CaBuf-dCai-Cai+Kbufc;
   cc          = Kbufc*(CaBuf+dCai+Cai);
   Cai         = (sqrt(bc*bc+4.*cc)-bc)/2.;
   dNai        =- (INa+IbNa+3.*INaK+3.*INaCa)*inverseVcF*CAPACITANCE;
   Nai         += dt*dNai;
   dKi         =- (IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
   Ki          += dt*dKi;

   //compute steady state values and time constants
   AM    = 1./(1.+exp((-60.-svolt)/5.));
   BM    = 0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
   TAU_M = AM*BM;
   M_INF = 1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
   if(svolt >= -40.) {
      AH_1  = 0.;
      BH_1  = (0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
      TAU_H = 1.0/(AH_1+BH_1);
   }
   else {
      AH_2  = (0.057*exp(-(svolt+80.)/6.8));
      BH_2  = (2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
      TAU_H = 1.0/(AH_2+BH_2);
   }
   H_INF = 1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
   if(svolt >= -40.) {
      AJ_1  = 0.;
      BJ_1  = (0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
      TAU_J = 1.0/(AJ_1+BJ_1);
   }
   else {
      AJ_2  = (((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*exp(-0.04391*svolt))*(svolt+37.78)/(1.+exp(0.311*(svolt+79.23))));
      BJ_2  = (0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
      TAU_J = 1.0/(AJ_2+BJ_2);
   }
   J_INF=H_INF;

   Xr1_INF = 1./(1.+exp((-26.-svolt)/7.));
   axr1    = 450./(1.+exp((-45.-svolt)/10.));
   bxr1    = 6./(1.+exp((svolt-(-30.))/11.5));
   TAU_Xr1 = axr1*bxr1;
   Xr2_INF = 1./(1.+exp((svolt-(-88.))/24.));
   axr2    = 3./(1.+exp((-60.-svolt)/20.));
   bxr2    = 1.12/(1.+exp((svolt-60.)/20.));
   TAU_Xr2 = axr2*bxr2;

   Xs_INF = 1./(1.+exp((-5.-svolt)/14.));
   Axs    = 1100./(sqrt(1.+exp((-10.-svolt)/6)));
   Bxs    = 1./(1.+exp((svolt-60.)/20.));
   TAU_Xs = Axs*Bxs;

   // if endo cell
   if(cell_type == 3) {
      R_INF = 1./(1.+exp((20-svolt)/6.));
      S_INF = 1./(1.+exp((svolt+28)/5.));
      TAU_R = 9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
      TAU_S = 1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
   }
   // if M cell
   else if(cell_type == 4) {
      R_INF = 1./(1.+exp((20-svolt)/6.));
      S_INF = 1./(1.+exp((svolt+20)/5.));
      TAU_R = 9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
      TAU_S = 85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
   }
   // if epi cell
   else {
      R_INF = 1./(1.+exp((20.-svolt)/6.));
      S_INF = 1./(1.+exp((svolt+20.)/5.));
      TAU_R = 9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
      TAU_S = 85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
   }

   D_INF = 1.0/(1.+exp((-5.-svolt)/7.5));
   Ad    = 1.4/(1.+exp((-35.-svolt)/13.))+0.25;
   Bd    = 1.4/(1.+exp((svolt+5.)/5.));
   Cd    = 1.0/(1.+exp((50.-svolt)/20.));
   TAU_D = Ad*Bd+Cd;
   F_INF = 1./(1.+exp((svolt+20.)/7.));
   TAU_F = 1125.*exp(-(svolt+27.)*(svolt+27.)/240.)+80.+165./(1.+exp((25.-svolt)/10.));

   FCa_INF=(1./(1.+pow((Cai/0.000325),8.))+0.1/(1.+exp((Cai-0.0005)/0.0001))+0.20/(1.+exp((Cai-0.00075)/0.0008))+0.23 )/1.46;
   if(Cai < 0.00035)
      G_INF = 1./(1.+pow((Cai/0.00035),6.));
   else
      G_INF = 1./(1.+pow((Cai/0.00035),16.));

   //Update gates
   sm     = M_INF-(M_INF-sm)*exp(-dt/TAU_M);
   sh     = H_INF-(H_INF-sh)*exp(-dt/TAU_H);
   sj     = J_INF-(J_INF-sj)*exp(-dt/TAU_J);
   sxr1   = Xr1_INF-(Xr1_INF-sxr1)*exp(-dt/TAU_Xr1);
   sxr2   = Xr2_INF-(Xr2_INF-sxr2)*exp(-dt/TAU_Xr2);
   sxs    = Xs_INF-(Xs_INF-sxs)*exp(-dt/TAU_Xs);
   ss     = S_INF-(S_INF-ss)*exp(-dt/TAU_S);
   sr     = R_INF-(R_INF-sr)*exp(-dt/TAU_R);
   sd     = D_INF-(D_INF-sd)*exp(-dt/TAU_D);
   sf     = F_INF-(F_INF-sf)*exp(-dt/TAU_F);
   fcaold = sfca;
   sfca   = FCa_INF-(FCa_INF-sfca)*exptaufca;
   if(sfca > fcaold && (svolt) > -60.)
      sfca = fcaold;
   gold = sg;
   sg = G_INF-(G_INF-sg)*exptaug;
   if(sg > gold && (svolt) > -60.)
      sg = gold;

   //update voltage
   svolt = svolt + dt*(-sItot);

   rDY[ 0] = svolt;
   rDY[ 1] = sm;
   rDY[ 2] = sh;
   rDY[ 3] = sj;
   rDY[ 4] = sxr1;
   rDY[ 5] = sxr2;
   rDY[ 6] = sxs;
   rDY[ 7] = ss;
   rDY[ 8] = sr;
   rDY[ 9] = sd;
   rDY[10] = sf;
   rDY[11] = sfca;
   rDY[12] = sg;
   rDY[13] = Cai;
   rDY[14] = CaSR;
   rDY[15] = Nai;
   rDY[16] = Ki;

   for(int i = 0; i < NEQ; i++)
      sv[cellID+i] = rDY[i];
}
