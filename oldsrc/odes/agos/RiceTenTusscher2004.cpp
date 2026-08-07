#include "RiceTenTusscher2004.hpp"

#define AGOS_NAN NAN
#define AGOS_INF INFINITY
#define __agos_xor(a,b) (!(a && b) && (a || b))
float __agos_factorial(int);

Solveode::Solveode(int type)
{
  if(type != 3 && type != 4 && type != 5)
  {
    std::cerr << "ERRO RICE TENTUSSCHER - type = " << type << std::endl;
    exit(1);
  }

  cell_type = type;

  std::cout << "Cell type = " << type << " ";

  time = 0.000000000000000e+00;
  R = 8.314472000000000e+03;
  T = 3.100000000000000e+02;
  F = 9.648534149999999e+04;
  Na_o = 1.400000000000000e+02;
  K_o = 5.400000000000000e+00;
  P_kna = 3.000000000000000e-02;
  Ca_o = 2.000000000000000e+00;
  g_K1 = 5.405000000000000e+00;
  g_Kr = 9.600000000000000e-02;
  if(cell_type ==4){
    g_Ks = 0.062;
  }
  else{
    g_Ks = 2.450000000000000e-01;
  }
  g_Na = 1.483800000000000e+01;
  g_bna = 2.900000000000000e-04;

  i_stimAmplitude = 0.0;
  i_stimStart = 0.0;
  i_stimEnd = 0.0;
  i_stimPeriod = 0.0;
  i_stimPulseDuration = 0.0;

  /*IstimAmplitude = -5.20000000e+01;
  IstimStart = 1.00000000e-04;
  IstimEnd = 5.00000000e+06;
  IstimPeriod = 1.00000000e+03;
  IstimPulseDuration = 1.00000000e+00;
  */

  /*
        i_stimAmplitude = -5.200000000000000e+01;
        i_stimStart = 0.000000000000000e+00;
        i_stimEnd = 5.000000000000000e+04;
        i_stimPeriod = 1.000000000000000e+03;
        i_stimPulseDuration = 1.000000000000000e+00;
   */

  g_CaL = 1.750000000000000e-04;
  g_bca = 5.920000000000000e-04;
  if(cell_type !=3){
    g_to = 0.294;
  }
  else{
    g_to = 7.300000000000000e-02;
  }
  P_NaK = 1.362000000000000e+00;
  K_mk = 1.000000000000000e+00;
  K_mNa = 4.000000000000000e+01;
  K_NaCa = 1.000000000000000e+03;
  gamma = 3.500000000000000e-01;
  alpha = 2.500000000000000e+00;
  Km_Nai = 8.750000000000000e+01;
  Km_Ca = 1.380000000000000e+00;
  K_sat = 1.000000000000000e-01;
  g_pCa = 8.250000000000000e-01;
  K_pCa = 5.000000000000000e-04;
  g_pK = 1.460000000000000e-02;
  a_rel = 1.646400000000000e-02;
  b_rel = 2.500000000000000e-01;
  c_rel = 8.232000000000000e-03;
  Vmax_up = 4.250000000000000e-04;
  K_up = 2.500000000000000e-04;
  V_leak = 8.000000000000001e-05;
  tau_g = 2.000000000000000e+00;
  Buf_sr = 1.000000000000000e+01;
  K_buf_sr = 3.000000000000000e-01;
  V_c = 1.640400000000000e-02;
  V_sr = 1.094000000000000e-03;
  Cm = 1.850000000000000e-01;
  K_buf_c = 1.000000000000000e-03;
  Buf_c = 0.131477;
  len_thick = 1.650000000000000e+00;
  len_thin = 1.200000000000000e+00;
  len_hbare = 1.000000000000000e-01;
  kon = 5.000000000000000e-02 * 0.992879;
  Qkon = 1.500000000000000e+00;
  TmpC = 3.700000000000000e+01;
  koffL = 2.500000000000000e-01;
  koffmod = 1.00634;
  Qkoff = 1.300000000000000e+00;
  koffH = 2.500000000000000e-02;
  kn_p = 5.000000000000000e-01 * 1.21556;
  Qkn_p = 1.600000000000000e+00;
  kp_n = 5.000000000000000e-02 * 0.320597;
  Qkp_n = 1.600000000000000e+00;
  perm50 = 5.000000000000000e-01;
  nperm = 1.500000000000000e+01;
  fapp = 5.000000000000000e-01 * 9.61566;
  xbmodsp = 0.20000000000000e+00;
  Qfapp = 6.250000000000000e+00;
  gslmod = 6.000000000000000e+00;
  gapp = 7.000000000000001e-02 * 1.3349;
  Qgapp = 2.500000000000000e+00;
  hfmdc = 5.000000000000000e+00;
  x_0 = 7.000000000000000e-03;
  hbmdc = 0.000000000000000e+00;
  hf = 2.000000000000000e+00 * 0.00513577;
  Qhf = 6.250000000000000e+00;
  hb = 4.000000000000000e-01 * 0.0850933;
  Qhb = 6.250000000000000e+00;
  sigmap = 8.000000000000000e+00;
  sigman = 1.000000000000000e+00;
  gxb = 7.000000000000001e-02 * 0.432424;
  Qgxb = 6.250000000000000e+00;
  xPsi = 2.000000000000000e+00;
  kxb = 1.200000000000000e+02;
  SLrest = 1.850000000000000e+00;
  PCon_t = 2.000000000000000e-03;
  PExp_t = 1.000000000000000e+01;
  PCon_c = 2.000000000000000e-02;
  PExp_c = 7.000000000000000e+01;
  SL_c = 2.250000000000000e+00;
  SLset = 1.900000000000000e+00;
  KSE = 1.000000000000000e+00;
  SEon = 1.000000000000000e+00;
  visc = 3.000000000000000e+00;
  massf = 5.000000000000000e+01;
  SLmax = 2.400000000000000e+00;
  SLmin = 1.400000000000000e+00;
  Trop_conc = 7.000000000000000e+01;
  stretch = 1.0;
  dtime = 0.0; time_vec__ = NULL;

  /*
        dependent_variable__ = N_VNew_Serial(28);
        if(check_flag((void *)dependent_variable__, "N_VNew_Serial", 0))
            exit(1);
        abstol__ = N_VNew_Serial(28);
        if (check_flag((void *)abstol__, "N_VNew_Serial", 0))
            exit(1);

        depvar__ = (double*)malloc(sizeof(double)*28);
        if(depvar__ == NULL){
            fprintf(stderr, "ERROR Cannot allocate memory for depvar__\n");
            exit(0);
        }
   */

   Xr1 = NULL;
   /*NV_Ith_S(dependent_variable__, 0) =*/
   Xr1_ini_ = 1.817573054244821e-04;
   //NV_Ith_S(abstol__, 0) = val_abstol__;
   Xr2 = NULL;
   /*NV_Ith_S(dependent_variable__, 1) =*/
   Xr2_ini_ = 4.822125371672458e-01;
   //NV_Ith_S(abstol__, 1) = val_abstol__;
   Xs = NULL;
   /*NV_Ith_S(dependent_variable__, 2) =*/
   Xs_ini_ = 2.998779317686131e-03;
   //NV_Ith_S(abstol__, 2) = val_abstol__;
   m = NULL;
   /*NV_Ith_S(dependent_variable__, 3) =*/
   m_ini_ = 1.368501861611859e-03;
   //NV_Ith_S(abstol__, 3) = val_abstol__;
   h = NULL;
   /*NV_Ith_S(dependent_variable__, 4) =*/
   h_ini_ = 7.728393292625330e-01;
   //NV_Ith_S(abstol__, 4) = val_abstol__;
   j = NULL;
   /*NV_Ith_S(dependent_variable__, 5) =*/
   j_ini_ = 3.694904454167577e-01;
   //NV_Ith_S(abstol__, 5) = val_abstol__;
   V = NULL;
   /*NV_Ith_S(dependent_variable__, 6) =*/
   V_ini_ = -8.629169987876476e+01;
   //NV_Ith_S(abstol__, 6) = val_abstol__;
   d = NULL;
   /*NV_Ith_S(dependent_variable__, 7) =*/
   d_ini_ = 1.962096961554300e-05;
   //NV_Ith_S(abstol__, 7) = val_abstol__;
   f = NULL;
   /*NV_Ith_S(dependent_variable__, 8) =*/
   f_ini_ = 9.998430021390249e-01;
   //NV_Ith_S(abstol__, 8) = val_abstol__;
   fCa = NULL;
   /*NV_Ith_S(dependent_variable__, 9) =*/
   fCa_ini_ = 1.006126715016802e+00;
   //NV_Ith_S(abstol__, 9) = val_abstol__;
   s = NULL;
   /*NV_Ith_S(dependent_variable__, 10) =*/
   s_ini_ = 6.699733924504092e-01;
   //NV_Ith_S(abstol__, 10) = val_abstol__;
   r = NULL;
   /*NV_Ith_S(dependent_variable__, 11) =*/
   r_ini_ = 2.024681561943644e-08;
   //NV_Ith_S(abstol__, 11) = val_abstol__;
   g = NULL;
   /*NV_Ith_S(dependent_variable__, 12) =*/
   g_ini_ = 9.999492474604431e-01;
   //NV_Ith_S(abstol__, 12) = val_abstol__;
   Ca_SR = NULL;
   /*NV_Ith_S(dependent_variable__, 13) =*/
   Ca_SR_ini_ =  5.205023026469432e-01;
   //NV_Ith_S(abstol__, 13) = val_abstol__;
   Ca_tot = NULL;
   /*NV_Ith_S(dependent_variable__, 14) =*/
   Ca_tot_ini_ = 9.281944866155925e-03;
   //NV_Ith_S(abstol__, 14) = val_abstol__;
   Na_i = NULL;
   /*NV_Ith_S(dependent_variable__, 15) =*/
   Na_i_ini_ = 1.114903806343864e+01;
   //NV_Ith_S(abstol__, 15) = val_abstol__;
   K_i = NULL;
   /*NV_Ith_S(dependent_variable__, 16) =*/
   K_i_ini_ = 1.383870498592104e+02;
   //NV_Ith_S(abstol__, 16) = val_abstol__;
   TRPNCaL = NULL;
   /*NV_Ith_S(dependent_variable__, 17) =*/
   TRPNCaL_ini_ =  1.311635148577327e-02;
   //NV_Ith_S(abstol__, 17) = val_abstol__;
   TRPNCaH = NULL;
   /*NV_Ith_S(dependent_variable__, 18) =*/
   TRPNCaH_ini_ =  1.175428809165677e-01;
   //NV_Ith_S(abstol__, 18) = val_abstol__;
   N_NoXB = NULL;
   /*NV_Ith_S(dependent_variable__, 19) =*/
   N_NoXB_ini_ =  9.999995709556512e-01;
   //NV_Ith_S(abstol__, 19) = val_abstol__;
   P_NoXB = NULL;
   /*NV_Ith_S(dependent_variable__, 20) =*/
   P_NoXB_ini_ =  4.290431881685032e-07;
   //NV_Ith_S(abstol__, 20) = val_abstol__;
   N = NULL;
   /*NV_Ith_S(dependent_variable__, 21) =*/
   N_ini_ =  9.996687602297853e-01;
   //NV_Ith_S(abstol__, 21) = val_abstol__;
   XBprer = NULL;
   /*NV_Ith_S(dependent_variable__, 22) =*/
   XBprer_ini_ =  2.211689363408625e-04;
   //NV_Ith_S(abstol__, 22) = val_abstol__;
   XBpostr = NULL;
   /*NV_Ith_S(dependent_variable__, 23) =*/
   XBpostr_ini_ =  1.079302937362542e-04;
   //NV_Ith_S(abstol__, 23) = val_abstol__;
   xXBprer = NULL;
   /*NV_Ith_S(dependent_variable__, 24) =*/
   xXBprer_ini_ =  1.007048444140799e-18;
   //NV_Ith_S(abstol__, 24) = val_abstol__;
   xXBpostr = NULL;
   /*NV_Ith_S(dependent_variable__, 25) =*/
   xXBpostr_ini_ =  7.000000000000143e-03;
   //NV_Ith_S(abstol__, 25) = val_abstol__;
   SL = NULL;
   /*NV_Ith_S(dependent_variable__, 26) =*/
   SL_ini_ = 1.9;
   //NV_Ith_S(abstol__, 26) = val_abstol__;
   intf = NULL;
   /*NV_Ith_S(dependent_variable__, 27) =*/
   intf_ini_ = -3.637088203841350e+07;
   //NV_Ith_S(abstol__, 27) = val_abstol__;

   Ca_i=NULL;
   active_force=NULL;
   //      reltol__ = val_reltol__;
   it_countx = 0;

   //      int nonlineariteration__;
   //if(method__ == CV_BDF) nonlineariteration__ = CV_NEWTON;
   //else nonlineariteration__ = CV_FUNCTIONAL;
   //cvode_mem_cvode__ = CVodeCreate(method__, nonlineariteration__);
   //if (check_flag((void *)cvode_mem_cvode__, "CVodeCreate", 0))
   //  exit(1);

   //flag__ = CVodeMalloc(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SV, reltol__, abstol__);
   //if (check_flag(&flag__, "CVodeMalloc", 1))
   //  exit(1);

   //flag__ = CVDense(cvode_mem_cvode__, 28);
   //if (check_flag(&flag__, "CVDense", 1))
   //  exit(1);
   //CVodeSetFdata(cvode_mem_cvode__, (void*)this);
}
Solveode::~Solveode()
{
  //if(depvar__ != NULL) free(depvar__);
  if(Xr1 != NULL) free(Xr1);
  if(Xr2 != NULL) free(Xr2);
  if(Xs != NULL) free(Xs);
  if(m != NULL) free(m);
  if(h != NULL) free(h);
  if(j != NULL) free(j);
  if(V != NULL) free(V);
  if(d != NULL) free(d);
  if(f != NULL) free(f);
  if(fCa != NULL) free(fCa);
  if(s != NULL) free(s);
  if(r != NULL) free(r);
  if(g != NULL) free(g);
  if(Ca_SR != NULL) free(Ca_SR);
  if(Ca_tot != NULL) free(Ca_tot);
  if(Na_i != NULL) free(Na_i);
  if(K_i != NULL) free(K_i);
  if(TRPNCaL != NULL) free(TRPNCaL);
  if(TRPNCaH != NULL) free(TRPNCaH);
  if(N_NoXB != NULL) free(N_NoXB);
  if(P_NoXB != NULL) free(P_NoXB);
  if(N != NULL) free(N);
  if(XBprer != NULL) free(XBprer);
  if(XBpostr != NULL) free(XBpostr);
  if(xXBprer != NULL) free(xXBprer);
  if(xXBpostr != NULL) free(xXBpostr);
  if(SL != NULL) free(SL);
  if(intf != NULL) free(intf);
  if(Ca_i != NULL) free(Ca_i);
  if(active_force != NULL) free(active_force);
  if(time_vec__ != NULL) free(time_vec__);

  //N_VDestroy_Serial(dependent_variable__);
  //N_VDestroy_Serial(abstol__);
  //CVodeFree(&cvode_mem_cvode__);
}

int Solveode::setVariables(int indVariable, double value_new)
{
  switch(indVariable)
  {
  case 0:     /*NV_Ith_S(dependent_variable__, 0) =*/ Xr1_old_ = Xr1_ini_ = value_new;    break;
  case 1:     /*NV_Ith_S(dependent_variable__, 1) =*/ Xr2_old_ = Xr2_ini_ = value_new;    break;
  case 2:     /*NV_Ith_S(dependent_variable__, 2) =*/ Xs_old_ = Xs_ini_ = value_new;    break;
  case 3:     /*NV_Ith_S(dependent_variable__, 3) =*/ m_old_ = m_ini_ = value_new;    break;
  case 4:     /*NV_Ith_S(dependent_variable__, 4) =*/ h_old_ = h_ini_ = value_new;    break;
  case 5:     /*NV_Ith_S(dependent_variable__, 5) =*/ j_old_ = j_ini_ = value_new;    break;
  case 6:     /*NV_Ith_S(dependent_variable__, 6) =*/ V_old_ = V_ini_ = value_new;    break;
  case 7:     /*NV_Ith_S(dependent_variable__, 7) =*/ d_old_ = d_ini_ = value_new;    break;
  case 8:     /*NV_Ith_S(dependent_variable__, 8) =*/ f_old_ = f_ini_ = value_new;    break;
  case 9:     /*NV_Ith_S(dependent_variable__, 9) =*/ fCa_old_ = fCa_ini_ = value_new;    break;
  case 10:        /*NV_Ith_S(dependent_variable__, 10) =*/ s_old_ = s_ini_ = value_new;    break;
  case 11:        /*NV_Ith_S(dependent_variable__, 11) =*/ r_old_ = r_ini_ = value_new;    break;
  case 12:        /*NV_Ith_S(dependent_variable__, 12) =*/ g_old_ = g_ini_ = value_new;    break;
  case 13:        /*NV_Ith_S(dependent_variable__, 13) =*/ Ca_SR_old_ = Ca_SR_ini_ = value_new;    break;
  case 14:        /*NV_Ith_S(dependent_variable__, 14) =*/ Ca_tot_old_ = Ca_tot_ini_ = value_new;    break;
  case 15:        /*NV_Ith_S(dependent_variable__, 15) =*/ Na_i_old_ = Na_i_ini_ = value_new;    break;
  case 16:        /*NV_Ith_S(dependent_variable__, 16) =*/ K_i_old_ = K_i_ini_ = value_new;    break;
  case 17:        /*NV_Ith_S(dependent_variable__, 17) =*/ TRPNCaL_old_ = TRPNCaL_ini_ = value_new;    break;
  case 18:        /*NV_Ith_S(dependent_variable__, 18) =*/ TRPNCaH_old_ = TRPNCaH_ini_ = value_new;    break;
  case 19:        /*NV_Ith_S(dependent_variable__, 19) =*/ N_NoXB_old_ = N_NoXB_ini_ = value_new;    break;
  case 20:        /*NV_Ith_S(dependent_variable__, 20) =*/ P_NoXB_old_ = P_NoXB_ini_ = value_new;    break;
  case 21:        /*NV_Ith_S(dependent_variable__, 21) =*/ N_old_ = N_ini_ = value_new;    break;
  case 22:        /*NV_Ith_S(dependent_variable__, 22) =*/ XBprer_old_ = XBprer_ini_ = value_new;    break;
  case 23:        /*NV_Ith_S(dependent_variable__, 23) =*/ XBpostr_old_ = XBpostr_ini_ = value_new;    break;
  case 24:        /*NV_Ith_S(dependent_variable__, 24) =*/ xXBprer_old_ = xXBprer_ini_ = value_new;    break;
  case 25:        /*NV_Ith_S(dependent_variable__, 25) =*/ xXBpostr_old_ = xXBpostr_ini_ = value_new;    break;
  case 26:        /*NV_Ith_S(dependent_variable__, 26) =*/ SL_old_ = SL_ini_ = value_new;    break;
  case 27:        /*NV_Ith_S(dependent_variable__, 27) =*/ intf_old_ = intf_ini_ = value_new;    break;
  default:    return 1;    break;
  }
  return 0;
}

int Solveode::setParameters(int indVariable, double value_new)
{
  switch(indVariable)
  {
  case 0:  time = time_new = value_new;   break;
  case 1:  R = value_new;   break;
  case 2:  T = value_new;   break;
  case 3:  F = value_new;   break;
  case 4:  Na_o = value_new;   break;
  case 5:  K_o = value_new;   break;
  case 6:  P_kna = value_new;   break;
  case 7:  Ca_o = value_new;   break;
  case 8:  g_K1 = value_new;   break;
  case 9:  g_Kr = value_new;   break;
  case 10: g_Ks = value_new;   break;
  case 11: g_Na = value_new;   break;
  case 12: g_bna = value_new;   break;
  case 13: i_stimAmplitude = value_new;   break;
  case 14: i_stimStart = value_new;   break;
  case 15: i_stimEnd = value_new;   break;
  case 16: i_stimPeriod = value_new;   break;
  case 17: i_stimPulseDuration = value_new;   break;
  case 18: g_CaL = value_new;   break;
  case 19: g_bca = value_new;   break;
  case 20: g_to = value_new;   break;
  case 21: P_NaK = value_new;   break;
  case 22: K_mk = value_new;   break;
  case 23: K_mNa = value_new;   break;
  case 24: K_NaCa = value_new;   break;
  case 25: gamma = value_new;   break;
  case 26: alpha = value_new;   break;
  case 27: Km_Nai = value_new;   break;
  case 28: Km_Ca = value_new;   break;
  case 29: K_sat = value_new;   break;
  case 30: g_pCa = value_new;   break;
  case 31: K_pCa = value_new;   break;
  case 32: g_pK = value_new;   break;
  case 33: a_rel = value_new;   break;
  case 34: b_rel = value_new;   break;
  case 35: c_rel = value_new;   break;
  case 36: Vmax_up = value_new;   break;
  case 37: K_up = value_new;   break;
  case 38: V_leak = value_new;   break;
  case 39: tau_g = value_new;   break;
  case 40: Buf_sr = value_new;   break;
  case 41: K_buf_sr = value_new;   break;
  case 42: V_c = value_new;   break;
  case 43: V_sr = value_new;   break;
  case 44: Cm = value_new;   break;
  case 45: K_buf_c = value_new;   break;
  case 46: Buf_c = value_new;   break;
  case 47: len_thick = value_new;   break;
  case 48: len_thin = value_new;   break;
  case 49: len_hbare = value_new;   break;
  case 50:       kon = value_new;   break;
  case 51:       Qkon = value_new;   break;
  case 52:       TmpC = value_new;   break;
  case 53:       koffL = value_new;   break;
  case 54:        koffmod = value_new;   break;
  case 55:        Qkoff = value_new;   break;
  case 56:        koffH = value_new;   break;
  case 57:        kn_p = value_new;   break;
  case 58:        Qkn_p = value_new;   break;
  case 59:        kp_n = value_new;   break;
  case 60:        Qkp_n = value_new;   break;
  case 61:        perm50 = value_new;   break;
  case 62:        nperm = value_new;   break;
  case 63:        fapp = value_new;   break;
  case 64:        xbmodsp = value_new;   break;
  case 65:        Qfapp = value_new;   break;
  case 66:        gslmod = value_new;   break;
  case 67:        gapp = value_new;   break;
  case 68:        Qgapp = value_new;   break;
  case 69:        hfmdc = value_new;   break;
  case 70:        x_0 = value_new;   break;
  case 71:        hbmdc = value_new;   break;
  case 72:        hf = value_new;   break;
  case 73:        Qhf = value_new;   break;
  case 74:        hb = value_new;   break;
  case 75:       Qhb = value_new;   break;
  case 76:       sigmap = value_new;   break;
  case 77:       sigman = value_new;   break;
  case 78:       gxb = value_new;   break;
  case 79:       Qgxb = value_new;   break;
  case 80:       xPsi = value_new;   break;
  case 81:       kxb = value_new;   break;
  case 82:       SLrest = value_new;   break;
  case 83:       PCon_t = value_new;   break;
  case 84:       PExp_t = value_new;   break;
  case 85:       PCon_c = value_new;   break;
  case 86:       PExp_c = value_new;   break;
  case 87:       SL_c = value_new;   break;
  case 88:       SLset = value_new;   break;
  case 89:       KSE = value_new;   break;
  case 90:       SEon = value_new;   break;
  case 91:       visc = value_new;   break;
  case 92:        massf = value_new;   break;
  case 93:       SLmax = value_new;   break;
  case 94:       SLmin = value_new;   break;
  case 95:       Trop_conc = value_new;   break;
  case 96:       cell_type = (int) value_new; break;
  case 97:       stretch = value_new; break;
  default:    return 1;    break;
  }
  return 0;
}

int Solveode::setFreeVariable(double value_new)
{
  dtime = value_new;
}

//Get Methods

double Solveode::getVariables(int indVariable)
{
  switch(indVariable)
  {
  case 0:  return Xr1_ini_;    break;
  case 1:  return Xr2_ini_;    break;
  case 2:  return Xs_ini_;    break;
  case 3:  return m_ini_;    break;
  case 4:  return h_ini_;    break;
  case 5:  return j_ini_;    break;
  case 6:  return V_ini_;    break;
  case 7:  return d_ini_;    break;
  case 8:  return f_ini_;    break;
  case 9:  return fCa_ini_;    break;
  case 10: return s_ini_;    break;
  case 11: return r_ini_;    break;
  case 12: return g_ini_;    break;
  case 13: return Ca_SR_ini_;    break;
  case 14: return Ca_tot_ini_;    break;
  case 15: return Na_i_ini_;    break;
  case 16: return K_i_ini_;    break;
  case 17: return TRPNCaL_ini_;    break;
  case 18: return TRPNCaH_ini_;    break;
  case 19: return N_NoXB_ini_;    break;
  case 20: return P_NoXB_ini_;    break;
  case 21: return N_ini_;    break;
  case 22: return XBprer_ini_;    break;
  case 23: return XBpostr_ini_;    break;
  case 24: return xXBprer_ini_;    break;
  case 25: return xXBpostr_ini_;    break;
  case 26: return SL_ini_;    break;
  case 27: return intf_ini_;    break;
  default: return 1;    break;
  }
}

double Solveode::getParameters(int indVariable)
{
  switch(indVariable)
  {
  case 0:     return time;    break;
  case 1:     return R;    break;
  case 2:     return T;    break;
  case 3:     return F;    break;
  case 4:     return Na_o;    break;
  case 5:     return K_o;    break;
  case 6:     return P_kna;    break;
  case 7:     return Ca_o;    break;
  case 8:     return g_K1;    break;
  case 9:     return g_Kr;    break;
  case 10:        return g_Ks;    break;
  case 11:        return g_Na;    break;
  case 12:        return g_bna;    break;
  case 13:        return i_stimAmplitude;    break;
  case 14:        return i_stimStart;    break;
  case 15:        return i_stimEnd;    break;
  case 16:        return i_stimPeriod;    break;
  case 17:        return i_stimPulseDuration;    break;
  case 18:        return g_CaL;    break;
  case 19:        return g_bca;    break;
  case 20:        return g_to;    break;
  case 21:        return P_NaK;    break;
  case 22:        return K_mk;    break;
  case 23:        return K_mNa;    break;
  case 24:        return K_NaCa;    break;
  case 25:        return gamma;    break;
  case 26:        return alpha;    break;
  case 27:        return Km_Nai;    break;
  case 28:        return Km_Ca;    break;
  case 29:        return K_sat;    break;
  case 30:        return g_pCa;    break;
  case 31:        return K_pCa;    break;
  case 32:        return g_pK;    break;
  case 33:        return a_rel;    break;
  case 34:        return b_rel;    break;
  case 35:        return c_rel;    break;
  case 36:        return Vmax_up;    break;
  case 37:        return K_up;    break;
  case 38:        return V_leak;    break;
  case 39:        return tau_g;    break;
  case 40:        return Buf_sr;    break;
  case 41:        return K_buf_sr;    break;
  case 42:        return V_c;    break;
  case 43:        return V_sr;    break;
  case 44:        return Cm;    break;
  case 45:        return K_buf_c;    break;
  case 46:        return Buf_c;    break;
  case 47:        return len_thick;    break;
  case 48:        return len_thin;    break;
  case 49:        return len_hbare;    break;
  case 50:        return kon;    break;
  case 51:        return Qkon;    break;
  case 52:        return TmpC;    break;
  case 53:        return koffL;    break;
  case 54:        return koffmod;    break;
  case 55:        return Qkoff;    break;
  case 56:        return koffH;    break;
  case 57:        return kn_p;    break;
  case 58:        return Qkn_p;    break;
  case 59:        return kp_n;    break;
  case 60:        return Qkp_n;    break;
  case 61:        return perm50;    break;
  case 62:        return nperm;    break;
  case 63:        return fapp;    break;
  case 64:        return xbmodsp;    break;
  case 65:        return Qfapp;    break;
  case 66:        return gslmod;    break;
  case 67:        return gapp;    break;
  case 68:        return Qgapp;    break;
  case 69:        return hfmdc;    break;
  case 70:        return x_0;    break;
  case 71:        return hbmdc;    break;
  case 72:        return hf;    break;
  case 73:        return Qhf;    break;
  case 74:        return hb;    break;
  case 75:        return Qhb;    break;
  case 76:        return sigmap;    break;
  case 77:        return sigman;    break;
  case 78:        return gxb;    break;
  case 79:        return Qgxb;    break;
  case 80:        return xPsi;    break;
  case 81:        return kxb;    break;
  case 82:        return SLrest;    break;
  case 83:        return PCon_t;    break;
  case 84:        return PExp_t;    break;
  case 85:        return PCon_c;    break;
  case 86:        return PExp_c;    break;
  case 87:        return SL_c;    break;
  case 88:        return SLset;    break;
  case 89:        return KSE;    break;
  case 90:        return SEon;    break;
  case 91:        return visc;    break;
  case 92:        return massf;    break;
  case 93:        return SLmax;    break;
  case 94:        return SLmin;    break;
  case 95:        return Trop_conc;    break;
  case 96:        return cell_type;    break;
  case 97:        return stretch; break;
  default:    break;
  }
}

double Solveode::getFreeVariable()
{
  return dtime;
}

//Get Methods - Variables

Variables Solveode::get_Variables()
{
  Variables v("|Xr1#|Xr2#|Xs#|m#|h#|j#|V#|d#|f#|fCa#|s#|r#|g#|Ca_SR#|Ca_tot#|Na_i#|K_i#|TRPNCaL#|TRPNCaH#|N_NoXB#|P_NoXB#|N#|XBprer#|XBpostr#|xXBprer#|xXBpostr#|SL#|intf#");
  return v;
}
Variables Solveode::get_Parameters()
{
  Variables v("|time#|R#|T#|F#|Na_o#|K_o#|P_kna#|Ca_o#|g_K1#|g_Kr#|g_Ks#|g_Na#|g_bna#|i_stimAmplitude#|i_stimStart#|i_stimEnd#|i_stimPeriod#|i_stimPulseDuration#|g_CaL#|g_bca#|g_to#|P_NaK#|K_mk#|K_mNa#|K_NaCa#|gamma#|alpha#|Km_Nai#|Km_Ca#|K_sat#|g_pCa#|K_pCa#|g_pK#|a_rel#|b_rel#|c_rel#|Vmax_up#|K_up#|V_leak#|tau_g#|Buf_sr#|K_buf_sr#|V_c#|V_sr#|Cm#|K_buf_c#|Buf_c#|len_thick#|len_thin#|len_hbare#|kon#|Qkon#|TmpC#|koffL#|koffmod#|Qkoff#|koffH#|kn_p#|Qkn_p#|kp_n#|Qkp_n#|perm50#|nperm#|fapp#|xbmodsp#|Qfapp#|gslmod#|gapp#|Qgapp#|hfmdc#|x_0#|hbmdc#|hf#|Qhf#|hb#|Qhb#|sigmap#|sigman#|gxb#|Qgxb#|xPsi#|kxb#|SLrest#|PCon_t#|PExp_t#|PCon_c#|PExp_c#|SL_c#|SLset#|KSE#|SEon#|visc#|massf#|SLmax#|SLmin#|Trop_conc#");
  return v;
}
Variables Solveode::get_FreeVariable()
{
  Variables v("|time#");
  return v;
}

void Solveode::setParametersFromFile(char *filename)
{
  FILE *file;
  if((file = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr,"ERROR - setParametersFromFile - Unable to open file %s\n", filename);
    exit(1);
  }
  double value;
  int k = 0;
  Variables v = get_Parameters();
  int s = v.getQuantity();
  for(;k<s;k++)
  {
    fscanf(file,"%lf", &value);
    setParameters(k, value);
  }
  fclose(file);
}

void Solveode::setVariablesFromFile(char *filename)
{
  FILE *file;
  if((file = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr,"ERROR - setVariablesFromFile - Unable to open file %s\n", filename);
    exit(1);
  }
  double value;
  int k = 0;
  Variables v = get_Variables();
  int s = v.getQuantity();
  for(;k<s;k++)
  {
    fscanf(file,"%lf", &value);
    setVariables(k, value);
  }
  fclose(file);
}

void Solveode::setFreeVariableFromFile(char *filename)
{
  FILE *file;
  if((file = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr,"ERROR - setFreeVariableFromFile - Unable to open file %s\n", filename);
    exit(1);
  }
  double value;
  fscanf(file,"%lf", &value);
  setFreeVariable(value);
  fclose(file);
}

/*
    double* Solveode::solveDiff()
    {

            calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));  //0
            calc_E_K = (((R*T)/F)*log((K_o/K_i_old_))); //1
            calc_E_Ks = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_))))); //2
            calc_xr1_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-2.600000000000000e+01)-V_old_)/7.000000000000000e+00))));  //9
            calc_alpha_xr1 = (4.500000000000000e+02/(1.000000000000000e+00+exp((((-4.500000000000000e+01)-V_old_)/1.000000000000000e+01))));    //10
            calc_beta_xr1 = (6.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.150000000000000e+01))));    //11
            calc_xr2_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+8.800000000000000e+01)/2.400000000000000e+01)))); //14
            calc_alpha_xr2 = (3.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/2.000000000000000e+01))));    //15
            calc_beta_xr2 = (1.120000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01))));    //16
            calc_xs_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/1.400000000000000e+01))));   //20
            calc_alpha_xs = (1.100000000000000e+03/pow((1.000000000000000e+00+exp((((-1.000000000000000e+01)-V_old_)/6.000000000000000e+00))),1.0/2.0));    //21
            calc_beta_xs = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01)))); //22
            calc_m_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp((((-5.686000000000000e+01)-V_old_)/9.029999999999999e+00))),2.000000000000000e+00)); //26
            calc_alpha_m = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/5.000000000000000e+00))));  //27
            calc_beta_m = ((1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_+3.500000000000000e+01)/5.000000000000000e+00))))+(1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_-5.000000000000000e+01)/2.000000000000000e+02)))));    //28
            calc_h_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //31
            calc_alpha_h = ifnumber_0();    //32
            calc_beta_h = ifnumber_1(); //33
            calc_j_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //36
            calc_alpha_j = ifnumber_2();    //37
            calc_beta_j = ifnumber_3(); //38
            calc_i_stim = ifnumber_4(); //42
            calc_d_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/7.500000000000000e+00))));    //45
            calc_alpha_d = ((1.400000000000000e+00/(1.000000000000000e+00+exp((((-3.500000000000000e+01)-V_old_)/1.300000000000000e+01))))+2.500000000000000e-01);  //46
            calc_beta_d = (1.400000000000000e+00/(1.000000000000000e+00+exp(((V_old_+5.000000000000000e+00)/5.000000000000000e+00))));  //47
            calc_gamma_d = (1.000000000000000e+00/(1.000000000000000e+00+exp(((5.000000000000000e+01-V_old_)/2.000000000000000e+01)))); //48
            calc_f_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/7.000000000000000e+00))));   //51
            calc_tau_f = ((1.125000000000000e+03*exp(((-pow((V_old_+2.700000000000000e+01),2.000000000000000e+00))/2.400000000000000e+02)))+8.000000000000000e+01+(1.650000000000000e+02/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/1.000000000000000e+01))))); //52
            calc_tau_fCa = 2.000000000000000e+00;   //58

            if(cell_type !=3){
                calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/5.000000000000000e+00))));   //63
                calc_tau_s = 85.*exp(-(V_old_+45.)*(V_old_+45.)/320.)+5./(1.+exp((V_old_-20.)/5.))+3.;  //64
            }
            else{
                calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00))));   //63
                calc_tau_s = ((1.000000000000000e+03*exp(((-pow((V_old_+6.700000000000000e+01),2.000000000000000e+00))/1.000000000000000e+03)))+8.000000000000000e+00); //64
            }

            calc_r_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((2.000000000000000e+01-V_old_)/6.000000000000000e+00))));   //66
            calc_tau_r = ((9.500000000000000e+00*exp(((-pow((V_old_+4.000000000000000e+01),2.000000000000000e+00))/1.800000000000000e+03)))+8.000000000000000e-01); //67
            calc_i_NaK = (((((P_NaK*K_o)/(K_o+K_mk))*Na_i_old_)/(Na_i_old_+K_mNa))/(1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*V_old_*F)/(R*T))))+(3.530000000000000e-02*exp((((-V_old_)*F)/(R*T))))));    //69
            calc_i_rel = ((((a_rel*pow(Ca_SR_old_,2.000000000000000e+00))/(pow(b_rel,2.000000000000000e+00)+pow(Ca_SR_old_,2.000000000000000e+00)))+c_rel)*d_old_*g_old_);  //73
            calc_Ca_sr_bufsr = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_sr*K_buf_sr)/pow((Ca_SR_old_+K_buf_sr),2.000000000000000e+00))));    //79
            calc_sovr_ze = ifnumber_8();    //88
            calc_sovr_cle = ifnumber_9();   //89
            calc_konT = (kon*pow(Qkon,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //93
            calc_koffLT = (koffL*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //94
            calc_koffHT = (koffH*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //95
            calc_fappT = (fapp*xbmodsp*pow(Qfapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //105
            calc_hfmd = exp((((-xXBprer_old_)/fabs(xXBprer_old_))*hfmdc*pow((xXBprer_old_/x_0),2.000000000000000e+00)));    //108
            calc_hbmd = exp((((xXBpostr_old_-x_0)/fabs((xXBpostr_old_-x_0)))*hbmdc*pow(((xXBpostr_old_-x_0)/x_0),2.000000000000000e+00)));  //109
            calc_gxbmd = ifnumber_11(); //112
            calc_P = (((1.000000000000000e+00-N_old_)-XBprer_old_)-XBpostr_old_);   //121
            calc_SSXBprer = (((hb*fapp)+(gxb*fapp))/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //128
            calc_SSXBpostr = ((fapp*hf)/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //129
            calc_ppforce_t = (((SL_old_-SLrest)/fabs((SL_old_-SLrest)))*PCon_t*(exp((PExp_t*fabs((SL_old_-SLrest))))-1.000000000000000e+00));   //133
            calc_ppforce_c = ifnumber_12(); //134
            calc_preload = ((fabs((SLset-SLrest))/(SLset-SLrest))*PCon_t*(exp((PExp_t*fabs((SLset-SLrest))))-1.000000000000000e+00));   //136
            calc_afterload = ifnumber_13(); //137



            calc_dSL =ifnumber_14();    //139   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ISOSARCOMETRIC  dSL = 0;



            calc_alpha_K1 = (1.000000000000000e-01/(1.000000000000000e+00+exp((6.000000000000000e-02*((V_old_-calc_E_K)-2.000000000000000e+02))))); //4
            calc_beta_K1 = (((3.000000000000000e+00*exp((2.000000000000000e-04*((V_old_-calc_E_K)+1.000000000000000e+02))))+exp((1.000000000000000e-01*((V_old_-calc_E_K)-1.000000000000000e+01))))/(1.000000000000000e+00+exp(((-5.000000000000000e-01)*(V_old_-calc_E_K))))); //5
            calc_i_Kr = (g_Kr*pow((K_o/5.400000000000000e+00),1.0/2.0)*Xr1_old_*Xr2_old_*(V_old_-calc_E_K));    //8
            calc_tau_xr1 = (1.000000000000000e+00*calc_alpha_xr1*calc_beta_xr1);    //12
            calc_tau_xr2 = (1.000000000000000e+00*calc_alpha_xr2*calc_beta_xr2);    //17
            calc_i_Ks = (g_Ks*pow(Xs_old_,2.000000000000000e+00)*(V_old_-calc_E_Ks));   //19
            calc_tau_xs = (1.000000000000000e+00*calc_alpha_xs*calc_beta_xs);   //23
            calc_i_Na = (g_Na*pow(m_old_,3.000000000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));  //25
            calc_tau_m = (1.000000000000000e+00*calc_alpha_m*calc_beta_m);  //29
            calc_tau_h = (1.000000000000000e+00/(calc_alpha_h+calc_beta_h));    //34
            calc_tau_j = (1.000000000000000e+00/(calc_alpha_j+calc_beta_j));    //39
            calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));   //41
            calc_tau_d = ((1.000000000000000e+00*calc_alpha_d*calc_beta_d)+calc_gamma_d);   //49
            calc_i_to = (g_to*r_old_*s_old_*(V_old_-calc_E_K)); //62
            calc_i_p_K = ((g_pK*(V_old_-calc_E_K))/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/5.980000000000000e+00))));    //72
            calc_len_sovr = (calc_sovr_ze-calc_sovr_cle);   //90
            calc_hfT = (hf*calc_hfmd*xbmodsp*pow(Qhf,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //110
            calc_hbT = (hb*calc_hbmd*xbmodsp*pow(Qhb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //111
            calc_gxbT = (gxb*calc_gxbmd*xbmodsp*pow(Qgxb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //113
            calc_Fnordv = (kxb*x_0*calc_SSXBpostr); //130
            calc_ppforce = (calc_ppforce_t+calc_ppforce_c); //135
            calc_FrSBXB = ((XBpostr_old_+XBprer_old_)/(calc_SSXBpostr+calc_SSXBprer));  //141
            calc_dsovr_ze = ifnumber_15();  //143
            calc_dsovr_cle = ifnumber_16(); //144
            calc_xK1_inf = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));    //6
            calc_SOVFThick = ((calc_len_sovr*2.000000000000000e+00)/(len_thick-len_hbare)); //91
            calc_SOVFThin = (calc_len_sovr/len_thin);   //92
            calc_dXBpostr = ((calc_hfT*XBprer_old_)-((calc_hbT*XBpostr_old_)+(calc_gxbT*XBpostr_old_)));    //120
            calc_dlen_sovr = (calc_dsovr_ze-calc_dsovr_cle);    //145
            calc_i_K1 = (g_K1*calc_xK1_inf*pow((K_o/5.400000000000000e+00),1.0/2.0)*(V_old_-calc_E_K)); //7
            calc_Tropreg = (((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*TRPNCaH_old_)); //104
            calc_gapslmd = (1.000000000000000e+00+((1.000000000000000e+00-calc_SOVFThick)*gslmod)); //106
            calc_force = (kxb*calc_SOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)));    //131
            calc_dSOVFThin = (calc_dlen_sovr/len_thin); //146
            calc_dSOVFThick = ((2.000000000000000e+00*calc_dlen_sovr)/(len_thick-len_hbare));   //147
            calc_TropTot = (Trop_conc*(((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))));  //148
            calc_gappT = (gapp*calc_gapslmd*xbmodsp*pow(Qgapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //107
            calc_active = ((1.000000000000000e+00*calc_force)/calc_Fnordv); //132
            calc_cc = (((calc_TropTot/1.000000000000000e+03)-Ca_tot_old_)*K_buf_c); //82
            calc_bc = (((K_buf_c+Buf_c)+(calc_TropTot/1.000000000000000e+03))-Ca_tot_old_); //83
            calc_permtot = pow(fabs((1.000000000000000e+00/(1.000000000000000e+00+pow((perm50/calc_Tropreg),nperm)))),1.0/2.0); //98
            calc_dXBprer = (((calc_fappT*calc_P)+(calc_hbT*XBpostr_old_))-((calc_gappT*XBprer_old_)+(calc_hfT*XBprer_old_)));   //118
            calc_dutyprer = (((calc_hbT*calc_fappT)+(calc_gxbT*calc_fappT))/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //126
            calc_dutypostr = ((calc_fappT*calc_hfT)/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //127
            calc_Ca_i = (((-calc_bc)+pow(((calc_bc*calc_bc)-(4.000000000000000e+00*calc_cc)),5.000000000000000e-01))/2.000000000000000e+00);    //84
            calc_inprmt = ifnumber_10();    //99
            calc_dFrSBXB = ((calc_dXBpostr+calc_dXBprer)/(calc_SSXBpostr+calc_SSXBprer));   //142
            calc_Cai = (1.000000000000000e+03*calc_Ca_i);   //85
            calc_kn_pT = (kn_p*calc_permtot*pow(Qkn_p,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //96
            calc_dxXBprer = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutyprer)*((calc_fappT*(-xXBprer_old_))+(calc_hbT*(xXBpostr_old_-(x_0+xXBprer_old_)))))); //123
            calc_dxXBpostr = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutypostr)*calc_hfT*((xXBprer_old_+x_0)-xXBpostr_old_)));    //125
            calc_E_Ca = (((5.000000000000000e-01*R*T)/F)*log((Ca_o/calc_Ca_i)));    //3
            calc_i_CaL = ((((g_CaL*d_old_*f_old_*fCa_old_*4.000000000000000e+00*V_old_*pow(F,2.000000000000000e+00))/(R*T))*((calc_Ca_i*exp(((2.000000000000000e+00*V_old_*F)/(R*T))))-(3.410000000000000e-01*Ca_o)))/(exp(((2.000000000000000e+00*V_old_*F)/(R*T)))-1.000000000000000e+00));   //44
            calc_alpha_fCa = (1.000000000000000e+00/(1.000000000000000e+00+pow((calc_Ca_i/3.250000000000000e-04),8.000000000000000e+00)));  //54
            calc_beta_fCa = (1.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-5.000000000000000e-04)/1.000000000000000e-04)))); //55
            calc_gama_fCa = (2.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-7.500000000000000e-04)/8.000000000000000e-04)))); //56
            calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F)/(R*T)))*pow(Na_i_old_,3.000000000000000e+00)*Ca_o)-(exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,3.000000000000000e+00)*calc_Ca_i*alpha)))/((pow(Km_Nai,3.000000000000000e+00)+pow(Na_o,3.000000000000000e+00))*(Km_Ca+Ca_o)*(1.000000000000000e+00+(K_sat*exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))))));  //70
            calc_i_p_Ca = ((g_pCa*calc_Ca_i)/(calc_Ca_i+K_pCa));    //71
            calc_i_up = (Vmax_up/(1.000000000000000e+00+(pow(K_up,2.000000000000000e+00)/pow(calc_Ca_i,2.000000000000000e+00))));   //74
            calc_i_leak = (V_leak*(Ca_SR_old_-calc_Ca_i));  //75
            calc_g_inf = ifnumber_6();  //76
            calc_kp_nT = (kp_n*calc_inprmt*pow(Qkp_n,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //97
            calc_dTRPNCaL = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaL_old_))-(calc_koffLT*TRPNCaL_old_)); //101
            calc_dTRPNCaH = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaH_old_))-(calc_koffHT*TRPNCaH_old_)); //103
            calc_dforce = ((kxb*calc_dSOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)))+(kxb*calc_SOVFThick*((calc_dxXBpostr*XBpostr_old_)+(xXBpostr_old_*calc_dXBpostr)+(calc_dxXBprer*XBprer_old_)+(xXBprer_old_*calc_dXBprer)))); //150
            calc_fCa_inf = ((calc_alpha_fCa+calc_beta_fCa+calc_gama_fCa+2.300000000000000e-01)/1.460000000000000e+00);  //57
            calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));   //61
            calc_d_g = ((calc_g_inf-g_old_)/tau_g); //77
            calc_dTropTot = (Trop_conc*(((-calc_dSOVFThin)*TRPNCaL_old_)+((1.000000000000000e+00-calc_SOVFThin)*calc_dTRPNCaL)+(calc_dSOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))+(calc_SOVFThin*(((calc_dFrSBXB*TRPNCaH_old_)+(calc_FrSBXB*calc_dTRPNCaH)+((1.000000000000000e+00-calc_FrSBXB)*calc_dTRPNCaL))-(calc_dFrSBXB*TRPNCaL_old_)))));  //149
            calc_d_fCa = ((calc_fCa_inf-fCa_old_)/calc_tau_fCa);    //59

        depvar__[0]= ((calc_xr1_inf-Xr1_old_)/calc_tau_xr1);    // 13
        depvar__[1]= ((calc_xr2_inf-Xr2_old_)/calc_tau_xr2);    // 18
        depvar__[2]= ((calc_xs_inf-Xs_old_)/calc_tau_xs);   // 24
        depvar__[3]= ((calc_m_inf-m_old_)/calc_tau_m);  // 30
        depvar__[4]= ((calc_h_inf-h_old_)/calc_tau_h);  // 35
        depvar__[5]= ((calc_j_inf-j_old_)/calc_tau_j);  // 40
        depvar__[6]= (((-1.000000000000000e+00)/1.000000000000000e+00)*(calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_CaL+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_NaCa+calc_i_b_Ca+calc_i_p_K+calc_i_p_Ca+calc_i_stim));   // 43
        depvar__[7]= ((calc_d_inf-d_old_)/calc_tau_d);  // 50
        depvar__[8]= ((calc_f_inf-f_old_)/calc_tau_f);  // 53
        depvar__[9]= ifnumber_5();  // 60
        depvar__[10]= ((calc_s_inf-s_old_)/calc_tau_s); // 65
        depvar__[11]= ((calc_r_inf-r_old_)/calc_tau_r); // 68
        depvar__[12]= ifnumber_7(); // 78
        depvar__[13]= (((calc_Ca_sr_bufsr*V_c)/V_sr)*(calc_i_up-(calc_i_rel+calc_i_leak))); // 80
        depvar__[14]= (((calc_i_leak-calc_i_up)+calc_i_rel)-(((1.000000000000000e+00*(((calc_i_CaL+calc_i_b_Ca)+calc_i_p_Ca)-(2.000000000000000e+00*calc_i_NaCa)))/(((2.000000000000000e+00*1.000000000000000e+00)*V_c)*F))*Cm));   // 81
        depvar__[15]= ((((-1.000000000000000e+00)*(calc_i_Na+calc_i_b_Na+(3.000000000000000e+00*calc_i_NaK)+(3.000000000000000e+00*calc_i_NaCa)))/(1.000000000000000e+00*V_c*F))*Cm);   // 86
        depvar__[16]= ((((-1.000000000000000e+00)*((calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_p_K+calc_i_stim)-(2.000000000000000e+00*calc_i_NaK)))/(1.000000000000000e+00*V_c*F))*Cm);    // 87
        depvar__[17]= calc_dTRPNCaL;    // 100
        depvar__[18]= calc_dTRPNCaH;    // 102
        depvar__[19]= ((calc_kp_nT*P_NoXB_old_)-(calc_kn_pT*N_NoXB_old_));  // 114
        depvar__[20]= ((calc_kn_pT*N_NoXB_old_)-(calc_kp_nT*P_NoXB_old_));  // 115
        depvar__[21]= ((calc_kp_nT*calc_P)-(calc_kn_pT*N_old_));    // 116
        depvar__[22]= calc_dXBprer; // 117
        depvar__[23]= calc_dXBpostr;    // 119
        depvar__[24]= calc_dxXBprer;    // 122
        depvar__[25]= calc_dxXBpostr;   // 124
        depvar__[26]= 0.0;//calc_dSL;   // 138
        depvar__[27]= ((calc_preload+calc_afterload)-(calc_ppforce+calc_active));   // 140
        return depvar__;

    }
 */

/*

    int Solveode::solveCVODE(int firstcall__, int steps__)
    {
      printf("SOLVING WITH CVODE\n");
        int iout = 0;
        realtype tout = time_new+ dtime;
        //printf("time_new = %f time = %f\n", (double) time_new,time);


//      printf("tout = %f dtime = %f\n",(double)tout, dtime );

        static int num_iterations_bak = 0;
        if(firstcall__){
            time_new = time;
            tout = time + dtime;

            if(steps__ <= 0)
            steps__ = 1;
            if(time_vec__ != NULL)free( time_vec__);
            time_vec__ = (double *)malloc(sizeof(double)*steps__);
            Xr1_old_ = Xr1_ini_;
            if(Xr1 != NULL)free( Xr1);
            Xr1 = (double *)malloc(sizeof(double)*steps__);
            Xr2_old_ = Xr2_ini_;
            if(Xr2 != NULL)free( Xr2);
            Xr2 = (double *)malloc(sizeof(double)*steps__);
            Xs_old_ = Xs_ini_;
            if(Xs != NULL)free( Xs);
            Xs = (double *)malloc(sizeof(double)*steps__);
            m_old_ = m_ini_;
            if(m != NULL)free( m);
            m = (double *)malloc(sizeof(double)*steps__);
            h_old_ = h_ini_;
            if(h != NULL)free( h);
            h = (double *)malloc(sizeof(double)*steps__);
            j_old_ = j_ini_;
            if(j != NULL)free( j);
            j = (double *)malloc(sizeof(double)*steps__);
            V_old_ = V_ini_;
            if(V != NULL)free( V);
            V = (double *)malloc(sizeof(double)*steps__);
            d_old_ = d_ini_;
            if(d != NULL)free( d);
            d = (double *)malloc(sizeof(double)*steps__);
            f_old_ = f_ini_;
            if(f != NULL)free( f);
            f = (double *)malloc(sizeof(double)*steps__);
            fCa_old_ = fCa_ini_;
            if(fCa != NULL)free( fCa);
            fCa = (double *)malloc(sizeof(double)*steps__);
            s_old_ = s_ini_;
            if(s != NULL)free( s);
            s = (double *)malloc(sizeof(double)*steps__);
            r_old_ = r_ini_;
            if(r != NULL)free( r);
            r = (double *)malloc(sizeof(double)*steps__);
            g_old_ = g_ini_;
            if(g != NULL)free( g);
            g = (double *)malloc(sizeof(double)*steps__);
            Ca_SR_old_ = Ca_SR_ini_;
            if(Ca_SR != NULL)free( Ca_SR);
            Ca_SR = (double *)malloc(sizeof(double)*steps__);
            Ca_tot_old_ = Ca_tot_ini_;
            if(Ca_tot != NULL)free( Ca_tot);
            Ca_tot = (double *)malloc(sizeof(double)*steps__);
            Na_i_old_ = Na_i_ini_;
            if(Na_i != NULL)free( Na_i);
            Na_i = (double *)malloc(sizeof(double)*steps__);
            K_i_old_ = K_i_ini_;
            if(K_i != NULL)free( K_i);
            K_i = (double *)malloc(sizeof(double)*steps__);
            TRPNCaL_old_ = TRPNCaL_ini_;
            if(TRPNCaL != NULL)free( TRPNCaL);
            TRPNCaL = (double *)malloc(sizeof(double)*steps__);
            TRPNCaH_old_ = TRPNCaH_ini_;
            if(TRPNCaH != NULL)free( TRPNCaH);
            TRPNCaH = (double *)malloc(sizeof(double)*steps__);
            N_NoXB_old_ = N_NoXB_ini_;
            if(N_NoXB != NULL)free( N_NoXB);
            N_NoXB = (double *)malloc(sizeof(double)*steps__);
            P_NoXB_old_ = P_NoXB_ini_;
            if(P_NoXB != NULL)free( P_NoXB);
            P_NoXB = (double *)malloc(sizeof(double)*steps__);
            N_old_ = N_ini_;
            if(N != NULL)free( N);
            N = (double *)malloc(sizeof(double)*steps__);
            XBprer_old_ = XBprer_ini_;
            if(XBprer != NULL)free( XBprer);
            XBprer = (double *)malloc(sizeof(double)*steps__);
            XBpostr_old_ = XBpostr_ini_;
            if(XBpostr != NULL)free( XBpostr);
            XBpostr = (double *)malloc(sizeof(double)*steps__);
            xXBprer_old_ = xXBprer_ini_;
            if(xXBprer != NULL)free( xXBprer);
            xXBprer = (double *)malloc(sizeof(double)*steps__);
            xXBpostr_old_ = xXBpostr_ini_;
            if(xXBpostr != NULL)free( xXBpostr);
            xXBpostr = (double *)malloc(sizeof(double)*steps__);
            SL_old_ = SL_ini_;
            if(SL != NULL)free( SL);
            SL = (double *)malloc(sizeof(double)*steps__);
            intf_old_ = intf_ini_;
            if(intf != NULL)free( intf);
            intf = (double *)malloc(sizeof(double)*steps__);

            if(Ca_i != NULL)free( Ca_i);
            Ca_i = (double *)malloc(sizeof(double)*steps__);

            if(active_force != NULL)free( active_force);
            active_force = (double *)malloc(sizeof(double)*steps__);

            num_iterations_bak = steps__;
        }
        flag__ = CVodeSetStopTime(cvode_mem_cvode__,tout);
        while(1){
        double time_new_bak;

            flag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new_bak, CV_NORMAL );
            //printf("time_new = %f time_new_bak = %f tout = %f\n", time_new,time_new_bak,tout);
            time_vec__[iout] = time_new_bak;
            Xr1_old_ = Xr1[iout] = NV_Ith_S(dependent_variable__, 0);
            Xr2_old_ = Xr2[iout] = NV_Ith_S(dependent_variable__, 1);
            Xs_old_ = Xs[iout] = NV_Ith_S(dependent_variable__, 2);
            m_old_ = m[iout] = NV_Ith_S(dependent_variable__, 3);
            h_old_ = h[iout] = NV_Ith_S(dependent_variable__, 4);
            j_old_ = j[iout] = NV_Ith_S(dependent_variable__, 5);
            V_old_ = V[iout] = NV_Ith_S(dependent_variable__, 6);
            d_old_ = d[iout] = NV_Ith_S(dependent_variable__, 7);
            f_old_ = f[iout] = NV_Ith_S(dependent_variable__, 8);
            fCa_old_ = fCa[iout] = NV_Ith_S(dependent_variable__, 9);
            s_old_ = s[iout] = NV_Ith_S(dependent_variable__, 10);
            r_old_ = r[iout] = NV_Ith_S(dependent_variable__, 11);
            g_old_ = g[iout] = NV_Ith_S(dependent_variable__, 12);
            Ca_SR_old_ = Ca_SR[iout] = NV_Ith_S(dependent_variable__, 13);
            Ca_tot_old_ = Ca_tot[iout] = NV_Ith_S(dependent_variable__, 14);
            Na_i_old_ = Na_i[iout] = NV_Ith_S(dependent_variable__, 15);
            K_i_old_ = K_i[iout] = NV_Ith_S(dependent_variable__, 16);
            TRPNCaL_old_ = TRPNCaL[iout] = NV_Ith_S(dependent_variable__, 17);
            TRPNCaH_old_ = TRPNCaH[iout] = NV_Ith_S(dependent_variable__, 18);
            N_NoXB_old_ = N_NoXB[iout] = NV_Ith_S(dependent_variable__, 19);
            P_NoXB_old_ = P_NoXB[iout] = NV_Ith_S(dependent_variable__, 20);
            N_old_ = N[iout] = NV_Ith_S(dependent_variable__, 21);
            XBprer_old_ = XBprer[iout] = NV_Ith_S(dependent_variable__, 22);
            XBpostr_old_ = XBpostr[iout] = NV_Ith_S(dependent_variable__, 23);
            xXBprer_old_ = xXBprer[iout] = NV_Ith_S(dependent_variable__, 24);
            xXBpostr_old_ = xXBpostr[iout] = NV_Ith_S(dependent_variable__, 25);
            SL_old_ = SL[iout] = NV_Ith_S(dependent_variable__, 26);
            intf_old_ = intf[iout] = NV_Ith_S(dependent_variable__, 27);
            time_new = time_new_bak;
            Ca_i[iout]=calc_Ca_i;
            active_force[iout]=calc_active;

            if (check_flag(&flag__, "CVode", 1)) break;
            if (flag__ == CV_SUCCESS){
                iout++;
                tout += dtime; // timestep
            }
            if (iout == num_iterations_bak ) break;
        }
        return flag__;
    }
 */

double Solveode::solve(int firstcall__, int num_iterations__, int num_results__)
{
  static int num_iterations_bak = 0;
  static int num_results_bak = 0;
  static int offset_step = 1;

  i_stimAmplitude = 0.0;
  i_stimStart = 0.0;
  i_stimEnd = 0.0;
  i_stimPeriod = 0.0;
  i_stimPulseDuration = 0.0;

  /*
        i_stimAmplitude = -5.200000000000000e+01;
        i_stimStart = 0.000000000000000e+00;
        i_stimEnd = 5.000000000000000e+04;
        i_stimPeriod = 1.000000000000000e+03;
        i_stimPulseDuration = 1.000000000000000e+00;
   */

  if(firstcall__==1){
    time_new = time;

    if(num_results__ <= 0)
      num_results__ = 1;
    if(num_iterations__ <= 0)
      num_iterations__ = 1;
    offset_step = num_iterations__ / num_results__;

    if(time_vec__ != NULL)free( time_vec__);
    time_vec__ = (double *)malloc(sizeof(double)*num_results__);

    Xr1_old_ = Xr1_ini_;
    if(Xr1 != NULL)free( Xr1);
    Xr1 = (double *)malloc(sizeof(double)*num_results__);

    Xr2_old_ = Xr2_ini_;
    if(Xr2 != NULL)free( Xr2);
    Xr2 = (double *)malloc(sizeof(double)*num_results__);
    Xs_old_ = Xs_ini_;
    if(Xs != NULL)free( Xs);
    Xs = (double *)malloc(sizeof(double)*num_results__);
    m_old_ = m_ini_;
    if(m != NULL)free( m);
    m = (double *)malloc(sizeof(double)*num_results__);
    h_old_ = h_ini_;
    if(h != NULL)free( h);
    h = (double *)malloc(sizeof(double)*num_results__);
    j_old_ = j_ini_;
    if(j != NULL)free( j);
    j = (double *)malloc(sizeof(double)*num_results__);
    V_old_ = V_ini_;
    if(V != NULL)free( V);
    V = (double *)malloc(sizeof(double)*num_results__);
    d_old_ = d_ini_;
    if(d != NULL)free( d);
    d = (double *)malloc(sizeof(double)*num_results__);
    f_old_ = f_ini_;
    if(f != NULL)free( f);
    f = (double *)malloc(sizeof(double)*num_results__);
    fCa_old_ = fCa_ini_;
    if(fCa != NULL)free( fCa);
    fCa = (double *)malloc(sizeof(double)*num_results__);
    s_old_ = s_ini_;
    if(s != NULL)free( s);
    s = (double *)malloc(sizeof(double)*num_results__);
    r_old_ = r_ini_;
    if(r != NULL)free( r);
    r = (double *)malloc(sizeof(double)*num_results__);
    g_old_ = g_ini_;
    if(g != NULL)free( g);
    g = (double *)malloc(sizeof(double)*num_results__);
    Ca_SR_old_ = Ca_SR_ini_;
    if(Ca_SR != NULL)free( Ca_SR);
    Ca_SR = (double *)malloc(sizeof(double)*num_results__);
    Ca_tot_old_ = Ca_tot_ini_;
    if(Ca_tot != NULL)free( Ca_tot);
    Ca_tot = (double *)malloc(sizeof(double)*num_results__);
    Na_i_old_ = Na_i_ini_;
    if(Na_i != NULL)free( Na_i);
    Na_i = (double *)malloc(sizeof(double)*num_results__);
    K_i_old_ = K_i_ini_;
    if(K_i != NULL)free( K_i);
    K_i = (double *)malloc(sizeof(double)*num_results__);
    TRPNCaL_old_ = TRPNCaL_ini_;
    if(TRPNCaL != NULL)free( TRPNCaL);
    TRPNCaL = (double *)malloc(sizeof(double)*num_results__);
    TRPNCaH_old_ = TRPNCaH_ini_;
    if(TRPNCaH != NULL)free( TRPNCaH);
    TRPNCaH = (double *)malloc(sizeof(double)*num_results__);
    N_NoXB_old_ = N_NoXB_ini_;
    if(N_NoXB != NULL)free( N_NoXB);
    N_NoXB = (double *)malloc(sizeof(double)*num_results__);
    P_NoXB_old_ = P_NoXB_ini_;
    if(P_NoXB != NULL)free( P_NoXB);
    P_NoXB = (double *)malloc(sizeof(double)*num_results__);
    N_old_ = N_ini_;
    if(N != NULL)free( N);
    N = (double *)malloc(sizeof(double)*num_results__);
    XBprer_old_ = XBprer_ini_;
    if(XBprer != NULL)free( XBprer);
    XBprer = (double *)malloc(sizeof(double)*num_results__);
    XBpostr_old_ = XBpostr_ini_;
    if(XBpostr != NULL)free( XBpostr);
    XBpostr = (double *)malloc(sizeof(double)*num_results__);
    xXBprer_old_ = xXBprer_ini_;
    if(xXBprer != NULL)free( xXBprer);
    xXBprer = (double *)malloc(sizeof(double)*num_results__);
    xXBpostr_old_ = xXBpostr_ini_;
    if(xXBpostr != NULL)free( xXBpostr);
    xXBpostr = (double *)malloc(sizeof(double)*num_results__);
    SL_old_ = SL_ini_;
    if(SL != NULL)free( SL);
    SL = (double *)malloc(sizeof(double)*num_results__);
    intf_old_ = intf_ini_;
    if(intf != NULL)free( intf);
    intf = (double *)malloc(sizeof(double)*num_results__);


    if(Ca_i != NULL)free( Ca_i);
    Ca_i = (double *)malloc(sizeof(double)*num_results__);

    if(active_force != NULL)free( active_force);
    active_force = (double *)malloc(sizeof(double)*num_results__);

    num_results_bak = num_results__;
    num_iterations_bak = num_iterations__;
  }

  int counter_it__ = 0;
  int aux = num_iterations_bak%num_results_bak;
  for(int it_countx = 1; it_countx<=num_iterations_bak; it_countx++){

    time_new += dtime;

    calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));  //0
    calc_E_K = (((R*T)/F)*log((K_o/K_i_old_))); //1
    calc_E_Ks = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_))))); //2
    calc_xr1_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-2.600000000000000e+01)-V_old_)/7.000000000000000e+00))));  //9
    calc_alpha_xr1 = (4.500000000000000e+02/(1.000000000000000e+00+exp((((-4.500000000000000e+01)-V_old_)/1.000000000000000e+01))));    //10
    calc_beta_xr1 = (6.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.150000000000000e+01))));    //11
    calc_xr2_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+8.800000000000000e+01)/2.400000000000000e+01)))); //14
    calc_alpha_xr2 = (3.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/2.000000000000000e+01))));    //15
    calc_beta_xr2 = (1.120000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01))));    //16
    calc_xs_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/1.400000000000000e+01))));   //20
    calc_alpha_xs = (1.100000000000000e+03/pow((1.000000000000000e+00+exp((((-1.000000000000000e+01)-V_old_)/6.000000000000000e+00))),1.0/2.0));    //21
    calc_beta_xs = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01)))); //22
    calc_m_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp((((-5.686000000000000e+01)-V_old_)/9.029999999999999e+00))),2.000000000000000e+00)); //26
    calc_alpha_m = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/5.000000000000000e+00))));  //27
    calc_beta_m = ((1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_+3.500000000000000e+01)/5.000000000000000e+00))))+(1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_-5.000000000000000e+01)/2.000000000000000e+02)))));    //28
    calc_h_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //31
    calc_alpha_h = ifnumber_0();    //32
    calc_beta_h = ifnumber_1(); //33
    calc_j_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //36
    calc_alpha_j = ifnumber_2();    //37
    calc_beta_j = ifnumber_3(); //38
    calc_i_stim = ifnumber_4(); //42
    calc_d_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/7.500000000000000e+00))));    //45
    calc_alpha_d = ((1.400000000000000e+00/(1.000000000000000e+00+exp((((-3.500000000000000e+01)-V_old_)/1.300000000000000e+01))))+2.500000000000000e-01);  //46
    calc_beta_d = (1.400000000000000e+00/(1.000000000000000e+00+exp(((V_old_+5.000000000000000e+00)/5.000000000000000e+00))));  //47
    calc_gamma_d = (1.000000000000000e+00/(1.000000000000000e+00+exp(((5.000000000000000e+01-V_old_)/2.000000000000000e+01)))); //48
    calc_f_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/7.000000000000000e+00))));   //51
    calc_tau_f = ((1.125000000000000e+03*exp(((-pow((V_old_+2.700000000000000e+01),2.000000000000000e+00))/2.400000000000000e+02)))+8.000000000000000e+01+(1.650000000000000e+02/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/1.000000000000000e+01))))); //52
    calc_tau_fCa = 2.000000000000000e+00;   //58


    if(cell_type !=3){
      calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/5.000000000000000e+00))));   //63
      calc_tau_s = 85.*exp(-(V_old_+45.)*(V_old_+45.)/320.)+5./(1.+exp((V_old_-20.)/5.))+3.;  //64
    }
    else{
      calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00))));   //63
      calc_tau_s = ((1.000000000000000e+03*exp(((-pow((V_old_+6.700000000000000e+01),2.000000000000000e+00))/1.000000000000000e+03)))+8.000000000000000e+00); //64
    }


    //calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00)))); //63
    //calc_tau_s = ((1.000000000000000e+03*exp(((-pow((V_old_+6.700000000000000e+01),2.000000000000000e+00))/1.000000000000000e+03)))+8.000000000000000e+00);   //64
    calc_r_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((2.000000000000000e+01-V_old_)/6.000000000000000e+00))));   //66
    calc_tau_r = ((9.500000000000000e+00*exp(((-pow((V_old_+4.000000000000000e+01),2.000000000000000e+00))/1.800000000000000e+03)))+8.000000000000000e-01); //67
    calc_i_NaK = (((((P_NaK*K_o)/(K_o+K_mk))*Na_i_old_)/(Na_i_old_+K_mNa))/(1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*V_old_*F)/(R*T))))+(3.530000000000000e-02*exp((((-V_old_)*F)/(R*T))))));    //69
    calc_i_rel = ((((a_rel*pow(Ca_SR_old_,2.000000000000000e+00))/(pow(b_rel,2.000000000000000e+00)+pow(Ca_SR_old_,2.000000000000000e+00)))+c_rel)*d_old_*g_old_);  //73
    calc_Ca_sr_bufsr = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_sr*K_buf_sr)/pow((Ca_SR_old_+K_buf_sr),2.000000000000000e+00))));    //79
    calc_sovr_ze = ifnumber_8();    //88
    calc_sovr_cle = ifnumber_9();   //89
    calc_konT = (kon*pow(Qkon,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //93
    calc_koffLT = (koffL*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //94
    calc_koffHT = (koffH*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //95
    calc_fappT = (fapp*xbmodsp*pow(Qfapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //105
    calc_hfmd = exp((((-xXBprer_old_)/fabs(xXBprer_old_))*hfmdc*pow((xXBprer_old_/x_0),2.000000000000000e+00)));    //108
    calc_hbmd = exp((((xXBpostr_old_-x_0)/fabs((xXBpostr_old_-x_0)))*hbmdc*pow(((xXBpostr_old_-x_0)/x_0),2.000000000000000e+00)));  //109
    calc_gxbmd = ifnumber_11(); //112
    calc_P = (((1.000000000000000e+00-N_old_)-XBprer_old_)-XBpostr_old_);   //121
    calc_SSXBprer = (((hb*fapp)+(gxb*fapp))/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //128
    calc_SSXBpostr = ((fapp*hf)/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //129
    calc_ppforce_t = (((SL_old_-SLrest)/fabs((SL_old_-SLrest)))*PCon_t*(exp((PExp_t*fabs((SL_old_-SLrest))))-1.000000000000000e+00));   //133
    calc_ppforce_c = ifnumber_12(); //134
    calc_preload = ((fabs((SLset-SLrest))/(SLset-SLrest))*PCon_t*(exp((PExp_t*fabs((SLset-SLrest))))-1.000000000000000e+00));   //136
    calc_afterload = ifnumber_13(); //137
    calc_dSL = ifnumber_14();   //139
    calc_alpha_K1 = (1.000000000000000e-01/(1.000000000000000e+00+exp((6.000000000000000e-02*((V_old_-calc_E_K)-2.000000000000000e+02))))); //4
    calc_beta_K1 = (((3.000000000000000e+00*exp((2.000000000000000e-04*((V_old_-calc_E_K)+1.000000000000000e+02))))+exp((1.000000000000000e-01*((V_old_-calc_E_K)-1.000000000000000e+01))))/(1.000000000000000e+00+exp(((-5.000000000000000e-01)*(V_old_-calc_E_K))))); //5
    calc_i_Kr = (g_Kr*pow((K_o/5.400000000000000e+00),1.0/2.0)*Xr1_old_*Xr2_old_*(V_old_-calc_E_K));    //8
    calc_tau_xr1 = (1.000000000000000e+00*calc_alpha_xr1*calc_beta_xr1);    //12
    calc_tau_xr2 = (1.000000000000000e+00*calc_alpha_xr2*calc_beta_xr2);    //17
    calc_i_Ks = (g_Ks*pow(Xs_old_,2.000000000000000e+00)*(V_old_-calc_E_Ks));   //19
    calc_tau_xs = (1.000000000000000e+00*calc_alpha_xs*calc_beta_xs);   //23
    calc_i_Na = (g_Na*pow(m_old_,3.000000000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));  //25
    calc_tau_m = (1.000000000000000e+00*calc_alpha_m*calc_beta_m);  //29
    calc_tau_h = (1.000000000000000e+00/(calc_alpha_h+calc_beta_h));    //34
    calc_tau_j = (1.000000000000000e+00/(calc_alpha_j+calc_beta_j));    //39
    calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));   //41
    calc_tau_d = ((1.000000000000000e+00*calc_alpha_d*calc_beta_d)+calc_gamma_d);   //49
    calc_i_to = (g_to*r_old_*s_old_*(V_old_-calc_E_K)); //62
    calc_i_p_K = ((g_pK*(V_old_-calc_E_K))/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/5.980000000000000e+00))));    //72
    calc_len_sovr = (calc_sovr_ze-calc_sovr_cle);   //90
    calc_hfT = (hf*calc_hfmd*xbmodsp*pow(Qhf,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //110
    calc_hbT = (hb*calc_hbmd*xbmodsp*pow(Qhb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //111
    calc_gxbT = (gxb*calc_gxbmd*xbmodsp*pow(Qgxb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //113
    calc_Fnordv = (kxb*x_0*calc_SSXBpostr); //130
    calc_ppforce = (calc_ppforce_t+calc_ppforce_c); //135
    calc_FrSBXB = ((XBpostr_old_+XBprer_old_)/(calc_SSXBpostr+calc_SSXBprer));  //141
    calc_dsovr_ze = ifnumber_15();  //143
    calc_dsovr_cle = ifnumber_16(); //144
    calc_xK1_inf = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));    //6
    calc_SOVFThick = ((calc_len_sovr*2.000000000000000e+00)/(len_thick-len_hbare)); //91
    calc_SOVFThin = (calc_len_sovr/len_thin);   //92
    calc_dXBpostr = ((calc_hfT*XBprer_old_)-((calc_hbT*XBpostr_old_)+(calc_gxbT*XBpostr_old_)));    //120
    calc_dlen_sovr = (calc_dsovr_ze-calc_dsovr_cle);    //145
    calc_i_K1 = (g_K1*calc_xK1_inf*pow((K_o/5.400000000000000e+00),1.0/2.0)*(V_old_-calc_E_K)); //7
    calc_Tropreg = (((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*TRPNCaH_old_)); //104
    calc_gapslmd = (1.000000000000000e+00+((1.000000000000000e+00-calc_SOVFThick)*gslmod)); //106
    calc_force = (kxb*calc_SOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)));    //131
    calc_dSOVFThin = (calc_dlen_sovr/len_thin); //146
    calc_dSOVFThick = ((2.000000000000000e+00*calc_dlen_sovr)/(len_thick-len_hbare));   //147
    calc_TropTot = (Trop_conc*(((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))));  //148
    calc_gappT = (gapp*calc_gapslmd*xbmodsp*pow(Qgapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //107
    calc_active = ((1.000000000000000e+00*calc_force)/calc_Fnordv); //132
    calc_cc = (((calc_TropTot/1.000000000000000e+03)-Ca_tot_old_)*K_buf_c); //82
    calc_bc = (((K_buf_c+Buf_c)+(calc_TropTot/1.000000000000000e+03))-Ca_tot_old_); //83
    calc_permtot = pow(fabs((1.000000000000000e+00/(1.000000000000000e+00+pow((perm50/calc_Tropreg),nperm)))),1.0/2.0); //98
    calc_dXBprer = (((calc_fappT*calc_P)+(calc_hbT*XBpostr_old_))-((calc_gappT*XBprer_old_)+(calc_hfT*XBprer_old_)));   //118
    calc_dutyprer = (((calc_hbT*calc_fappT)+(calc_gxbT*calc_fappT))/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //126
    calc_dutypostr = ((calc_fappT*calc_hfT)/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //127
    calc_Ca_i = (((-calc_bc)+pow(((calc_bc*calc_bc)-(4.000000000000000e+00*calc_cc)),5.000000000000000e-01))/2.000000000000000e+00);    //84

    //printf("calc_Ca_i=%e \t SL_old_ = %e\n",calc_Ca_i, SL_old_);

    calc_inprmt = ifnumber_10();    //99
    calc_dFrSBXB = ((calc_dXBpostr+calc_dXBprer)/(calc_SSXBpostr+calc_SSXBprer));   //142
    calc_Cai = (1.000000000000000e+03*calc_Ca_i);   //85
    calc_kn_pT = (kn_p*calc_permtot*pow(Qkn_p,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //96
    calc_dxXBprer = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutyprer)*((calc_fappT*(-xXBprer_old_))+(calc_hbT*(xXBpostr_old_-(x_0+xXBprer_old_)))))); //123
    calc_dxXBpostr = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutypostr)*calc_hfT*((xXBprer_old_+x_0)-xXBpostr_old_)));    //125
    calc_E_Ca = (((5.000000000000000e-01*R*T)/F)*log((Ca_o/calc_Ca_i)));    //3
    calc_i_CaL = ((((g_CaL*d_old_*f_old_*fCa_old_*4.000000000000000e+00*V_old_*pow(F,2.000000000000000e+00))/(R*T))*((calc_Ca_i*exp(((2.000000000000000e+00*V_old_*F)/(R*T))))-(3.410000000000000e-01*Ca_o)))/(exp(((2.000000000000000e+00*V_old_*F)/(R*T)))-1.000000000000000e+00));   //44
    calc_alpha_fCa = (1.000000000000000e+00/(1.000000000000000e+00+pow((calc_Ca_i/3.250000000000000e-04),8.000000000000000e+00)));  //54
    calc_beta_fCa = (1.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-5.000000000000000e-04)/1.000000000000000e-04)))); //55
    calc_gama_fCa = (2.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-7.500000000000000e-04)/8.000000000000000e-04)))); //56
    calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F)/(R*T)))*pow(Na_i_old_,3.000000000000000e+00)*Ca_o)-(exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,3.000000000000000e+00)*calc_Ca_i*alpha)))/((pow(Km_Nai,3.000000000000000e+00)+pow(Na_o,3.000000000000000e+00))*(Km_Ca+Ca_o)*(1.000000000000000e+00+(K_sat*exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))))));  //70
    calc_i_p_Ca = ((g_pCa*calc_Ca_i)/(calc_Ca_i+K_pCa));    //71
    calc_i_up = (Vmax_up/(1.000000000000000e+00+(pow(K_up,2.000000000000000e+00)/pow(calc_Ca_i,2.000000000000000e+00))));   //74
    calc_i_leak = (V_leak*(Ca_SR_old_-calc_Ca_i));  //75
    calc_g_inf = ifnumber_6();  //76
    calc_kp_nT = (kp_n*calc_inprmt*pow(Qkp_n,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //97
    calc_dTRPNCaL = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaL_old_))-(calc_koffLT*TRPNCaL_old_)); //101
    calc_dTRPNCaH = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaH_old_))-(calc_koffHT*TRPNCaH_old_)); //103
    calc_dforce = ((kxb*calc_dSOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)))+(kxb*calc_SOVFThick*((calc_dxXBpostr*XBpostr_old_)+(xXBpostr_old_*calc_dXBpostr)+(calc_dxXBprer*XBprer_old_)+(xXBprer_old_*calc_dXBprer)))); //150
    calc_fCa_inf = ((calc_alpha_fCa+calc_beta_fCa+calc_gama_fCa+2.300000000000000e-01)/1.460000000000000e+00);  //57
    calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));   //61
    calc_d_g = ((calc_g_inf-g_old_)/tau_g); //77
    calc_dTropTot = (Trop_conc*(((-calc_dSOVFThin)*TRPNCaL_old_)+((1.000000000000000e+00-calc_SOVFThin)*calc_dTRPNCaL)+(calc_dSOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))+(calc_SOVFThin*(((calc_dFrSBXB*TRPNCaH_old_)+(calc_FrSBXB*calc_dTRPNCaH)+((1.000000000000000e+00-calc_FrSBXB)*calc_dTRPNCaL))-(calc_dFrSBXB*TRPNCaL_old_)))));  //149
    calc_d_fCa = ((calc_fCa_inf-fCa_old_)/calc_tau_fCa);    //59

    Xr1_new_=dtime*(((calc_xr1_inf-Xr1_old_)/calc_tau_xr1))+Xr1_old_;   // 13
    Xr2_new_=dtime*(((calc_xr2_inf-Xr2_old_)/calc_tau_xr2))+Xr2_old_;   // 18
    Xs_new_=dtime*(((calc_xs_inf-Xs_old_)/calc_tau_xs))+Xs_old_;    // 24
    m_new_=dtime*(((calc_m_inf-m_old_)/calc_tau_m))+m_old_; // 30
    h_new_=dtime*(((calc_h_inf-h_old_)/calc_tau_h))+h_old_; // 35
    j_new_=dtime*(((calc_j_inf-j_old_)/calc_tau_j))+j_old_; // 40
    V_new_=dtime*((((-1.000000000000000e+00)/1.000000000000000e+00)*(calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_CaL+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_NaCa+calc_i_b_Ca+calc_i_p_K+calc_i_p_Ca+calc_i_stim)))+V_old_;  // 43
    d_new_=dtime*(((calc_d_inf-d_old_)/calc_tau_d))+d_old_; // 50
    f_new_=dtime*(((calc_f_inf-f_old_)/calc_tau_f))+f_old_; // 53
    fCa_new_=dtime*(ifnumber_5())+fCa_old_; // 60
    s_new_=dtime*(((calc_s_inf-s_old_)/calc_tau_s))+s_old_; // 65
    r_new_=dtime*(((calc_r_inf-r_old_)/calc_tau_r))+r_old_; // 68
    g_new_=dtime*(ifnumber_7())+g_old_; // 78
    Ca_SR_new_=dtime*((((calc_Ca_sr_bufsr*V_c)/V_sr)*(calc_i_up-(calc_i_rel+calc_i_leak))))+Ca_SR_old_; // 80
    Ca_tot_new_=dtime*((((calc_i_leak-calc_i_up)+calc_i_rel)-(((1.000000000000000e+00*(((calc_i_CaL+calc_i_b_Ca)+calc_i_p_Ca)-(2.000000000000000e+00*calc_i_NaCa)))/(((2.000000000000000e+00*1.000000000000000e+00)*V_c)*F))*Cm)))+Ca_tot_old_; // 81
    Na_i_new_=dtime*(((((-1.000000000000000e+00)*(calc_i_Na+calc_i_b_Na+(3.000000000000000e+00*calc_i_NaK)+(3.000000000000000e+00*calc_i_NaCa)))/(1.000000000000000e+00*V_c*F))*Cm))+Na_i_old_; // 86
    K_i_new_=dtime*(((((-1.000000000000000e+00)*((calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_p_K+calc_i_stim)-(2.000000000000000e+00*calc_i_NaK)))/(1.000000000000000e+00*V_c*F))*Cm))+K_i_old_;    // 87
    TRPNCaL_new_=dtime*(calc_dTRPNCaL)+TRPNCaL_old_;    // 100
    TRPNCaH_new_=dtime*(calc_dTRPNCaH)+TRPNCaH_old_;    // 102
    N_NoXB_new_=dtime*(((calc_kp_nT*P_NoXB_old_)-(calc_kn_pT*N_NoXB_old_)))+N_NoXB_old_;    // 114
    P_NoXB_new_=dtime*(((calc_kn_pT*N_NoXB_old_)-(calc_kp_nT*P_NoXB_old_)))+P_NoXB_old_;    // 115
    N_new_=dtime*(((calc_kp_nT*calc_P)-(calc_kn_pT*N_old_)))+N_old_;    // 116
    XBprer_new_=dtime*(calc_dXBprer)+XBprer_old_;   // 117
    XBpostr_new_=dtime*(calc_dXBpostr)+XBpostr_old_;    // 119
    xXBprer_new_=dtime*(calc_dxXBprer)+xXBprer_old_;    // 122
    xXBpostr_new_=dtime*(calc_dxXBpostr)+xXBpostr_old_; // 124
    SL_new_=SL_old_;//dtime*(calc_dSL)+SL_old_; // 138
    intf_new_=dtime*(((calc_preload+calc_afterload)-(calc_ppforce+calc_active)))+intf_old_; // 140

    if(it_countx != aux && (it_countx-aux)%offset_step == 0)
    {
      Xr1[counter_it__] = Xr1_new_;
      Xr2[counter_it__] = Xr2_new_;
      Xs[counter_it__] = Xs_new_;
      m[counter_it__] = m_new_;
      h[counter_it__] = h_new_;
      j[counter_it__] = j_new_;
      V[counter_it__] = V_new_;
      d[counter_it__] = d_new_;
      f[counter_it__] = f_new_;
      fCa[counter_it__] = fCa_new_;
      s[counter_it__] = s_new_;
      r[counter_it__] = r_new_;
      g[counter_it__] = g_new_;
      Ca_SR[counter_it__] = Ca_SR_new_;
      Ca_tot[counter_it__] = Ca_tot_new_;
      Na_i[counter_it__] = Na_i_new_;
      K_i[counter_it__] = K_i_new_;
      TRPNCaL[counter_it__] = TRPNCaL_new_;
      TRPNCaH[counter_it__] = TRPNCaH_new_;
      N_NoXB[counter_it__] = N_NoXB_new_;
      P_NoXB[counter_it__] = P_NoXB_new_;
      N[counter_it__] = N_new_;
      XBprer[counter_it__] = XBprer_new_;
      XBpostr[counter_it__] = XBpostr_new_;
      xXBprer[counter_it__] = xXBprer_new_;
      xXBpostr[counter_it__] = xXBpostr_new_;
      SL[counter_it__] = SL_new_;
      intf[counter_it__] = intf_new_;
      time_vec__[counter_it__] = time_new;
      Ca_i[counter_it__]=calc_Ca_i;
      active_force[counter_it__]=calc_active;

      counter_it__++;
    }
    Xr1_old_ = Xr1_new_;
    Xr2_old_ = Xr2_new_;
    Xs_old_ = Xs_new_;
    m_old_ = m_new_;
    h_old_ = h_new_;
    j_old_ = j_new_;
    V_old_ = V_new_;
    d_old_ = d_new_;
    f_old_ = f_new_;
    fCa_old_ = fCa_new_;
    s_old_ = s_new_;
    r_old_ = r_new_;
    g_old_ = g_new_;
    Ca_SR_old_ = Ca_SR_new_;
    Ca_tot_old_ = Ca_tot_new_;
    Na_i_old_ = Na_i_new_;
    K_i_old_ = K_i_new_;
    TRPNCaL_old_ = TRPNCaL_new_;
    TRPNCaH_old_ = TRPNCaH_new_;
    N_NoXB_old_ = N_NoXB_new_;
    P_NoXB_old_ = P_NoXB_new_;
    N_old_ = N_new_;
    XBprer_old_ = XBprer_new_;
    XBpostr_old_ = XBpostr_new_;
    xXBprer_old_ = xXBprer_new_;
    xXBpostr_old_ = xXBpostr_new_;
    SL_old_ = SL_new_;
    intf_old_ = intf_new_;
  }
  return (num_iterations_bak%offset_step)*dtime;
}


double Solveode::solve(int firstcall__, int num_iterations__, int num_results__,
    double stim, double start, double end)
{
  static int num_iterations_bak = 0;
  static int num_results_bak = 0;
  static int offset_step = 1;

  //i_stimAmplitude = -5.200000000000000e+01;
  //i_stimStart = 1.000000000000000e+00;
  //i_stimEnd = 5.000000000000000e+04;
  i_stimAmplitude = stim;
  i_stimStart = start;
  i_stimEnd = end;
  i_stimPeriod = 1.000000000000000e+03;
  i_stimPulseDuration = 1.000000000000000e+00;

  if(firstcall__==1)
  {
    time_new = time;

    if(num_results__ <= 0)
      num_results__ = 1;
    if(num_iterations__ <= 0)
      num_iterations__ = 1;
    offset_step = num_iterations__ / num_results__;

    if(time_vec__ != NULL)free( time_vec__);
    time_vec__ = (double *)malloc(sizeof(double)*num_results__);

    Xr1_old_ = Xr1_ini_;
    if(Xr1 != NULL)free( Xr1);
    Xr1 = (double *)malloc(sizeof(double)*num_results__);

    Xr2_old_ = Xr2_ini_;
    if(Xr2 != NULL)free( Xr2);
    Xr2 = (double *)malloc(sizeof(double)*num_results__);
    Xs_old_ = Xs_ini_;
    if(Xs != NULL)free( Xs);
    Xs = (double *)malloc(sizeof(double)*num_results__);
    m_old_ = m_ini_;
    if(m != NULL)free( m);
    m = (double *)malloc(sizeof(double)*num_results__);
    h_old_ = h_ini_;
    if(h != NULL)free( h);
    h = (double *)malloc(sizeof(double)*num_results__);
    j_old_ = j_ini_;
    if(j != NULL)free( j);
    j = (double *)malloc(sizeof(double)*num_results__);
    V_old_ = V_ini_;
    if(V != NULL)free( V);
    V = (double *)malloc(sizeof(double)*num_results__);
    d_old_ = d_ini_;
    if(d != NULL)free( d);
    d = (double *)malloc(sizeof(double)*num_results__);
    f_old_ = f_ini_;
    if(f != NULL)free( f);
    f = (double *)malloc(sizeof(double)*num_results__);
    fCa_old_ = fCa_ini_;
    if(fCa != NULL)free( fCa);
    fCa = (double *)malloc(sizeof(double)*num_results__);
    s_old_ = s_ini_;
    if(s != NULL)free( s);
    s = (double *)malloc(sizeof(double)*num_results__);
    r_old_ = r_ini_;
    if(r != NULL)free( r);
    r = (double *)malloc(sizeof(double)*num_results__);
    g_old_ = g_ini_;
    if(g != NULL)free( g);
    g = (double *)malloc(sizeof(double)*num_results__);
    Ca_SR_old_ = Ca_SR_ini_;
    if(Ca_SR != NULL)free( Ca_SR);
    Ca_SR = (double *)malloc(sizeof(double)*num_results__);
    Ca_tot_old_ = Ca_tot_ini_;
    if(Ca_tot != NULL)free( Ca_tot);
    Ca_tot = (double *)malloc(sizeof(double)*num_results__);
    Na_i_old_ = Na_i_ini_;
    if(Na_i != NULL)free( Na_i);
    Na_i = (double *)malloc(sizeof(double)*num_results__);
    K_i_old_ = K_i_ini_;
    if(K_i != NULL)free( K_i);
    K_i = (double *)malloc(sizeof(double)*num_results__);
    TRPNCaL_old_ = TRPNCaL_ini_;
    if(TRPNCaL != NULL)free( TRPNCaL);
    TRPNCaL = (double *)malloc(sizeof(double)*num_results__);
    TRPNCaH_old_ = TRPNCaH_ini_;
    if(TRPNCaH != NULL)free( TRPNCaH);
    TRPNCaH = (double *)malloc(sizeof(double)*num_results__);
    N_NoXB_old_ = N_NoXB_ini_;
    if(N_NoXB != NULL)free( N_NoXB);
    N_NoXB = (double *)malloc(sizeof(double)*num_results__);
    P_NoXB_old_ = P_NoXB_ini_;
    if(P_NoXB != NULL)free( P_NoXB);
    P_NoXB = (double *)malloc(sizeof(double)*num_results__);
    N_old_ = N_ini_;
    if(N != NULL)free( N);
    N = (double *)malloc(sizeof(double)*num_results__);
    XBprer_old_ = XBprer_ini_;
    if(XBprer != NULL)free( XBprer);
    XBprer = (double *)malloc(sizeof(double)*num_results__);
    XBpostr_old_ = XBpostr_ini_;
    if(XBpostr != NULL)free( XBpostr);
    XBpostr = (double *)malloc(sizeof(double)*num_results__);
    xXBprer_old_ = xXBprer_ini_;
    if(xXBprer != NULL)free( xXBprer);
    xXBprer = (double *)malloc(sizeof(double)*num_results__);
    xXBpostr_old_ = xXBpostr_ini_;
    if(xXBpostr != NULL)free( xXBpostr);
    xXBpostr = (double *)malloc(sizeof(double)*num_results__);
    SL_old_ = SL_ini_;
    if(SL != NULL)free( SL);
    SL = (double *)malloc(sizeof(double)*num_results__);
    intf_old_ = intf_ini_;
    if(intf != NULL)free( intf);
    intf = (double *)malloc(sizeof(double)*num_results__);


    if(Ca_i != NULL)free( Ca_i);
    Ca_i = (double *)malloc(sizeof(double)*num_results__);

    if(active_force != NULL)free( active_force);
    active_force = (double *)malloc(sizeof(double)*num_results__);

    num_results_bak = num_results__;
    num_iterations_bak = num_iterations__;
  }

  int counter_it__ = 0;
  int aux = num_iterations_bak%num_results_bak;
  for(int it_countx = 1; it_countx<=num_iterations_bak; it_countx++)
  {
    time_new += dtime;

    calc_E_Na      = (((R*T)/F)*log((Na_o/Na_i_old_)));  //0
    calc_E_K       = (((R*T)/F)*log((K_o/K_i_old_))); //1
    calc_E_Ks      = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_))))); //2
    calc_xr1_inf   = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-2.600000000000000e+01)-V_old_)/7.000000000000000e+00))));  //9
    calc_alpha_xr1 = (4.500000000000000e+02/(1.000000000000000e+00+exp((((-4.500000000000000e+01)-V_old_)/1.000000000000000e+01))));    //10
    calc_beta_xr1  = (6.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.150000000000000e+01))));    //11
    calc_xr2_inf   = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+8.800000000000000e+01)/2.400000000000000e+01)))); //14
    calc_alpha_xr2 = (3.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/2.000000000000000e+01))));    //15
    calc_beta_xr2  = (1.120000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01))));    //16
    calc_xs_inf    = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/1.400000000000000e+01))));   //20
    calc_alpha_xs  = (1.100000000000000e+03/pow((1.000000000000000e+00+exp((((-1.000000000000000e+01)-V_old_)/6.000000000000000e+00))),1.0/2.0));    //21
    calc_beta_xs   = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01)))); //22
    calc_m_inf     = (1.000000000000000e+00/pow((1.000000000000000e+00+exp((((-5.686000000000000e+01)-V_old_)/9.029999999999999e+00))),2.000000000000000e+00)); //26
    calc_alpha_m   = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/5.000000000000000e+00))));  //27
    calc_beta_m    = ((1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_+3.500000000000000e+01)/5.000000000000000e+00))))+(1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_-5.000000000000000e+01)/2.000000000000000e+02)))));    //28
    calc_h_inf     = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //31
    calc_alpha_h   = ifnumber_0();    //32
    calc_beta_h    = ifnumber_1(); //33
    calc_j_inf     = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //36
    calc_alpha_j   = ifnumber_2();    //37
    calc_beta_j    = ifnumber_3(); //38
    calc_i_stim    = ifnumber_4(); //42
    calc_d_inf     = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/7.500000000000000e+00))));    //45
    calc_alpha_d   = ((1.400000000000000e+00/(1.000000000000000e+00+exp((((-3.500000000000000e+01)-V_old_)/1.300000000000000e+01))))+2.500000000000000e-01);  //46
    calc_beta_d    = (1.400000000000000e+00/(1.000000000000000e+00+exp(((V_old_+5.000000000000000e+00)/5.000000000000000e+00))));  //47
    calc_gamma_d   = (1.000000000000000e+00/(1.000000000000000e+00+exp(((5.000000000000000e+01-V_old_)/2.000000000000000e+01)))); //48
    calc_f_inf     = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/7.000000000000000e+00))));   //51
    calc_tau_f     = ((1.125000000000000e+03*exp(((-pow((V_old_+2.700000000000000e+01),2.000000000000000e+00))/2.400000000000000e+02)))+8.000000000000000e+01+(1.650000000000000e+02/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/1.000000000000000e+01))))); //52
    calc_tau_fCa   = 2.000000000000000e+00;   //58


    if(cell_type !=3){
      calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/5.000000000000000e+00))));   //63
      calc_tau_s = 85.*exp(-(V_old_+45.)*(V_old_+45.)/320.)+5./(1.+exp((V_old_-20.)/5.))+3.;  //64
    }
    else{
      calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00))));   //63
      calc_tau_s = ((1.000000000000000e+03*exp(((-pow((V_old_+6.700000000000000e+01),2.000000000000000e+00))/1.000000000000000e+03)))+8.000000000000000e+00); //64
    }


    //calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00)))); //63
    //calc_tau_s = ((1.000000000000000e+03*exp(((-pow((V_old_+6.700000000000000e+01),2.000000000000000e+00))/1.000000000000000e+03)))+8.000000000000000e+00);   //64
    calc_r_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((2.000000000000000e+01-V_old_)/6.000000000000000e+00))));   //66
    calc_tau_r = ((9.500000000000000e+00*exp(((-pow((V_old_+4.000000000000000e+01),2.000000000000000e+00))/1.800000000000000e+03)))+8.000000000000000e-01); //67
    calc_i_NaK = (((((P_NaK*K_o)/(K_o+K_mk))*Na_i_old_)/(Na_i_old_+K_mNa))/(1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*V_old_*F)/(R*T))))+(3.530000000000000e-02*exp((((-V_old_)*F)/(R*T))))));    //69
    calc_i_rel = ((((a_rel*pow(Ca_SR_old_,2.000000000000000e+00))/(pow(b_rel,2.000000000000000e+00)+pow(Ca_SR_old_,2.000000000000000e+00)))+c_rel)*d_old_*g_old_);  //73
    calc_Ca_sr_bufsr = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_sr*K_buf_sr)/pow((Ca_SR_old_+K_buf_sr),2.000000000000000e+00))));    //79
    calc_sovr_ze = ifnumber_8();    //88
    calc_sovr_cle = ifnumber_9();   //89
    calc_konT = (kon*pow(Qkon,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //93
    calc_koffLT = (koffL*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //94
    calc_koffHT = (koffH*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //95
    calc_fappT = (fapp*xbmodsp*pow(Qfapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //105
    calc_hfmd = exp((((-xXBprer_old_)/fabs(xXBprer_old_))*hfmdc*pow((xXBprer_old_/x_0),2.000000000000000e+00)));    //108
    calc_hbmd = exp((((xXBpostr_old_-x_0)/fabs((xXBpostr_old_-x_0)))*hbmdc*pow(((xXBpostr_old_-x_0)/x_0),2.000000000000000e+00)));  //109
    calc_gxbmd = ifnumber_11(); //112
    calc_P = (((1.000000000000000e+00-N_old_)-XBprer_old_)-XBpostr_old_);   //121
    calc_SSXBprer = (((hb*fapp)+(gxb*fapp))/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //128
    calc_SSXBpostr = ((fapp*hf)/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //129
    calc_ppforce_t = (((SL_old_-SLrest)/fabs((SL_old_-SLrest)))*PCon_t*(exp((PExp_t*fabs((SL_old_-SLrest))))-1.000000000000000e+00));   //133
    calc_ppforce_c = ifnumber_12(); //134
    calc_preload = ((fabs((SLset-SLrest))/(SLset-SLrest))*PCon_t*(exp((PExp_t*fabs((SLset-SLrest))))-1.000000000000000e+00));   //136
    calc_afterload = ifnumber_13(); //137
    calc_dSL = ifnumber_14();   //139
    calc_alpha_K1 = (1.000000000000000e-01/(1.000000000000000e+00+exp((6.000000000000000e-02*((V_old_-calc_E_K)-2.000000000000000e+02))))); //4
    calc_beta_K1 = (((3.000000000000000e+00*exp((2.000000000000000e-04*((V_old_-calc_E_K)+1.000000000000000e+02))))+exp((1.000000000000000e-01*((V_old_-calc_E_K)-1.000000000000000e+01))))/(1.000000000000000e+00+exp(((-5.000000000000000e-01)*(V_old_-calc_E_K))))); //5
    calc_i_Kr = (g_Kr*pow((K_o/5.400000000000000e+00),1.0/2.0)*Xr1_old_*Xr2_old_*(V_old_-calc_E_K));    //8
    calc_tau_xr1 = (1.000000000000000e+00*calc_alpha_xr1*calc_beta_xr1);    //12
    calc_tau_xr2 = (1.000000000000000e+00*calc_alpha_xr2*calc_beta_xr2);    //17
    calc_i_Ks = (g_Ks*pow(Xs_old_,2.000000000000000e+00)*(V_old_-calc_E_Ks));   //19
    calc_tau_xs = (1.000000000000000e+00*calc_alpha_xs*calc_beta_xs);   //23
    calc_i_Na = (g_Na*pow(m_old_,3.000000000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));  //25
    calc_tau_m = (1.000000000000000e+00*calc_alpha_m*calc_beta_m);  //29
    calc_tau_h = (1.000000000000000e+00/(calc_alpha_h+calc_beta_h));    //34
    calc_tau_j = (1.000000000000000e+00/(calc_alpha_j+calc_beta_j));    //39
    calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));   //41
    calc_tau_d = ((1.000000000000000e+00*calc_alpha_d*calc_beta_d)+calc_gamma_d);   //49
    calc_i_to = (g_to*r_old_*s_old_*(V_old_-calc_E_K)); //62
    calc_i_p_K = ((g_pK*(V_old_-calc_E_K))/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/5.980000000000000e+00))));    //72
    calc_len_sovr = (calc_sovr_ze-calc_sovr_cle);   //90
    calc_hfT = (hf*calc_hfmd*xbmodsp*pow(Qhf,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //110
    calc_hbT = (hb*calc_hbmd*xbmodsp*pow(Qhb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //111
    calc_gxbT = (gxb*calc_gxbmd*xbmodsp*pow(Qgxb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //113
    calc_Fnordv = (kxb*x_0*calc_SSXBpostr); //130
    calc_ppforce = (calc_ppforce_t+calc_ppforce_c); //135
    calc_FrSBXB = ((XBpostr_old_+XBprer_old_)/(calc_SSXBpostr+calc_SSXBprer));  //141
    calc_dsovr_ze = ifnumber_15();  //143
    calc_dsovr_cle = ifnumber_16(); //144
    calc_xK1_inf = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));    //6
    calc_SOVFThick = ((calc_len_sovr*2.000000000000000e+00)/(len_thick-len_hbare)); //91
    calc_SOVFThin = (calc_len_sovr/len_thin);   //92
    calc_dXBpostr = ((calc_hfT*XBprer_old_)-((calc_hbT*XBpostr_old_)+(calc_gxbT*XBpostr_old_)));    //120
    calc_dlen_sovr = (calc_dsovr_ze-calc_dsovr_cle);    //145
    calc_i_K1 = (g_K1*calc_xK1_inf*pow((K_o/5.400000000000000e+00),1.0/2.0)*(V_old_-calc_E_K)); //7
    calc_Tropreg = (((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*TRPNCaH_old_)); //104
    calc_gapslmd = (1.000000000000000e+00+((1.000000000000000e+00-calc_SOVFThick)*gslmod)); //106
    calc_force = (kxb*calc_SOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)));    //131
    calc_dSOVFThin = (calc_dlen_sovr/len_thin); //146
    calc_dSOVFThick = ((2.000000000000000e+00*calc_dlen_sovr)/(len_thick-len_hbare));   //147
    calc_TropTot = (Trop_conc*(((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))));  //148
    calc_gappT = (gapp*calc_gapslmd*xbmodsp*pow(Qgapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //107
    calc_active = ((1.000000000000000e+00*calc_force)/calc_Fnordv); //132
    calc_cc = (((calc_TropTot/1.000000000000000e+03)-Ca_tot_old_)*K_buf_c); //82
    calc_bc = (((K_buf_c+Buf_c)+(calc_TropTot/1.000000000000000e+03))-Ca_tot_old_); //83
    calc_permtot = pow(fabs((1.000000000000000e+00/(1.000000000000000e+00+pow((perm50/calc_Tropreg),nperm)))),1.0/2.0); //98
    calc_dXBprer = (((calc_fappT*calc_P)+(calc_hbT*XBpostr_old_))-((calc_gappT*XBprer_old_)+(calc_hfT*XBprer_old_)));   //118
    calc_dutyprer = (((calc_hbT*calc_fappT)+(calc_gxbT*calc_fappT))/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //126
    calc_dutypostr = ((calc_fappT*calc_hfT)/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //127
    calc_Ca_i = (((-calc_bc)+pow(((calc_bc*calc_bc)-(4.000000000000000e+00*calc_cc)),5.000000000000000e-01))/2.000000000000000e+00);    //84
    calc_inprmt = ifnumber_10();    //99
    calc_dFrSBXB = ((calc_dXBpostr+calc_dXBprer)/(calc_SSXBpostr+calc_SSXBprer));   //142
    calc_Cai = (1.000000000000000e+03*calc_Ca_i);   //85
    calc_kn_pT = (kn_p*calc_permtot*pow(Qkn_p,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //96
    calc_dxXBprer = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutyprer)*((calc_fappT*(-xXBprer_old_))+(calc_hbT*(xXBpostr_old_-(x_0+xXBprer_old_)))))); //123
    calc_dxXBpostr = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutypostr)*calc_hfT*((xXBprer_old_+x_0)-xXBpostr_old_)));    //125
    calc_E_Ca = (((5.000000000000000e-01*R*T)/F)*log((Ca_o/calc_Ca_i)));    //3
    calc_i_CaL = ((((g_CaL*d_old_*f_old_*fCa_old_*4.000000000000000e+00*V_old_*pow(F,2.000000000000000e+00))/(R*T))*((calc_Ca_i*exp(((2.000000000000000e+00*V_old_*F)/(R*T))))-(3.410000000000000e-01*Ca_o)))/(exp(((2.000000000000000e+00*V_old_*F)/(R*T)))-1.000000000000000e+00));   //44
    calc_alpha_fCa = (1.000000000000000e+00/(1.000000000000000e+00+pow((calc_Ca_i/3.250000000000000e-04),8.000000000000000e+00)));  //54
    calc_beta_fCa = (1.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-5.000000000000000e-04)/1.000000000000000e-04)))); //55
    calc_gama_fCa = (2.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-7.500000000000000e-04)/8.000000000000000e-04)))); //56
    calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F)/(R*T)))*pow(Na_i_old_,3.000000000000000e+00)*Ca_o)-(exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,3.000000000000000e+00)*calc_Ca_i*alpha)))/((pow(Km_Nai,3.000000000000000e+00)+pow(Na_o,3.000000000000000e+00))*(Km_Ca+Ca_o)*(1.000000000000000e+00+(K_sat*exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))))));  //70
    calc_i_p_Ca = ((g_pCa*calc_Ca_i)/(calc_Ca_i+K_pCa));    //71
    calc_i_up = (Vmax_up/(1.000000000000000e+00+(pow(K_up,2.000000000000000e+00)/pow(calc_Ca_i,2.000000000000000e+00))));   //74
    calc_i_leak = (V_leak*(Ca_SR_old_-calc_Ca_i));  //75
    calc_g_inf = ifnumber_6();  //76
    calc_kp_nT = (kp_n*calc_inprmt*pow(Qkp_n,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //97
    calc_dTRPNCaL = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaL_old_))-(calc_koffLT*TRPNCaL_old_)); //101
    calc_dTRPNCaH = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaH_old_))-(calc_koffHT*TRPNCaH_old_)); //103
    calc_dforce = ((kxb*calc_dSOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)))+(kxb*calc_SOVFThick*((calc_dxXBpostr*XBpostr_old_)+(xXBpostr_old_*calc_dXBpostr)+(calc_dxXBprer*XBprer_old_)+(xXBprer_old_*calc_dXBprer)))); //150
    calc_fCa_inf = ((calc_alpha_fCa+calc_beta_fCa+calc_gama_fCa+2.300000000000000e-01)/1.460000000000000e+00);  //57
    calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));   //61
    calc_d_g = ((calc_g_inf-g_old_)/tau_g); //77
    calc_dTropTot = (Trop_conc*(((-calc_dSOVFThin)*TRPNCaL_old_)+((1.000000000000000e+00-calc_SOVFThin)*calc_dTRPNCaL)+(calc_dSOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))+(calc_SOVFThin*(((calc_dFrSBXB*TRPNCaH_old_)+(calc_FrSBXB*calc_dTRPNCaH)+((1.000000000000000e+00-calc_FrSBXB)*calc_dTRPNCaL))-(calc_dFrSBXB*TRPNCaL_old_)))));  //149
    calc_d_fCa = ((calc_fCa_inf-fCa_old_)/calc_tau_fCa);    //59

    Xr1_new_=dtime*(((calc_xr1_inf-Xr1_old_)/calc_tau_xr1))+Xr1_old_;   // 13
    Xr2_new_=dtime*(((calc_xr2_inf-Xr2_old_)/calc_tau_xr2))+Xr2_old_;   // 18
    Xs_new_=dtime*(((calc_xs_inf-Xs_old_)/calc_tau_xs))+Xs_old_;    // 24
    m_new_=dtime*(((calc_m_inf-m_old_)/calc_tau_m))+m_old_; // 30
    h_new_=dtime*(((calc_h_inf-h_old_)/calc_tau_h))+h_old_; // 35
    j_new_=dtime*(((calc_j_inf-j_old_)/calc_tau_j))+j_old_; // 40
    V_new_=dtime*((((-1.000000000000000e+00)/1.000000000000000e+00)*(calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_CaL+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_NaCa+calc_i_b_Ca+calc_i_p_K+calc_i_p_Ca+calc_i_stim)))+V_old_;  // 43
    d_new_=dtime*(((calc_d_inf-d_old_)/calc_tau_d))+d_old_; // 50
    f_new_=dtime*(((calc_f_inf-f_old_)/calc_tau_f))+f_old_; // 53
    fCa_new_=dtime*(ifnumber_5())+fCa_old_; // 60
    s_new_=dtime*(((calc_s_inf-s_old_)/calc_tau_s))+s_old_; // 65
    r_new_=dtime*(((calc_r_inf-r_old_)/calc_tau_r))+r_old_; // 68
    g_new_=dtime*(ifnumber_7())+g_old_; // 78
    Ca_SR_new_=dtime*((((calc_Ca_sr_bufsr*V_c)/V_sr)*(calc_i_up-(calc_i_rel+calc_i_leak))))+Ca_SR_old_; // 80
    Ca_tot_new_=dtime*((((calc_i_leak-calc_i_up)+calc_i_rel)-(((1.000000000000000e+00*(((calc_i_CaL+calc_i_b_Ca)+calc_i_p_Ca)-(2.000000000000000e+00*calc_i_NaCa)))/(((2.000000000000000e+00*1.000000000000000e+00)*V_c)*F))*Cm)))+Ca_tot_old_; // 81
    Na_i_new_=dtime*(((((-1.000000000000000e+00)*(calc_i_Na+calc_i_b_Na+(3.000000000000000e+00*calc_i_NaK)+(3.000000000000000e+00*calc_i_NaCa)))/(1.000000000000000e+00*V_c*F))*Cm))+Na_i_old_; // 86
    K_i_new_=dtime*(((((-1.000000000000000e+00)*((calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_p_K+calc_i_stim)-(2.000000000000000e+00*calc_i_NaK)))/(1.000000000000000e+00*V_c*F))*Cm))+K_i_old_;    // 87
    TRPNCaL_new_=dtime*(calc_dTRPNCaL)+TRPNCaL_old_;    // 100
    TRPNCaH_new_=dtime*(calc_dTRPNCaH)+TRPNCaH_old_;    // 102
    N_NoXB_new_=dtime*(((calc_kp_nT*P_NoXB_old_)-(calc_kn_pT*N_NoXB_old_)))+N_NoXB_old_;    // 114
    P_NoXB_new_=dtime*(((calc_kn_pT*N_NoXB_old_)-(calc_kp_nT*P_NoXB_old_)))+P_NoXB_old_;    // 115
    N_new_=dtime*(((calc_kp_nT*calc_P)-(calc_kn_pT*N_old_)))+N_old_;    // 116
    XBprer_new_=dtime*(calc_dXBprer)+XBprer_old_;   // 117
    XBpostr_new_=dtime*(calc_dXBpostr)+XBpostr_old_;    // 119
    xXBprer_new_=dtime*(calc_dxXBprer)+xXBprer_old_;    // 122
    xXBpostr_new_=dtime*(calc_dxXBpostr)+xXBpostr_old_; // 124
    SL_new_=SL_old_;//dtime*(calc_dSL)+SL_old_; // 138
    intf_new_=dtime*(((calc_preload+calc_afterload)-(calc_ppforce+calc_active)))+intf_old_; // 140

    if(it_countx != aux && (it_countx-aux)%offset_step == 0)
    {
      Xr1[counter_it__] = Xr1_new_;
      Xr2[counter_it__] = Xr2_new_;
      Xs[counter_it__] = Xs_new_;
      m[counter_it__] = m_new_;
      h[counter_it__] = h_new_;
      j[counter_it__] = j_new_;
      V[counter_it__] = V_new_;
      d[counter_it__] = d_new_;
      f[counter_it__] = f_new_;
      fCa[counter_it__] = fCa_new_;
      s[counter_it__] = s_new_;
      r[counter_it__] = r_new_;
      g[counter_it__] = g_new_;
      Ca_SR[counter_it__] = Ca_SR_new_;
      Ca_tot[counter_it__] = Ca_tot_new_;
      Na_i[counter_it__] = Na_i_new_;
      K_i[counter_it__] = K_i_new_;
      TRPNCaL[counter_it__] = TRPNCaL_new_;
      TRPNCaH[counter_it__] = TRPNCaH_new_;
      N_NoXB[counter_it__] = N_NoXB_new_;
      P_NoXB[counter_it__] = P_NoXB_new_;
      N[counter_it__] = N_new_;
      XBprer[counter_it__] = XBprer_new_;
      XBpostr[counter_it__] = XBpostr_new_;
      xXBprer[counter_it__]      = xXBprer_new_;
      xXBpostr[counter_it__]     = xXBpostr_new_;
      SL[counter_it__]           = SL_new_;
      intf[counter_it__]         = intf_new_;
      time_vec__[counter_it__]   = time_new;
      Ca_i[counter_it__]         = calc_Ca_i;
      active_force[counter_it__] = calc_active;

      counter_it__++;
    }

    Xr1_old_ = Xr1_new_;
    Xr2_old_ = Xr2_new_;
    Xs_old_ = Xs_new_;
    m_old_ = m_new_;
    h_old_ = h_new_;
    j_old_ = j_new_;
    V_old_ = V_new_;
    d_old_ = d_new_;
    f_old_ = f_new_;
    fCa_old_ = fCa_new_;
    s_old_ = s_new_;
    r_old_ = r_new_;
    g_old_ = g_new_;
    Ca_SR_old_ = Ca_SR_new_;
    Ca_tot_old_ = Ca_tot_new_;
    Na_i_old_ = Na_i_new_;
    K_i_old_ = K_i_new_;
    TRPNCaL_old_ = TRPNCaL_new_;
    TRPNCaH_old_ = TRPNCaH_new_;
    N_NoXB_old_ = N_NoXB_new_;
    P_NoXB_old_ = P_NoXB_new_;
    N_old_ = N_new_;
    XBprer_old_ = XBprer_new_;
    XBpostr_old_ = XBpostr_new_;
    xXBprer_old_ = xXBprer_new_;
    xXBpostr_old_ = xXBpostr_new_;
    SL_old_ = SL_new_;
    intf_old_ = intf_new_;
  }
  return (num_iterations_bak%offset_step)*dtime;
}


double* Solveode::getIndependentVar()
{
  return time_vec__;
}

double* Solveode::getSolution(int indVariable)
{
  switch(indVariable)
  {
  case 0:  return Xr1;    break;
  case 1:  return Xr2;    break;
  case 2:  return Xs;    break;
  case 3:  return m;    break;
  case 4:  return h;    break;
  case 5:  return j;    break;
  case 6:  return V;    break;
  case 7:  return d;    break;
  case 8:  return f;    break;
  case 9:  return fCa;    break;
  case 10: return s;    break;
  case 11: return r;    break;
  case 12: return g;    break;
  case 13: return Ca_SR;    break;
  case 14: return Ca_tot;    break;
  case 15: return Na_i;    break;
  case 16: return K_i;    break;
  case 17: return TRPNCaL;    break;
  case 18: return TRPNCaH;    break;
  case 19: return N_NoXB;    break;
  case 20: return P_NoXB;    break;
  case 21: return N;    break;
  case 22: return XBprer;    break;
  case 23: return XBpostr;    break;
  case 24: return xXBprer;    break;
  case 25: return xXBpostr;    break;
  case 26: return SL;    break;
  case 27: return intf;    break;
  case 28: return Ca_i;    break;
  case 29: return active_force;    break;
  default: return NULL;    break;
  }
}
double Solveode::solveToFile(char *filename, char *fileaccess, int firstcall__, int num_iterations__, int num_results__)
{

  static int num_iterations_bak = 0;
  static int num_results_bak = 0;
  static int offset_step = 1;
  static char *fileaccess_bak = "";
  if(firstcall__){
    time_new = time;

    if(num_results__ <= 0)
      num_results__ = 1;
    if(num_iterations__ <= 0)
      num_iterations__ = 1;
    offset_step = num_iterations__ / num_results__;
    Xr1_old_ = Xr1_ini_;
    if(Xr1 != NULL)free( Xr1);
    Xr1 = (double *)malloc(sizeof(double)*num_results__);
    Xr2_old_ = Xr2_ini_;
    if(Xr2 != NULL)free( Xr2);
    Xr2 = (double *)malloc(sizeof(double)*num_results__);
    Xs_old_ = Xs_ini_;
    if(Xs != NULL)free( Xs);
    Xs = (double *)malloc(sizeof(double)*num_results__);
    m_old_ = m_ini_;
    if(m != NULL)free( m);
    m = (double *)malloc(sizeof(double)*num_results__);
    h_old_ = h_ini_;
    if(h != NULL)free( h);
    h = (double *)malloc(sizeof(double)*num_results__);
    j_old_ = j_ini_;
    if(j != NULL)free( j);
    j = (double *)malloc(sizeof(double)*num_results__);
    V_old_ = V_ini_;
    if(V != NULL)free( V);
    V = (double *)malloc(sizeof(double)*num_results__);
    d_old_ = d_ini_;
    if(d != NULL)free( d);
    d = (double *)malloc(sizeof(double)*num_results__);
    f_old_ = f_ini_;
    if(f != NULL)free( f);
    f = (double *)malloc(sizeof(double)*num_results__);
    fCa_old_ = fCa_ini_;
    if(fCa != NULL)free( fCa);
    fCa = (double *)malloc(sizeof(double)*num_results__);
    s_old_ = s_ini_;
    if(s != NULL)free( s);
    s = (double *)malloc(sizeof(double)*num_results__);
    r_old_ = r_ini_;
    if(r != NULL)free( r);
    r = (double *)malloc(sizeof(double)*num_results__);
    g_old_ = g_ini_;
    if(g != NULL)free( g);
    g = (double *)malloc(sizeof(double)*num_results__);
    Ca_SR_old_ = Ca_SR_ini_;
    if(Ca_SR != NULL)free( Ca_SR);
    Ca_SR = (double *)malloc(sizeof(double)*num_results__);
    Ca_tot_old_ = Ca_tot_ini_;
    if(Ca_tot != NULL)free( Ca_tot);
    Ca_tot = (double *)malloc(sizeof(double)*num_results__);
    Na_i_old_ = Na_i_ini_;
    if(Na_i != NULL)free( Na_i);
    Na_i = (double *)malloc(sizeof(double)*num_results__);
    K_i_old_ = K_i_ini_;
    if(K_i != NULL)free( K_i);
    K_i = (double *)malloc(sizeof(double)*num_results__);
    TRPNCaL_old_ = TRPNCaL_ini_;
    if(TRPNCaL != NULL)free( TRPNCaL);
    TRPNCaL = (double *)malloc(sizeof(double)*num_results__);
    TRPNCaH_old_ = TRPNCaH_ini_;
    if(TRPNCaH != NULL)free( TRPNCaH);
    TRPNCaH = (double *)malloc(sizeof(double)*num_results__);
    N_NoXB_old_ = N_NoXB_ini_;
    if(N_NoXB != NULL)free( N_NoXB);
    N_NoXB = (double *)malloc(sizeof(double)*num_results__);
    P_NoXB_old_ = P_NoXB_ini_;
    if(P_NoXB != NULL)free( P_NoXB);
    P_NoXB = (double *)malloc(sizeof(double)*num_results__);
    N_old_ = N_ini_;
    if(N != NULL)free( N);
    N = (double *)malloc(sizeof(double)*num_results__);
    XBprer_old_ = XBprer_ini_;
    if(XBprer != NULL)free( XBprer);
    XBprer = (double *)malloc(sizeof(double)*num_results__);
    XBpostr_old_ = XBpostr_ini_;
    if(XBpostr != NULL)free( XBpostr);
    XBpostr = (double *)malloc(sizeof(double)*num_results__);
    xXBprer_old_ = xXBprer_ini_;
    if(xXBprer != NULL)free( xXBprer);
    xXBprer = (double *)malloc(sizeof(double)*num_results__);
    xXBpostr_old_ = xXBpostr_ini_;
    if(xXBpostr != NULL)free( xXBpostr);
    xXBpostr = (double *)malloc(sizeof(double)*num_results__);
    SL_old_ = SL_ini_;
    if(SL != NULL)free( SL);
    SL = (double *)malloc(sizeof(double)*num_results__);
    intf_old_ = intf_ini_;
    if(intf != NULL)free( intf);
    intf = (double *)malloc(sizeof(double)*num_results__);
    num_results_bak = num_results__;
    num_iterations_bak = num_iterations__;
    fileaccess_bak = fileaccess;
  }
  FILE *file = fopen(filename, fileaccess_bak);
  if(!file){
    fprintf(stderr,"ERROR - solveToFile - Unable to open file %s\n",filename);
    exit(1);
  }
  int counter_it__ = 0;
  int aux = num_iterations_bak%num_results_bak;
  for(int it_countx = 1; it_countx<=num_iterations_bak; it_countx++){
    time_new += dtime;

    calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));  //0
    calc_E_K = (((R*T)/F)*log((K_o/K_i_old_))); //1
    calc_E_Ks = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_))))); //2
    calc_xr1_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-2.600000000000000e+01)-V_old_)/7.000000000000000e+00))));  //9
    calc_alpha_xr1 = (4.500000000000000e+02/(1.000000000000000e+00+exp((((-4.500000000000000e+01)-V_old_)/1.000000000000000e+01))));    //10
    calc_beta_xr1 = (6.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+3.000000000000000e+01)/1.150000000000000e+01))));    //11
    calc_xr2_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+8.800000000000000e+01)/2.400000000000000e+01)))); //14
    calc_alpha_xr2 = (3.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/2.000000000000000e+01))));    //15
    calc_beta_xr2 = (1.120000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01))));    //16
    calc_xs_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/1.400000000000000e+01))));   //20
    calc_alpha_xs = (1.100000000000000e+03/pow((1.000000000000000e+00+exp((((-1.000000000000000e+01)-V_old_)/6.000000000000000e+00))),1.0/2.0));    //21
    calc_beta_xs = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_-6.000000000000000e+01)/2.000000000000000e+01)))); //22
    calc_m_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp((((-5.686000000000000e+01)-V_old_)/9.029999999999999e+00))),2.000000000000000e+00)); //26
    calc_alpha_m = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-6.000000000000000e+01)-V_old_)/5.000000000000000e+00))));  //27
    calc_beta_m = ((1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_+3.500000000000000e+01)/5.000000000000000e+00))))+(1.000000000000000e-01/(1.000000000000000e+00+exp(((V_old_-5.000000000000000e+01)/2.000000000000000e+02)))));    //28
    calc_h_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //31
    calc_alpha_h = ifnumber_0();    //32
    calc_beta_h = ifnumber_1(); //33
    calc_j_inf = (1.000000000000000e+00/pow((1.000000000000000e+00+exp(((V_old_+7.155000000000000e+01)/7.430000000000000e+00))),2.000000000000000e+00));    //36
    calc_alpha_j = ifnumber_2();    //37
    calc_beta_j = ifnumber_3(); //38
    calc_i_stim = ifnumber_4(); //42
    calc_d_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp((((-5.000000000000000e+00)-V_old_)/7.500000000000000e+00))));    //45
    calc_alpha_d = ((1.400000000000000e+00/(1.000000000000000e+00+exp((((-3.500000000000000e+01)-V_old_)/1.300000000000000e+01))))+2.500000000000000e-01);  //46
    calc_beta_d = (1.400000000000000e+00/(1.000000000000000e+00+exp(((V_old_+5.000000000000000e+00)/5.000000000000000e+00))));  //47
    calc_gamma_d = (1.000000000000000e+00/(1.000000000000000e+00+exp(((5.000000000000000e+01-V_old_)/2.000000000000000e+01)))); //48
    calc_f_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.000000000000000e+01)/7.000000000000000e+00))));   //51
    calc_tau_f = ((1.125000000000000e+03*exp(((-pow((V_old_+2.700000000000000e+01),2.000000000000000e+00))/2.400000000000000e+02)))+8.000000000000000e+01+(1.650000000000000e+02/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/1.000000000000000e+01))))); //52
    calc_tau_fCa = 2.000000000000000e+00;   //58
    calc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00))));   //63
    calc_tau_s = ((1.000000000000000e+03*exp(((-pow((V_old_+6.700000000000000e+01),2.000000000000000e+00))/1.000000000000000e+03)))+8.000000000000000e+00); //64
    calc_r_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((2.000000000000000e+01-V_old_)/6.000000000000000e+00))));   //66
    calc_tau_r = ((9.500000000000000e+00*exp(((-pow((V_old_+4.000000000000000e+01),2.000000000000000e+00))/1.800000000000000e+03)))+8.000000000000000e-01); //67
    calc_i_NaK = (((((P_NaK*K_o)/(K_o+K_mk))*Na_i_old_)/(Na_i_old_+K_mNa))/(1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*V_old_*F)/(R*T))))+(3.530000000000000e-02*exp((((-V_old_)*F)/(R*T))))));    //69
    calc_i_rel = ((((a_rel*pow(Ca_SR_old_,2.000000000000000e+00))/(pow(b_rel,2.000000000000000e+00)+pow(Ca_SR_old_,2.000000000000000e+00)))+c_rel)*d_old_*g_old_);  //73
    calc_Ca_sr_bufsr = (1.000000000000000e+00/(1.000000000000000e+00+((Buf_sr*K_buf_sr)/pow((Ca_SR_old_+K_buf_sr),2.000000000000000e+00))));    //79
    calc_sovr_ze = ifnumber_8();    //88
    calc_sovr_cle = ifnumber_9();   //89
    calc_konT = (kon*pow(Qkon,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //93
    calc_koffLT = (koffL*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //94
    calc_koffHT = (koffH*koffmod*pow(Qkoff,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));  //95
    calc_fappT = (fapp*xbmodsp*pow(Qfapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //105
    calc_hfmd = exp((((-xXBprer_old_)/fabs(xXBprer_old_))*hfmdc*pow((xXBprer_old_/x_0),2.000000000000000e+00)));    //108
    calc_hbmd = exp((((xXBpostr_old_-x_0)/fabs((xXBpostr_old_-x_0)))*hbmdc*pow(((xXBpostr_old_-x_0)/x_0),2.000000000000000e+00)));  //109
    calc_gxbmd = ifnumber_11(); //112
    calc_P = (((1.000000000000000e+00-N_old_)-XBprer_old_)-XBpostr_old_);   //121
    calc_SSXBprer = (((hb*fapp)+(gxb*fapp))/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //128
    calc_SSXBpostr = ((fapp*hf)/((fapp*hf)+(gxb*hf)+(gxb*gapp)+(hb*fapp)+(hb*gapp)+(gxb*fapp)));    //129
    calc_ppforce_t = (((SL_old_-SLrest)/fabs((SL_old_-SLrest)))*PCon_t*(exp((PExp_t*fabs((SL_old_-SLrest))))-1.000000000000000e+00));   //133
    calc_ppforce_c = ifnumber_12(); //134
    calc_preload = ((fabs((SLset-SLrest))/(SLset-SLrest))*PCon_t*(exp((PExp_t*fabs((SLset-SLrest))))-1.000000000000000e+00));   //136
    calc_afterload = ifnumber_13(); //137
    calc_dSL = ifnumber_14();   //139
    calc_alpha_K1 = (1.000000000000000e-01/(1.000000000000000e+00+exp((6.000000000000000e-02*((V_old_-calc_E_K)-2.000000000000000e+02))))); //4
    calc_beta_K1 = (((3.000000000000000e+00*exp((2.000000000000000e-04*((V_old_-calc_E_K)+1.000000000000000e+02))))+exp((1.000000000000000e-01*((V_old_-calc_E_K)-1.000000000000000e+01))))/(1.000000000000000e+00+exp(((-5.000000000000000e-01)*(V_old_-calc_E_K))))); //5
    calc_i_Kr = (g_Kr*pow((K_o/5.400000000000000e+00),1.0/2.0)*Xr1_old_*Xr2_old_*(V_old_-calc_E_K));    //8
    calc_tau_xr1 = (1.000000000000000e+00*calc_alpha_xr1*calc_beta_xr1);    //12
    calc_tau_xr2 = (1.000000000000000e+00*calc_alpha_xr2*calc_beta_xr2);    //17
    calc_i_Ks = (g_Ks*pow(Xs_old_,2.000000000000000e+00)*(V_old_-calc_E_Ks));   //19
    calc_tau_xs = (1.000000000000000e+00*calc_alpha_xs*calc_beta_xs);   //23
    calc_i_Na = (g_Na*pow(m_old_,3.000000000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));  //25
    calc_tau_m = (1.000000000000000e+00*calc_alpha_m*calc_beta_m);  //29
    calc_tau_h = (1.000000000000000e+00/(calc_alpha_h+calc_beta_h));    //34
    calc_tau_j = (1.000000000000000e+00/(calc_alpha_j+calc_beta_j));    //39
    calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));   //41
    calc_tau_d = ((1.000000000000000e+00*calc_alpha_d*calc_beta_d)+calc_gamma_d);   //49
    calc_i_to = (g_to*r_old_*s_old_*(V_old_-calc_E_K)); //62
    calc_i_p_K = ((g_pK*(V_old_-calc_E_K))/(1.000000000000000e+00+exp(((2.500000000000000e+01-V_old_)/5.980000000000000e+00))));    //72
    calc_len_sovr = (calc_sovr_ze-calc_sovr_cle);   //90
    calc_hfT = (hf*calc_hfmd*xbmodsp*pow(Qhf,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //110
    calc_hbT = (hb*calc_hbmd*xbmodsp*pow(Qhb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //111
    calc_gxbT = (gxb*calc_gxbmd*xbmodsp*pow(Qgxb,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //113
    calc_Fnordv = (kxb*x_0*calc_SSXBpostr); //130
    calc_ppforce = (calc_ppforce_t+calc_ppforce_c); //135
    calc_FrSBXB = ((XBpostr_old_+XBprer_old_)/(calc_SSXBpostr+calc_SSXBprer));  //141
    calc_dsovr_ze = ifnumber_15();  //143
    calc_dsovr_cle = ifnumber_16(); //144
    calc_xK1_inf = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));    //6
    calc_SOVFThick = ((calc_len_sovr*2.000000000000000e+00)/(len_thick-len_hbare)); //91
    calc_SOVFThin = (calc_len_sovr/len_thin);   //92
    calc_dXBpostr = ((calc_hfT*XBprer_old_)-((calc_hbT*XBpostr_old_)+(calc_gxbT*XBpostr_old_)));    //120
    calc_dlen_sovr = (calc_dsovr_ze-calc_dsovr_cle);    //145
    calc_i_K1 = (g_K1*calc_xK1_inf*pow((K_o/5.400000000000000e+00),1.0/2.0)*(V_old_-calc_E_K)); //7
    calc_Tropreg = (((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*TRPNCaH_old_)); //104
    calc_gapslmd = (1.000000000000000e+00+((1.000000000000000e+00-calc_SOVFThick)*gslmod)); //106
    calc_force = (kxb*calc_SOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)));    //131
    calc_dSOVFThin = (calc_dlen_sovr/len_thin); //146
    calc_dSOVFThick = ((2.000000000000000e+00*calc_dlen_sovr)/(len_thick-len_hbare));   //147
    calc_TropTot = (Trop_conc*(((1.000000000000000e+00-calc_SOVFThin)*TRPNCaL_old_)+(calc_SOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))));  //148
    calc_gappT = (gapp*calc_gapslmd*xbmodsp*pow(Qgapp,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //107
    calc_active = ((1.000000000000000e+00*calc_force)/calc_Fnordv); //132
    calc_cc = (((calc_TropTot/1.000000000000000e+03)-Ca_tot_old_)*K_buf_c); //82
    calc_bc = (((K_buf_c+Buf_c)+(calc_TropTot/1.000000000000000e+03))-Ca_tot_old_); //83
    calc_permtot = pow(fabs((1.000000000000000e+00/(1.000000000000000e+00+pow((perm50/calc_Tropreg),nperm)))),1.0/2.0); //98
    calc_dXBprer = (((calc_fappT*calc_P)+(calc_hbT*XBpostr_old_))-((calc_gappT*XBprer_old_)+(calc_hfT*XBprer_old_)));   //118
    calc_dutyprer = (((calc_hbT*calc_fappT)+(calc_gxbT*calc_fappT))/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //126
    calc_dutypostr = ((calc_fappT*calc_hfT)/((calc_fappT*calc_hfT)+(calc_gxbT*calc_hfT)+(calc_gxbT*calc_gappT)+(calc_hbT*calc_fappT)+(calc_hbT*calc_gappT)+(calc_gxbT*calc_fappT)));    //127
    calc_Ca_i = (((-calc_bc)+pow(((calc_bc*calc_bc)-(4.000000000000000e+00*calc_cc)),5.000000000000000e-01))/2.000000000000000e+00);    //84
    calc_inprmt = ifnumber_10();    //99
    calc_dFrSBXB = ((calc_dXBpostr+calc_dXBprer)/(calc_SSXBpostr+calc_SSXBprer));   //142
    calc_Cai = (1.000000000000000e+03*calc_Ca_i);   //85
    calc_kn_pT = (kn_p*calc_permtot*pow(Qkn_p,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));   //96
    calc_dxXBprer = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutyprer)*((calc_fappT*(-xXBprer_old_))+(calc_hbT*(xXBpostr_old_-(x_0+xXBprer_old_)))))); //123
    calc_dxXBpostr = ((calc_dSL/2.000000000000000e+00)+((xPsi/calc_dutypostr)*calc_hfT*((xXBprer_old_+x_0)-xXBpostr_old_)));    //125
    calc_E_Ca = (((5.000000000000000e-01*R*T)/F)*log((Ca_o/calc_Ca_i)));    //3
    calc_i_CaL = ((((g_CaL*d_old_*f_old_*fCa_old_*4.000000000000000e+00*V_old_*pow(F,2.000000000000000e+00))/(R*T))*((calc_Ca_i*exp(((2.000000000000000e+00*V_old_*F)/(R*T))))-(3.410000000000000e-01*Ca_o)))/(exp(((2.000000000000000e+00*V_old_*F)/(R*T)))-1.000000000000000e+00));   //44
    calc_alpha_fCa = (1.000000000000000e+00/(1.000000000000000e+00+pow((calc_Ca_i/3.250000000000000e-04),8.000000000000000e+00)));  //54
    calc_beta_fCa = (1.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-5.000000000000000e-04)/1.000000000000000e-04)))); //55
    calc_gama_fCa = (2.000000000000000e-01/(1.000000000000000e+00+exp(((calc_Ca_i-7.500000000000000e-04)/8.000000000000000e-04)))); //56
    calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F)/(R*T)))*pow(Na_i_old_,3.000000000000000e+00)*Ca_o)-(exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,3.000000000000000e+00)*calc_Ca_i*alpha)))/((pow(Km_Nai,3.000000000000000e+00)+pow(Na_o,3.000000000000000e+00))*(Km_Ca+Ca_o)*(1.000000000000000e+00+(K_sat*exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))))));  //70
    calc_i_p_Ca = ((g_pCa*calc_Ca_i)/(calc_Ca_i+K_pCa));    //71
    calc_i_up = (Vmax_up/(1.000000000000000e+00+(pow(K_up,2.000000000000000e+00)/pow(calc_Ca_i,2.000000000000000e+00))));   //74
    calc_i_leak = (V_leak*(Ca_SR_old_-calc_Ca_i));  //75
    calc_g_inf = ifnumber_6();  //76
    calc_kp_nT = (kp_n*calc_inprmt*pow(Qkp_n,((TmpC-3.700000000000000e+01)/1.000000000000000e+01)));    //97
    calc_dTRPNCaL = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaL_old_))-(calc_koffLT*TRPNCaL_old_)); //101
    calc_dTRPNCaH = ((calc_konT*calc_Cai*(1.000000000000000e+00-TRPNCaH_old_))-(calc_koffHT*TRPNCaH_old_)); //103
    calc_dforce = ((kxb*calc_dSOVFThick*((xXBpostr_old_*XBpostr_old_)+(xXBprer_old_*XBprer_old_)))+(kxb*calc_SOVFThick*((calc_dxXBpostr*XBpostr_old_)+(xXBpostr_old_*calc_dXBpostr)+(calc_dxXBprer*XBprer_old_)+(xXBprer_old_*calc_dXBprer)))); //150
    calc_fCa_inf = ((calc_alpha_fCa+calc_beta_fCa+calc_gama_fCa+2.300000000000000e-01)/1.460000000000000e+00);  //57
    calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));   //61
    calc_d_g = ((calc_g_inf-g_old_)/tau_g); //77
    calc_dTropTot = (Trop_conc*(((-calc_dSOVFThin)*TRPNCaL_old_)+((1.000000000000000e+00-calc_SOVFThin)*calc_dTRPNCaL)+(calc_dSOVFThin*((calc_FrSBXB*TRPNCaH_old_)+((1.000000000000000e+00-calc_FrSBXB)*TRPNCaL_old_)))+(calc_SOVFThin*(((calc_dFrSBXB*TRPNCaH_old_)+(calc_FrSBXB*calc_dTRPNCaH)+((1.000000000000000e+00-calc_FrSBXB)*calc_dTRPNCaL))-(calc_dFrSBXB*TRPNCaL_old_)))));  //149
    calc_d_fCa = ((calc_fCa_inf-fCa_old_)/calc_tau_fCa);    //59

    Xr1_new_=dtime*(((calc_xr1_inf-Xr1_old_)/calc_tau_xr1))+Xr1_old_;   // 13
    Xr2_new_=dtime*(((calc_xr2_inf-Xr2_old_)/calc_tau_xr2))+Xr2_old_;   // 18
    Xs_new_=dtime*(((calc_xs_inf-Xs_old_)/calc_tau_xs))+Xs_old_;    // 24
    m_new_=dtime*(((calc_m_inf-m_old_)/calc_tau_m))+m_old_; // 30
    h_new_=dtime*(((calc_h_inf-h_old_)/calc_tau_h))+h_old_; // 35
    j_new_=dtime*(((calc_j_inf-j_old_)/calc_tau_j))+j_old_; // 40
    V_new_=dtime*((((-1.000000000000000e+00)/1.000000000000000e+00)*(calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_CaL+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_NaCa+calc_i_b_Ca+calc_i_p_K+calc_i_p_Ca+calc_i_stim)))+V_old_;  // 43
    d_new_=dtime*(((calc_d_inf-d_old_)/calc_tau_d))+d_old_; // 50
    f_new_=dtime*(((calc_f_inf-f_old_)/calc_tau_f))+f_old_; // 53
    fCa_new_=dtime*(ifnumber_5())+fCa_old_; // 60
    s_new_=dtime*(((calc_s_inf-s_old_)/calc_tau_s))+s_old_; // 65
    r_new_=dtime*(((calc_r_inf-r_old_)/calc_tau_r))+r_old_; // 68
    g_new_=dtime*(ifnumber_7())+g_old_; // 78
    Ca_SR_new_=dtime*((((calc_Ca_sr_bufsr*V_c)/V_sr)*(calc_i_up-(calc_i_rel+calc_i_leak))))+Ca_SR_old_; // 80
    Ca_tot_new_=dtime*((((calc_i_leak-calc_i_up)+calc_i_rel)-(((1.000000000000000e+00*(((calc_i_CaL+calc_i_b_Ca)+calc_i_p_Ca)-(2.000000000000000e+00*calc_i_NaCa)))/(((2.000000000000000e+00*1.000000000000000e+00)*V_c)*F))*Cm)))+Ca_tot_old_; // 81
    Na_i_new_=dtime*(((((-1.000000000000000e+00)*(calc_i_Na+calc_i_b_Na+(3.000000000000000e+00*calc_i_NaK)+(3.000000000000000e+00*calc_i_NaCa)))/(1.000000000000000e+00*V_c*F))*Cm))+Na_i_old_; // 86
    K_i_new_=dtime*(((((-1.000000000000000e+00)*((calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_p_K+calc_i_stim)-(2.000000000000000e+00*calc_i_NaK)))/(1.000000000000000e+00*V_c*F))*Cm))+K_i_old_;    // 87
    TRPNCaL_new_=dtime*(calc_dTRPNCaL)+TRPNCaL_old_;    // 100
    TRPNCaH_new_=dtime*(calc_dTRPNCaH)+TRPNCaH_old_;    // 102
    N_NoXB_new_=dtime*(((calc_kp_nT*P_NoXB_old_)-(calc_kn_pT*N_NoXB_old_)))+N_NoXB_old_;    // 114
    P_NoXB_new_=dtime*(((calc_kn_pT*N_NoXB_old_)-(calc_kp_nT*P_NoXB_old_)))+P_NoXB_old_;    // 115
    N_new_=dtime*(((calc_kp_nT*calc_P)-(calc_kn_pT*N_old_)))+N_old_;    // 116
    XBprer_new_=dtime*(calc_dXBprer)+XBprer_old_;   // 117
    XBpostr_new_=dtime*(calc_dXBpostr)+XBpostr_old_;    // 119
    xXBprer_new_=dtime*(calc_dxXBprer)+xXBprer_old_;    // 122
    xXBpostr_new_=dtime*(calc_dxXBpostr)+xXBpostr_old_; // 124
    SL_new_=dtime*(calc_dSL)+SL_old_;   // 138
    intf_new_=dtime*(((calc_preload+calc_afterload)-(calc_ppforce+calc_active)))+intf_old_; // 140

    if(it_countx != aux && (it_countx-aux)%offset_step == 0){
      fprintf(file,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",time_new,Xr1_new_,Xr2_new_,Xs_new_,m_new_,h_new_,j_new_,V_new_,d_new_,f_new_,fCa_new_,s_new_,r_new_,g_new_,Ca_SR_new_,Ca_tot_new_,Na_i_new_,K_i_new_,TRPNCaL_new_,TRPNCaH_new_,N_NoXB_new_,P_NoXB_new_,N_new_,XBprer_new_,XBpostr_new_,xXBprer_new_,xXBpostr_new_,SL_new_,intf_new_,calc_Ca_i,calc_active);
      counter_it__++;
    }
    Xr1_old_ = Xr1_new_;
    Xr2_old_ = Xr2_new_;
    Xs_old_ = Xs_new_;
    m_old_ = m_new_;
    h_old_ = h_new_;
    j_old_ = j_new_;
    V_old_ = V_new_;
    d_old_ = d_new_;
    f_old_ = f_new_;
    fCa_old_ = fCa_new_;
    s_old_ = s_new_;
    r_old_ = r_new_;
    g_old_ = g_new_;
    Ca_SR_old_ = Ca_SR_new_;
    Ca_tot_old_ = Ca_tot_new_;
    Na_i_old_ = Na_i_new_;
    K_i_old_ = K_i_new_;
    TRPNCaL_old_ = TRPNCaL_new_;
    TRPNCaH_old_ = TRPNCaH_new_;
    N_NoXB_old_ = N_NoXB_new_;
    P_NoXB_old_ = P_NoXB_new_;
    N_old_ = N_new_;
    XBprer_old_ = XBprer_new_;
    XBpostr_old_ = XBpostr_new_;
    xXBprer_old_ = xXBprer_new_;
    xXBpostr_old_ = xXBpostr_new_;
    SL_old_ = SL_new_;
    intf_old_ = intf_new_;
  }
  fclose(file);
  return (num_iterations_bak%offset_step)*dtime;
}
/*
    void Solveode::solveCVODEToFile(char *filename, char *fileaccess, int firstcall__, int steps__)
    {

        static int num_iterations_bak = 0;
        static char *fileaccess_bak = "";
        if(firstcall__){
            if(steps__ <= 0)
                steps__ = 1;
            num_iterations_bak = steps__;
            fileaccess_bak = fileaccess;
        }
        FILE *file = fopen(filename, fileaccess_bak);
        if(!file){
            fprintf(stderr,"ERROR - solveCVODEToFile - Unable to open file %s\n",filename);
            exit(1);
        }
        solveCVODE(firstcall__, num_iterations_bak);
        for(int i=0;i<num_iterations_bak; i++){
            fprintf(file,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",time_vec__[i],Xr1[i],Xr2[i],Xs[i],m[i],h[i],j[i],V[i],d[i],f[i],fCa[i],s[i],r[i],g[i],Ca_SR[i],Ca_tot[i],Na_i[i],K_i[i],TRPNCaL[i],TRPNCaH[i],N_NoXB[i],P_NoXB[i],N[i],XBprer[i],XBpostr[i],xXBprer[i],xXBpostr[i],SL[i],intf[i]);
        }
        fclose(file);
    }

    void Solveode::reInitCVODE()
    {

        flag__ = CVodeReInit(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SV, reltol__, abstol__);
        if (check_flag(&flag__, "CVodeReInit", 1))
            exit(1);
    }

    void Solveode::setCVODEMaxStep(double maxstep)
    {

            CVodeSetMaxNumSteps(cvode_mem_cvode__, 1000000);

            CVodeSetInitStep(cvode_mem_cvode__, maxstep/10000.0);

        flag__ = CVodeSetMaxStep(cvode_mem_cvode__, maxstep);
        if (check_flag(&flag__, "CVodeSetMaxStep", 1))
            exit(1);

    }
 */

/*
    void Solveode::setReltol(double value_new)
    {
        reltol__ = value_new;
    }
 */
/*
    void Solveode::setAbstol(int index, double value_new)
    {
        switch(index){
        case 0:     NV_Ith_S(abstol__, 0) = value_new;    break;
        case 1:     NV_Ith_S(abstol__, 1) = value_new;    break;
        case 2:     NV_Ith_S(abstol__, 2) = value_new;    break;
        case 3:     NV_Ith_S(abstol__, 3) = value_new;    break;
        case 4:     NV_Ith_S(abstol__, 4) = value_new;    break;
        case 5:     NV_Ith_S(abstol__, 5) = value_new;    break;
        case 6:     NV_Ith_S(abstol__, 6) = value_new;    break;
        case 7:     NV_Ith_S(abstol__, 7) = value_new;    break;
        case 8:     NV_Ith_S(abstol__, 8) = value_new;    break;
        case 9:     NV_Ith_S(abstol__, 9) = value_new;    break;
        case 10:        NV_Ith_S(abstol__, 10) = value_new;    break;
        case 11:        NV_Ith_S(abstol__, 11) = value_new;    break;
        case 12:        NV_Ith_S(abstol__, 12) = value_new;    break;
        case 13:        NV_Ith_S(abstol__, 13) = value_new;    break;
        case 14:        NV_Ith_S(abstol__, 14) = value_new;    break;
        case 15:        NV_Ith_S(abstol__, 15) = value_new;    break;
        case 16:        NV_Ith_S(abstol__, 16) = value_new;    break;
        case 17:        NV_Ith_S(abstol__, 17) = value_new;    break;
        case 18:        NV_Ith_S(abstol__, 18) = value_new;    break;
        case 19:        NV_Ith_S(abstol__, 19) = value_new;    break;
        case 20:        NV_Ith_S(abstol__, 20) = value_new;    break;
        case 21:        NV_Ith_S(abstol__, 21) = value_new;    break;
        case 22:        NV_Ith_S(abstol__, 22) = value_new;    break;
        case 23:        NV_Ith_S(abstol__, 23) = value_new;    break;
        case 24:        NV_Ith_S(abstol__, 24) = value_new;    break;
        case 25:        NV_Ith_S(abstol__, 25) = value_new;    break;
        case 26:        NV_Ith_S(abstol__, 26) = value_new;    break;
        case 27:        NV_Ith_S(abstol__, 27) = value_new;    break;
        default: fprintf(stderr,"ERROR - setAbstol - index = %d out of bounds\n",index);    break;
        }
    }
 */
double Solveode::ifnumber_0(){
  if((V_old_<(-4.000000000000000e+01))){
    return ((5.700000000000000e-02*exp(((-(V_old_+8.000000000000000e+01))/6.800000000000000e+00))));
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_1(){
  if((V_old_<(-4.000000000000000e+01))){
    return (((2.700000000000000e+00*exp((7.900000000000000e-02*V_old_)))+(3.100000000000000e+05*exp((3.485000000000000e-01*V_old_)))));
  }else{
    return ((7.700000000000000e-01/(1.300000000000000e-01*(1.000000000000000e+00+exp(((V_old_+1.066000000000000e+01)/(-1.110000000000000e+01)))))));
  }
}
double Solveode::ifnumber_2(){
  if((V_old_<(-4.000000000000000e+01))){
    return (((((((-2.542800000000000e+04)*exp((2.444000000000000e-01*V_old_)))-(0.000000000000000e+00*exp(((-4.391000000000000e-02)*V_old_))))*(V_old_+3.778000000000000e+01))/1.000000000000000e+00)/(1.000000000000000e+00+exp((3.110000000000000e-01*(V_old_+7.923000000000000e+01))))));
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_3(){
  if((V_old_<(-4.000000000000000e+01))){
    return (((2.424000000000000e-02*exp(((-1.052000000000000e-02)*V_old_)))/(1.000000000000000e+00+exp(((-1.378000000000000e-01)*(V_old_+4.014000000000000e+01))))));
  }else{
    return (((6.000000000000000e-01*exp((5.700000000000000e-02*V_old_)))/(1.000000000000000e+00+exp(((-1.000000000000000e-01)*(V_old_+3.200000000000000e+01))))));
  }
}
double Solveode::ifnumber_4(){
  if(((time_new>=i_stimStart)&&(time_new<=i_stimEnd)&&(((time_new-i_stimStart)-(floor(((time_new-i_stimStart)/i_stimPeriod))*i_stimPeriod))<=i_stimPulseDuration))){
    return (i_stimAmplitude);
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_5(){
  if((((1.000000000000000e-02*calc_d_fCa)>0.000000000000000e+00)&&(V_old_>(-6.000000000000000e+01)))){
    return (0.000000000000000e+00);
  }else{
    return (calc_d_fCa);
  }
}
double Solveode::ifnumber_6(){
  if((calc_Ca_i<3.500000000000000e-04)){
    return ((1.000000000000000e+00/(1.000000000000000e+00+pow((calc_Ca_i/3.500000000000000e-04),6.000000000000000e+00))));
  }else{
    return ((1.000000000000000e+00/(1.000000000000000e+00+pow((calc_Ca_i/3.500000000000000e-04),1.600000000000000e+01))));
  }
}
double Solveode::ifnumber_7(){
  if((((1.000000000000000e-02*calc_d_g)>0.000000000000000e+00)&&(V_old_>(-6.000000000000000e+01)))){
    return (0.000000000000000e+00);
  }else{
    return (calc_d_g);
  }
}
double Solveode::ifnumber_8(){
  if(((len_thick/2.000000000000000e+00)<(SL_old_/2.000000000000000e+00))){
    return ((len_thick/2.000000000000000e+00));
  }else{
    return ((SL_old_/2.000000000000000e+00));
  }
}
double Solveode::ifnumber_9(){
  if((((SL_old_/2.000000000000000e+00)-(SL_old_-len_thin))>(len_hbare/2.000000000000000e+00))){
    return (((SL_old_/2.000000000000000e+00)-(SL_old_-len_thin)));
  }else{
    return ((len_hbare/2.000000000000000e+00));
  }
}
double Solveode::ifnumber_10(){
  if(((1.000000000000000e+00/calc_permtot)<1.000000000000000e+02)){
    return ((1.000000000000000e+00/calc_permtot));
  }else{
    return (1.000000000000000e+02);
  }
}
double Solveode::ifnumber_11(){
  if((xXBpostr_old_<x_0)){
    return (exp((sigmap*pow(((x_0-xXBpostr_old_)/x_0),2.000000000000000e+00))));
  }else{
    return (exp((sigman*pow(((xXBpostr_old_-x_0)/x_0),2.000000000000000e+00))));
  }
}
double Solveode::ifnumber_12(){
  if((SL_old_>SL_c)){
    return ((PCon_c*(exp((PExp_c*fabs((SL_old_-SL_c))))-1.000000000000000e+00)));
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_13(){
  if((SEon==1.000000000000000e+00)){
    return ((KSE*(SLset-SL_old_)));
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_14(){
  if(((SL_old_<=SLmax)&&(SL_old_>SLmin))){
    return SL_old_*(stretch-1)/dtime;
    //return (((intf_old_+((SLset-SL_old_)*visc))/massf));
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_15(){
  if((SL_old_<len_thick)){
    return (((-5.000000000000000e-01)*calc_dSL));
  }else{
    return (0.000000000000000e+00);
  }
}
double Solveode::ifnumber_16(){
  if((((2.000000000000000e+00*len_thin)-SL_old_)>len_hbare)){
    return (((-5.000000000000000e-01)*calc_dSL));
  }else{
    return (0.000000000000000e+00);
  }
}

static int check_flag(void *flagvalue, char *funcname, int opt){
  int *errflag;
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
    return(1);}
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      //  fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",funcname, *errflag);
      return(1); }}
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
    return(1); }
  return 0;
}

/*
static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__){
    Solveode *ode = (Solveode *) f_data__;
    ode->setVariables( 0 ,NV_Ith_S(dependent_variable__, 0));
    ode->setVariables( 1 ,NV_Ith_S(dependent_variable__, 1));
    ode->setVariables( 2 ,NV_Ith_S(dependent_variable__, 2));
    ode->setVariables( 3 ,NV_Ith_S(dependent_variable__, 3));
    ode->setVariables( 4 ,NV_Ith_S(dependent_variable__, 4));
    ode->setVariables( 5 ,NV_Ith_S(dependent_variable__, 5));
    ode->setVariables( 6 ,NV_Ith_S(dependent_variable__, 6));
    ode->setVariables( 7 ,NV_Ith_S(dependent_variable__, 7));
    ode->setVariables( 8 ,NV_Ith_S(dependent_variable__, 8));
    ode->setVariables( 9 ,NV_Ith_S(dependent_variable__, 9));
    ode->setVariables( 10 ,NV_Ith_S(dependent_variable__, 10));
    ode->setVariables( 11 ,NV_Ith_S(dependent_variable__, 11));
    ode->setVariables( 12 ,NV_Ith_S(dependent_variable__, 12));
    ode->setVariables( 13 ,NV_Ith_S(dependent_variable__, 13));
    ode->setVariables( 14 ,NV_Ith_S(dependent_variable__, 14));
    ode->setVariables( 15 ,NV_Ith_S(dependent_variable__, 15));
    ode->setVariables( 16 ,NV_Ith_S(dependent_variable__, 16));
    ode->setVariables( 17 ,NV_Ith_S(dependent_variable__, 17));
    ode->setVariables( 18 ,NV_Ith_S(dependent_variable__, 18));
    ode->setVariables( 19 ,NV_Ith_S(dependent_variable__, 19));
    ode->setVariables( 20 ,NV_Ith_S(dependent_variable__, 20));
    ode->setVariables( 21 ,NV_Ith_S(dependent_variable__, 21));
    ode->setVariables( 22 ,NV_Ith_S(dependent_variable__, 22));
    ode->setVariables( 23 ,NV_Ith_S(dependent_variable__, 23));
    ode->setVariables( 24 ,NV_Ith_S(dependent_variable__, 24));
    ode->setVariables( 25 ,NV_Ith_S(dependent_variable__, 25));
    ode->setVariables( 26 ,NV_Ith_S(dependent_variable__, 26));
    ode->setVariables( 27 ,NV_Ith_S(dependent_variable__, 27));
    ode->setParameters(0,time);
    double *t = ode->solveDiff();
    NV_Ith_S(dep_var_dot__, 0) = t[0];
    NV_Ith_S(dep_var_dot__, 1) = t[1];
    NV_Ith_S(dep_var_dot__, 2) = t[2];
    NV_Ith_S(dep_var_dot__, 3) = t[3];
    NV_Ith_S(dep_var_dot__, 4) = t[4];
    NV_Ith_S(dep_var_dot__, 5) = t[5];
    NV_Ith_S(dep_var_dot__, 6) = t[6];
    NV_Ith_S(dep_var_dot__, 7) = t[7];
    NV_Ith_S(dep_var_dot__, 8) = t[8];
    NV_Ith_S(dep_var_dot__, 9) = t[9];
    NV_Ith_S(dep_var_dot__, 10) = t[10];
    NV_Ith_S(dep_var_dot__, 11) = t[11];
    NV_Ith_S(dep_var_dot__, 12) = t[12];
    NV_Ith_S(dep_var_dot__, 13) = t[13];
    NV_Ith_S(dep_var_dot__, 14) = t[14];
    NV_Ith_S(dep_var_dot__, 15) = t[15];
    NV_Ith_S(dep_var_dot__, 16) = t[16];
    NV_Ith_S(dep_var_dot__, 17) = t[17];
    NV_Ith_S(dep_var_dot__, 18) = t[18];
    NV_Ith_S(dep_var_dot__, 19) = t[19];
    NV_Ith_S(dep_var_dot__, 20) = t[20];
    NV_Ith_S(dep_var_dot__, 21) = t[21];
    NV_Ith_S(dep_var_dot__, 22) = t[22];
    NV_Ith_S(dep_var_dot__, 23) = t[23];
    NV_Ith_S(dep_var_dot__, 24) = t[24];
    NV_Ith_S(dep_var_dot__, 25) = t[25];
    NV_Ith_S(dep_var_dot__, 26) = t[26];
    NV_Ith_S(dep_var_dot__, 27) = t[27];
    return 0;
}
 */


float __agos_factorial(int f){
  if(f>=0 & f<2)
    return 1.0;
  else if(f < 0)
    return 0.0/0.0;
  for(int i=f-1; i>=2; i--)
    f *= i;
  return (float)f;
}
