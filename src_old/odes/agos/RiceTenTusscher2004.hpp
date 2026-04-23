#ifndef RICETT2_HPP
#define RICETT2_HPP

// tentusscher2004 coupled to Rice
// Author: Bernardo Lino de Oliveira
// Date:   10/08/2010

//Cell Types:
//3 - Endo
//4 - M
//5 - Epi
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "MCutil.hpp"

#define ABSTOL 1.0E-8
#define RELTOL 1.0E-6

//static int check_flag(void *flagvalue, char *funcname, int opt);
//static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__);

class Solveode
{
  //PARAMETERS
private:
  int it_countx;
  double time; 	 // millisecond
  double R; 	 // joule_per_mole_kelvin
  double T; 	 // kelvin
  double F; 	 // coulomb_per_millimole
  double Na_o; 	 // millimolar
  double K_o; 	 // millimolar
  double P_kna; 	 // dimensionless
  double Ca_o; 	 // millimolar
  double g_K1; 	 // nanoS_per_picoF
  double g_Kr; 	 // nanoS_per_picoF
  double g_Ks; 	 // nanoS_per_picoF
  double g_Na; 	 // nanoS_per_picoF
  double g_bna; 	 // nanoS_per_picoF
  double i_stimAmplitude; 	 // picoA_per_picoF
  double i_stimStart; 	 // millisecond
  double i_stimEnd; 	 // millisecond
  double i_stimPeriod; 	 // millisecond
  double i_stimPulseDuration; 	 // millisecond
  double g_CaL; 	 // litre_per_farad_second
  double g_bca; 	 // nanoS_per_picoF
  double g_to; 	 // nanoS_per_picoF
  double P_NaK; 	 // picoA_per_picoF
  double K_mk; 	 // millimolar
  double K_mNa; 	 // millimolar
  double K_NaCa; 	 // picoA_per_picoF
  double gamma; 	 // dimensionless
  double alpha; 	 // dimensionless
  double Km_Nai; 	 // millimolar
  double Km_Ca; 	 // millimolarcalc_s_inf = (1.000000000000000e+00/(1.000000000000000e+00+exp(((V_old_+2.800000000000000e+01)/5.000000000000000e+00))));	//63
  double K_sat; 	 // dimensionless
  double g_pCa; 	 // picoA_per_picoF
  double K_pCa; 	 // millimolar
  double g_pK; 	 // nanoS_per_picoF
  double a_rel; 	 // millimolar_per_millisecond
  double b_rel; 	 // millimolar
  double c_rel; 	 // millimolar_per_millisecond
  double Vmax_up; 	 // millimolar_per_millisecond
  double K_up; 	 // millimolar
  double V_leak; 	 // per_millisecond
  double tau_g; 	 // millisecond
  double Buf_sr; 	 // millimolar
  double K_buf_sr; 	 // millimolar
  double V_c; 	 // micrometre3
  double V_sr; 	 // micrometre3
  double Cm; 	 // microF
  double K_buf_c; 	 // millimolar
  double Buf_c; 	 // millimolar
  double len_thick; 	 // micrometre
  double len_thin; 	 // micrometre
  double len_hbare; 	 // micrometre
  double kon; 	 // second_order_rate_constant
  double Qkon; 	 // dimensionless
  double TmpC; 	 // celsius
  double koffL; 	 // first_order_rate_constant
  double koffmod; 	 // dimensionless
  double Qkoff; 	 // dimensionless
  double koffH; 	 // first_order_rate_constant
  double kn_p; 	 // first_order_rate_constant
  double Qkn_p; 	 // dimensionless
  double kp_n; 	 // first_order_rate_constant
  double Qkp_n; 	 // dimensionless
  double perm50; 	 // dimensionless
  double nperm; 	 // dimensionless
  double fapp; 	 // first_order_rate_constant
  double xbmodsp; 	 // dimensionless
  double Qfapp; 	 // dimensionless
  double gslmod; 	 // dimensionless
  double gapp; 	 // first_order_rate_constant
  double Qgapp; 	 // dimensionless
  double hfmdc; 	 // dimensionless
  double x_0; 	 // micrometre
  double hbmdc; 	 // dimensionless
  double hf; 	 // first_order_rate_constant
  double Qhf; 	 // dimensionless
  double hb; 	 // first_order_rate_constant
  double Qhb; 	 // dimensionless
  double sigmap; 	 // dimensionless
  double sigman; 	 // dimensionless
  double gxb; 	 // first_order_rate_constant
  double Qgxb; 	 // dimensionless
  double xPsi; 	 // dimensionless
  double kxb; 	 // millinewton_per_millimetre2
  double SLrest; 	 // micrometre
  double PCon_t; 	 // unit_normalised_force
  double PExp_t; 	 // per_micrometre
  double PCon_c; 	 // unit_normalised_force
  double PExp_c; 	 // per_micrometre
  double SL_c; 	 // micrometre
  double SLset; 	 // micrometre
  double KSE; 	 // unit_normalised_force_per_micrometre
  double SEon; 	 // dimensionless
  double visc; 	 // unit_normalised_force_millisecond_per_micrometre
  double massf; 	 // unit_normalised_force_millisecond2_per_micrometre
  double SLmax; 	 // micrometre
  double SLmin; 	 // micrometre
  double Trop_conc; 	 // micromolar
  double calc_E_Na; 	 // millivolt
  double calc_E_K; 	 // millivolt
  double calc_E_Ks; 	 // millivolt
  double calc_E_Ca; 	 // millivolt
  double calc_alpha_K1; 	 // dimensionless
  double calc_beta_K1; 	 // dimensionless
  double calc_xK1_inf; 	 // dimensionless
  double calc_i_K1; 	 // picoA_per_picoF
  double calc_i_Kr; 	 // picoA_per_picoF
  double calc_xr1_inf; 	 // dimensionless
  double calc_alpha_xr1; 	 // dimensionless
  double calc_beta_xr1; 	 // dimensionless
  double calc_tau_xr1; 	 // millisecond
  double calc_xr2_inf; 	 // dimensionless
  double calc_alpha_xr2; 	 // dimensionless
  double calc_beta_xr2; 	 // dimensionless
  double calc_tau_xr2; 	 // millisecond
  double calc_i_Ks; 	 // picoA_per_picoF
  double calc_xs_inf; 	 // dimensionless
  double calc_alpha_xs; 	 // dimensionless
  double calc_beta_xs; 	 // dimensionless
  double calc_tau_xs; 	 // millisecond
  double calc_i_Na; 	 // picoA_per_picoF
  double calc_m_inf; 	 // dimensionless
  double calc_alpha_m; 	 // dimensionless
  double calc_beta_m; 	 // dimensionless
  double calc_tau_m; 	 // millisecond
  double calc_h_inf; 	 // dimensionless
  double calc_alpha_h; 	 // per_millisecond
  double calc_beta_h; 	 // per_millisecond
  double calc_tau_h; 	 // millisecond
  double calc_j_inf; 	 // dimensionless
  double calc_alpha_j; 	 // per_millisecond
  double calc_beta_j; 	 // per_millisecond
  double calc_tau_j; 	 // millisecond
  double calc_i_b_Na; 	 // picoA_per_picoF
  double calc_i_stim; 	 // picoA_per_picoF
  double calc_i_CaL; 	 // picoA_per_picoF
  double calc_d_inf; 	 // dimensionless/state/partition2/belino/simulador_assincrono_backup_12_11/tissue_test.stl

  double calc_alpha_d; 	 // dimensionless
  double calc_beta_d; 	 // dimensionless
  double calc_gamma_d; 	 // millisecond
  double calc_tau_d; 	 // millisecond
  double calc_f_inf; 	 // dimensionless
  double calc_tau_f; 	 // millisecond
  double calc_alpha_fCa; 	 // dimensionless
  double calc_beta_fCa; 	 // dimensionless
  double calc_gama_fCa; 	 // dimensionless
  double calc_fCa_inf; 	 // dimensionless
  double calc_tau_fCa; 	 // millisecond
  double calc_d_fCa; 	 // per_millisecond
  double calc_i_b_Ca; 	 // picoA_per_picoF
  double calc_i_to; 	 // picoA_per_picoF
  double calc_s_inf; 	 // dimensionless
  double calc_tau_s; 	 // millisecond
  double calc_r_inf; 	 // dimensionless
  double calc_tau_r; 	 // millisecond
  double calc_i_NaK; 	 // picoA_per_picoF
  double calc_i_NaCa; 	 // picoA_per_picoF
  double calc_i_p_Ca; 	 // picoA_per_picoF
  double calc_i_p_K; 	 // picoA_per_picoF
  double calc_i_rel; 	 // millimolar_per_millisecond
  double calc_i_up; 	 // millimolar_per_millisecond
  double calc_i_leak; 	 // millimolar_per_millisecond
  double calc_g_inf; 	 // dimensionless
  double calc_d_g; 	 // per_millisecond
  double calc_Ca_sr_bufsr; 	 // dimensionless
  double calc_cc; 	 // (null)
  double calc_bc; 	 // (null)
  double calc_Ca_i; 	 // millimolar
  double calc_Cai; 	 // micromolar
  double calc_sovr_ze; 	 // micrometre
  double calc_sovr_cle; 	 // micrometre
  double calc_len_sovr; 	 // micrometre
  double calc_SOVFThick; 	 // dimensionless
  double calc_SOVFThin; 	 // dimensionless
  double calc_konT; 	 // second_order_rate_constant
  double calc_koffLT; 	 // first_order_rate_constant
  double calc_koffHT; 	 // first_order_rate_constant
  double calc_kn_pT; 	 // first_order_rate_constant
  double calc_kp_nT; 	 // first_order_rate_constant
  double calc_permtot; 	 // dimensionless
  double calc_inprmt; 	 // dimensionless
  double calc_dTRPNCaL; 	 // first_order_rate_constant
  double calc_dTRPNCaH; 	 // first_order_rate_constant
  double calc_Tropreg; 	 // dimensionless
  double calc_fappT; 	 // first_order_rate_constant
  double calc_gapslmd; 	 // dimensionless
  double calc_gappT; 	 // first_order_rate_constant
  double calc_hfmd; 	 // dimensionless
  double calc_hbmd; 	 // dimensionless
  double calc_hfT; 	 // first_order_rate_constant
  double calc_hbT; 	 // first_order_rate_constant
  double calc_gxbmd; 	 // dimensionless
  double calc_gxbT; 	 // first_order_rate_constant
  double calc_dXBprer; 	 // first_order_rate_constant
  double calc_dXBpostr; 	 // first_order_rate_constant
  double calc_P; 	 // dimensionless
  double calc_dxXBprer; 	 // micrometre_per_millisecond
  double calc_dxXBpostr; 	 // micrometre_per_millisecond
  double calc_dutyprer; 	 // dimensionless
  double calc_dutypostr; 	 // dimensionless
  double calc_SSXBprer; 	 // dimensionless
  double calc_SSXBpostr; 	 // dimensionless
  double calc_Fnordv; 	 // millinewton_micrometre_per_millimetre2
  double calc_force; 	 // millinewton_micrometre_per_millimetre2
  double calc_active; 	 // unit_normalised_force
  double calc_ppforce_t; 	 // unit_normalised_force
  double calc_ppforce_c; 	 // unit_normalised_force
  double calc_ppforce; 	 // unit_normalised_force
  double calc_preload; 	 // unit_normalised_force
  double calc_afterload; 	 // unit_normalised_force
  double calc_dSL; 	 // micrometre_per_millisecond
  double calc_FrSBXB; 	 // dimensionless
  double calc_dFrSBXB; 	 // first_order_rate_constant
  double calc_dsovr_ze; 	 // micrometre_per_millisecond
  double calc_dsovr_cle; 	 // micrometre_per_millisecond
  double calc_dlen_sovr; 	 // micrometre_per_millisecond
  double calc_dSOVFThin; 	 // first_order_rate_constant
  double calc_dSOVFThick; 	 // first_order_rate_constant
  double calc_TropTot; 	 // micromolar
  double calc_dTropTot; 	 // micromolar_per_millisecond
  double calc_dforce; 	 // millinewton_micrometre_per_millimetre2_per_millise
  double dtime, *time_vec__;
  double time_new;

  double stretch;
  int cell_type;

  //DEPENDENT VARIABLES
  double *Xr1;
  double Xr1_new_, Xr1_old_, Xr1_ini_;
  double *Xr2;
  double Xr2_new_, Xr2_old_, Xr2_ini_;
  double *Xs;
  double Xs_new_, Xs_old_, Xs_ini_;
  double *m;
  double m_new_, m_old_, m_ini_;
  double *h;
  double h_new_, h_old_, h_ini_;
  double *j;
  double j_new_, j_old_, j_ini_;
  double *V;
  double V_new_, V_old_, V_ini_;
  double *d;
  double d_new_, d_old_, d_ini_;
  double *f;
  double f_new_, f_old_, f_ini_;
  double *fCa;
  double fCa_new_, fCa_old_, fCa_ini_;
  double *s;
  double s_new_, s_old_, s_ini_;
  double *r;
  double r_new_, r_old_, r_ini_;
  double *g;
  double g_new_, g_old_, g_ini_;
  double *Ca_SR;
  double Ca_SR_new_, Ca_SR_old_, Ca_SR_ini_;
  double *Ca_tot;
  double Ca_tot_new_, Ca_tot_old_, Ca_tot_ini_;
  double *Na_i;
  double Na_i_new_, Na_i_old_, Na_i_ini_;
  double *K_i;
  double K_i_new_, K_i_old_, K_i_ini_;
  double *TRPNCaL;
  double TRPNCaL_new_, TRPNCaL_old_, TRPNCaL_ini_;
  double *TRPNCaH;
  double TRPNCaH_new_, TRPNCaH_old_, TRPNCaH_ini_;
  double *N_NoXB;
  double N_NoXB_new_, N_NoXB_old_, N_NoXB_ini_;
  double *P_NoXB;
  double P_NoXB_new_, P_NoXB_old_, P_NoXB_ini_;
  double *N;
  double N_new_, N_old_, N_ini_;
  double *XBprer;
  double XBprer_new_, XBprer_old_, XBprer_ini_;
  double *XBpostr;
  double XBpostr_new_, XBpostr_old_, XBpostr_ini_;
  double *xXBprer;
  double xXBprer_new_, xXBprer_old_, xXBprer_ini_;
  double *xXBpostr;
  double xXBpostr_new_, xXBpostr_old_, xXBpostr_ini_;
  double *SL;
  double SL_new_, SL_old_, SL_ini_;
  double *intf;
  double intf_new_, intf_old_, intf_ini_;

  double *Ca_i;
  double *active_force;

  // CVODE VARIABLES
  /*
	realtype reltol__;
	void *cvode_mem_cvode__;
	N_Vector dependent_variable__, abstol__;
	int flag__, flagr__;
	double *depvar__;
   */

public:
  Solveode(int type);
  //Solveode(int type, double val_abstol__ = ABSTOL, double val_reltol__= RELTOL, int method__ = CV_BDF);
  ~Solveode();
  int setVariables(int, double);
  int setParameters(int, double);
  int setFreeVariable(double);
  void setAbstol(int index, double value);
  void setReltol(double value);
  //void reInitCVODE();
  //void setCVODEMaxStep(double maxstep);
  double getVariables(int);
  double getParameters(int);
  double getFreeVariable();
  Variables get_Parameters();
  Variables get_Variables();
  Variables get_FreeVariable();
  void setParametersFromFile(char*);
  void setVariablesFromFile(char*);
  void setFreeVariableFromFile(char*);
  double* solveDiff();
  //int solveCVODE(int firstcall__ = 0, int steps__ = 0);
  double solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0);
  double solve(int firstcall__, int num_iterations, int num_results__, double stim, double start, double end);
  double* getSolution(int indVariable);
  double* getIndependentVar();
  double solveToFile(char *filename, char *fileaccess = "", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);
  //void solveCVODEToFile(char *filename, char *fileaccess = "", int firstcall__ = 0, int steps__ = 0);
private:
  inline double ifnumber_0();
  inline double ifnumber_1();
  inline double ifnumber_2();
  inline double ifnumber_3();
  inline double ifnumber_4();
  inline double ifnumber_5();
  inline double ifnumber_6();
  inline double ifnumber_7();
  inline double ifnumber_8();
  inline double ifnumber_9();
  inline double ifnumber_10();
  inline double ifnumber_11();
  inline double ifnumber_12();
  inline double ifnumber_13();
  inline double ifnumber_14();
  inline double ifnumber_15();
  inline double ifnumber_16();
};

#endif
