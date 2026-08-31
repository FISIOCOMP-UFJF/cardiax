#ifndef RICE_TEN_TUSSCHER_HPP
#define	RICE_TEN_TUSSCHER_HPP

#include <cmath>
#include "cellmodel.hpp"
#include "ode_solver.hpp"

///////////////////////////////////#define HF

using namespace std;

/*! Coupled electromechanical cell model developed by B. L. de Oliveira.
    The model coupled the TNNP electrophysiological cell model to the
    Rice et al. myofilament model of contraction.

    References: 
    - "Approximate Model of Cooperative Activation and Crossbridge Cycling
       in Cardiac Muscle Using Ordinary Differential Equations"
       John Jeremy Rice, Fei Wang, Donald M. Bers and Pieter P. de Tombe

    - "A model for human ventricular tissue"
       K. H. W. J. ten Tusscher, D. Noble, P. J. Noble and A. V. Panfilov
       Am J Physiol Heart Circ Physiol 286, H1573-H1589, 2004
*/

class RiceTenTusscher : public CellModel
{

public:

  //! Default constructor (initializes constants ans some other vars)
  RiceTenTusscher();

  //! Set initial conditions for RiceTT2 model
  virtual void init(double * values) const;

  //! Compute RHS equations for the RiceTT2 model
  virtual void equation(const double t, const double * sv, double * values);

  //! TNNP vars
  double V, Xr1, Xr2, Xs, m, h, j, d, f, fCa, s, r, g, Ca_SR, Ca_i_total, Na_i, K_i;
  
  //! Rice vars
  double TRPNCaL, TRPNCaH, N_NoXB, P_NoXB, N, XBprer, XBpostr;
  double xXBprer, xXBpostr, SL, intf;

  //! Rate variables
  double d_dt_V, d_dt_Xr1, d_dt_Xr2, d_dt_Xs, d_dt_m, d_dt_h, d_dt_j;
  double d_dt_d, d_dt_f, d_dt_fCa, d_dt_s, d_dt_r, d_dt_g, d_dt_Ca_SR;
  double d_dt_Ca_i_total, d_dt_Na_i, d_dt_K_i, d_dt_N_NoXB, d_dt_P_NoXB, d_dt_N;

  //! Is isometric twitch?
  bool isometric;
        
private:   
      
  //! Unit normalised force    
  double active;
  
  //! Some variables
  double TropTot;
  double Ca_i_bufc; // millimolar
  double Ca_sr_bufsr; // dimensionless
  double Ca_i; // millimolar
  double dTRPNCaL; // first_order_rate_constant
  double dTRPNCaH; // first_order_rate_constant  
  double dSL; // micrometre_per_millisecond
  double Cai; // micromolar
  
  //! Currents in picoA_per_picoF
  double IbNa, ICaL, IKs, INa, IK1, IKr, IbCa, Ito, INaK, INaCa, IpCa, IpK, Istim;
  
  //! Currents in millimolar_per_millisecond
  double Irel, Iup, Ileak;  

};

#endif

