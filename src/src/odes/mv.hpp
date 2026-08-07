#ifndef MV_HPP
#define	MV_HPP

#include "cellmodel.hpp"
#include "ode_solver.hpp"

class MinimalVentricular  : public CellModel
{
public:
  
  //! Default constructor
  MinimalVentricular();
    
  //! Set initial conditions for LRI model
  virtual void init(double * values) const; 
  
  //! Compute RHS equations for the LRI model
  virtual void equation(const double time, const double * sv, double * values);
  
  //! Initialize constants and parameters for cell type
  void init_params();
  
private:
  
  bool paramset;
  double uo, uu, thetav, thetaw, thetavb, thetao;
  double tauv1b, tauv2b, tauvp, tauw1b, tauw2b, kwb, uwb;
  double tauwp, tau_fi, tauo1, tauo2, tauso1, tauso2, kso, uso;
  double taus1, taus2, ks, us, tau_si, tau_winf, w_infs;
    
};


#endif	/* MV_HPP */

