#include "nash_panfilov.hpp"

NashPanfilov::NashPanfilov() : CellModel(3)
{
  var_names.insert( std::pair<int, std::string>(0, "Vm") );
  var_names.insert( std::pair<int, std::string>(1, "r") );
  var_names.insert( std::pair<int, std::string>(2, "Ta") );
}

void NashPanfilov::init(double * values) const
{
  assert(values != 0);
  values[0] = 0.0; //-80; //0.0
  values[1] = 0.1;
  values[2] = 0.0;
}

void NashPanfilov::equation(const double time, const double *rY, double *rDY)
{
  // State variables
  const double V  = rY[0];
  const double r  = rY[1];
  const double Ta = rY[2];
 
  // Parameters
  const double k   = 8.0;     
  const double a   = 0.01;    
  const double b   = 0.15;
  const double eps = 0.002; 
  const double mu1 = 0.2;
  const double mu2 = 0.3;
  const double kTa = 13.7; // 47.9 ou 10 kPa
    
  // Calculations
  double d_dt_v  = k*V*(V-a)*(1.0-V) - r*V + i_stim;
  double d_dt_r  = (eps + (mu1 * r)/(mu2 + V)) * (-r - k*V*(V-b-1.0));
  double d_dt_Ta = rateswitchEr(V) * (kTa * (V) - Ta);
  
  rDY[0] = d_dt_v;
  rDY[1] = d_dt_r;
  rDY[2] = d_dt_Ta;
}

double NashPanfilov::rateswitchNP(double v)
{
  if (v < 0.05)
    return 1.0;
  else
    return 0.1;
}

double NashPanfilov::rateswitchGK(double v)
{
  // My adjustments
  const double e0   = 0.1;
  const double einf = 1.0;
  const double psi  = 1.0;
  const double vbar = 0.0;
  return e0 + (einf-e0) * exp(-exp(-psi*(v-vbar))); 
}

double NashPanfilov::rateswitchEr(double v)
{
  const double e0   = 1.0;
  const double einf = 0.1;
  const double psi  = 0.3;
  const double vbar = 0.2;
  return e0 + (einf-e0)*exp(-exp(-psi*(v-vbar)));
}