#include "mitchell_schaeffer.hpp"

MitchellSchaeffer::MitchellSchaeffer() : CellModel(3)
{
  var_names.insert( std::pair<int, std::string>(0, "V") );
  var_names.insert( std::pair<int, std::string>(1, "h") );
  var_names.insert( std::pair<int, std::string>(2, "Ta") );
}

void MitchellSchaeffer::init(double * values) const
{  
  assert(values != 0);
  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 0.0;
}

void MitchellSchaeffer::equation(const double time,
				 const double * statevars,
				 double * values)
{
  // State variables
  const double v  = statevars[0];
  const double h  = statevars[1];
  const double Ta = statevars[2];
  const double kTa = 1.0; 

  // Some constants and parameters
  const double v_gate    = 0.13;
  const double tau_in    = 0.3; 
  const double tau_out   = 6.0; 
  const double tau_open  = 120.0; 
  const double tau_close = 150.0; 

  // Calculations
  double j_in   = (h * (v*v*(1.-v)))/tau_in;
  double j_out  = -v/tau_out;
  double j_stim = i_stim;

  values[0] = j_in + j_out + j_stim;
  values[1] = (v < v_gate) ? ((1.-h)/tau_open) : (-h/tau_close);
  values[2] = rateswitchGK(v, h) * (kTa * v - 0.8*Ta);
}

double MitchellSchaeffer::rateswitchGK2(double v)
{
  // My adjustments
  const double e0   = 0.1;
  const double einf = 1.0;
  const double psi  = 1.0;
  const double vbar = 0.0;
  return e0 + (einf-e0) * exp(-exp(-psi*(v-vbar))); 
}

double MitchellSchaeffer::rateswitchGK(double v, double h)
{
  if(h < 0.9)
    return 0.02;
  else if (h >= 0.9 && v < 0.01)
    return 0.02; 
  else
    return 0.0;
}
