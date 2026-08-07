#include "kerkoff2003.hpp"

Kerkoff2003::Kerkoff2003() : CellModel(2)
{
  var_names.insert( std::pair<int, std::string>(0, "lc") );
  
  monitored.push_back( &active_stress);
  active_stress = 0.0; 
}

void Kerkoff2003::init(double * values) const
{  
  assert(values != 0);
  values[0] = 1.9;
  values[1] = 0.136; 
}

void Kerkoff2003::equation(const double time,
				 const double * statevars,
				 double * values)
{
  // State variables
  const double lc  = statevars[0];
  const double lat = statevars[1];

  // Some constants and parameters
  const double a6 = 2.0;
  const double a7 = 1.5;
  const double T0 = 180.0;
  const double Ea = 20.0;
  const double v0 = 7.5;
  const double ls0 = 1.9;
  const double tr = 0.075;
  const double td = 0.075;
  const double b = 0.21; 
  const double ld = -0.4;  

  // Calculations 
  values[0] = ((Ea * (ls0 - lc) - 1.0) * v0);
  values[1] = 0.0; //the lat is constant

  double _time = time-lat; //
  double tmax = b * (ls0 - ld);
  if (_time < 0.0 || _time > tmax)
    active_stress = 0.0;
  else
    active_stress = ((lc <= a7) ? 0.0 : T0 * pow(tanh(a6 * (lc - a7)), 2)) * pow(tanh(_time / tr), 2) * pow(tanh((tmax - _time) / td), 2) * (ls0-lc) *Ea;

}
