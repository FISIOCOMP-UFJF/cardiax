#include "fitz_hugh_nagumo.hpp"

FitzHughNagumo::FitzHughNagumo() : CellModel(2)
{
  var_names.insert( std::pair<int, std::string>(0, "v") );
  var_names.insert( std::pair<int, std::string>(1, "w") );
}

void FitzHughNagumo::init(double * values) const
{  
  assert(values != 0);
  values[0] = -1.19940803524;
  values[1] = -0.624260044055;
}

void FitzHughNagumo::equation(const double time,
			      const double * statevars,
			      double * values)
{
  // State variables
  const double v = statevars[0];
  const double w = statevars[1];

  // Some constants and parameters
  const double eps   = 0.2;
  const double gamma = 0.8;
  const double beta  = 0.7;

  // Calculations
  double v3 = pow(v,3);
  double dVdt = ((1.0/eps) * (v-((1.0/3.0)*v3)-w)) + i_stim;
  double dWdt = eps * (v - gamma*w + beta);

  values[0] = dVdt;
  values[1] = dWdt;
}

