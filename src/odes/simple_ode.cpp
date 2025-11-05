#include "simple_ode.hpp"

SimpleODE::SimpleODE() : CellModel(1)
{
  var_names.insert( std::pair<int, std::string>(0, "variable_v") );
}

void SimpleODE::init(double * values) const
{
  values[0] = 5.0;
}

void SimpleODE::equation(const double t, const double * statevars,
			 double * values)
{
  // State variable
  const double u = statevars[0];
  
  // Calculations
  values[0] = -2.0 * u;
}

void SimpleODE::dump_exact()
{
  ofstream outfile("exact");
  int N = ts->number_of_time_steps();
  double t = 0.0;
  double step = ts->timestep();
  double * y  = new double[N];
  
  y[0] = 5.0;
  outfile << t << "\t" << y[0] << endl;

  for(int i=1; i<=N; i++)
  {
    t = i*step;
    y[i] = 5.0 * exp(-2.0*t);
    outfile << t << "\t" << y[i] << endl;
  }    
  
  outfile.close();
  delete [] y;
}
