#include "timestepper.h"

TimeStepper::TimeStepper() 
  : it(0), start(0.0), end(0.0), dt(0), time_rate(0)
{
  nsteps = 0;
  init();
  info();
}

TimeStepper::TimeStepper(double timestep, double tend)
  : it(0), start(0.0), end(tend), dt(timestep), time_rate(1.0)
{
  if (end>start)
  {
    nsteps = (int) ((end-start)/dt);
    print_rate = (int) ceil((time_rate/dt));
    size = (nsteps/print_rate);
  }
  else
    nsteps=0;

  init();
  info();
}

TimeStepper::TimeStepper(double timestep, double tend, double tr)
  : it(0), start(0.0), end(tend), dt(timestep), time_rate(tr)
{
  if (end>start)
  {
    nsteps = (int) ((end-start)/dt);
    print_rate = (int) ceil((time_rate/dt));
    size = (nsteps/print_rate);
  }
  else
    nsteps=0;

  init();
  info();
}

void TimeStepper::info()
{
  cout << std::setprecision(3);
  cout << fixed;
#ifdef DEBUG
  cout << " Solving ODE on time interval [0.0," << end << "]" << endl;
  cout << " with step size " << dt << endl;
#endif
}

void TimeStepper::init ()
{
  it = 0;
  cur_time = 0.0;
}

void TimeStepper::increase_time ()
{
  it++;
  cur_time = it * dt;
}

double TimeStepper::timestep() const
{
  return dt; 
}

int TimeStepper::number_of_time_steps() const
{
  return nsteps; 
}

double TimeStepper::time() const
{
  return cur_time; 
}

double TimeStepper::starttime() const
{
  return start;
}

double TimeStepper::endtime() const
{
  return end; 
}

