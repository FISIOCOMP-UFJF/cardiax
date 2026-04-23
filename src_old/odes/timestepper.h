#ifndef TIME_PRM_H
#define TIME_PRM_H

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/** March in time to compute the solution of an ODE problem */
class TimeStepper
{
 private:

  //! Time step counter
  int it;

  //! Number of time steps
  int  nsteps;

  //! Start time
  double start;

  //! End time
  double end;

  //! Time step size
  double dt;

  //! Current time of the time stepping
  double cur_time;

  //! Time to print
  double time_rate;

  //! Time to save
  double tts;
  
  //! Print rate in iterations
  int  print_rate;

  //! Dont remember exactly!
  int size;

public:

  //! Default constructor
  TimeStepper();

  //! Create a time stepper with given time step and duration
  TimeStepper(double timestep, double tend);

  //! Create a time stepper with given time step, duration and time to print
  TimeStepper(double timestep, double tend, double tr);

  //! Default destructor
  ~TimeStepper(){}
  
  //! Print information about time stepping
  void info();

  //! Prepare to loop in time
  void init ();

  //! Increment current time
  void increase_time ();

  //! Return time step
  double timestep() const;

  //! Return the current time value
  double time() const;

  //! Return the start time [t0,t1] -> t0
  double starttime() const;

  //! Return the end time [t0,t1] -> t1
  double endtime() const;

  //! Return the nsteps
  int number_of_time_steps() const;

  //! Check if it has finished
  inline bool finished() {return (cur_time>=end) ? true : false;}

  //! Check if it is time to print or save to file
  inline bool time_to_print() {return (it%print_rate==0) ? true : false;}

  //! Check if it is time to save to file
  inline double time_to_save() const { return time_rate; }

};

#endif
