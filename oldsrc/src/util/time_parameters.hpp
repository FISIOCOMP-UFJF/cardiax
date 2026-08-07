#ifndef TIME_PARAMETERS_HPP
#define TIME_PARAMETERS_HPP

#include <iostream>
#include <cmath>

using namespace std;

/** 
    Simple class to handle time parameters for use in time discretization 
*/
class TimeParameters
{
public:

  //! Default constructor
  TimeParameters() {};

  //! Constructor
  TimeParameters(double dt, double tend, double tr) : 
    iteration(0), start(0.0), stop(tend), 
    time_step(dt), time_rate(tr), cur_time(0.0)
  {
    // compute some values
    nsteps = (int) ((stop-start)/time_step);
    
    print_rate = (int) ceil((time_rate/time_step));
    size = (nsteps/print_rate);

    // print information
    info();
  };

  ~TimeParameters(){}

  inline int it()          { return iteration;  }
  inline int pr()          { return print_rate; }
  inline int get_size()    { return size;       }
  inline int get_nsteps()  { return nsteps;     }
  inline double get_dt()   { return time_step;  }
  inline double get_stop() { return stop;       }

  inline bool finished ()
  {
    return (cur_time < stop) ? false : true;
  }

  inline void increase_time ()
  {
    iteration++;
    cur_time = iteration * time_step;
  }

  void info();

  void parse_xml(std::string filename, std::string nodename);

  inline void reset()
  {
    iteration = 0;
    cur_time = 0.0;
  }

  inline void set_time(double a)
  {
    cur_time = cur_time * a;
  }

  inline void show_time()
  {
    if (time2print())
      cout << "Time  " << scientific << cur_time << " ";
  }

  inline double time()
  {
    return cur_time;
  }

  inline bool time2print()
  {
    return (iteration % print_rate == 0) ? true : false;
  }

private:

  //! Number of the current iteration
  int iteration;

  //! Total number of time steps
  int nsteps;

  //! Output time rate in iterations (e.g. each 10 iterations)
  int print_rate;

  //! Num time steps divided by print rate
  int size;

  //! Time to start simulation
  double start;

  //! Time to stop simulation (final time)
  double stop;

  //! Time step
  double time_step;

  //! Output time rate (e.g. each 1 ms)
  double time_rate;

  //! Current time
  double cur_time;

};

#endif
