#ifndef TIMER_H
#define TIMER_H

#include <list>
#include <algorithm>
#include <string>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <map>
/*!
 *  A simple class for measuring time in seconds.
 *  
 *   TODO: it is said that gettimeofday is not good to measure time.
 *   Need to improve this!
 */

class Timer
{
public:

  //! Constructor
  Timer();

  //! Start to measure time
  void start();

  //! Stop to measure and accumulate to total time
  void stop();

  //! Zero the total time
  void reset();

  //! Returns the accumulated total time
  double get_total_time();

  //! Returns the number of calls to measure time
  int get_num_calls();

private:

  int num_calls;
  double total_time;
  struct timeval start_time;
  struct timeval stop_time;

};

/*!
 *  Class for measuring time (in seconds) and number of calls
 *   of several sections of the code. 
 *
 *   Example of usage:
 * 
 *      timer.enter("Assembly");
 *      assemble_stiffness();
 *      timer.leave();
 *
 *      timer.enter("Solution");
 *      solve(A,x,b,max,tol);
 *      timer.leave();
 *
 *      timer.summary();
*/

class TimerSection
{
public:

  //! Constructor
  TimerSection();

  //! Destructor
  ~TimerSection();

  //! Enter section
  void enter(const std::string & section_name);
  
  //! Return the current time of a section
  double time(const std::string & section_name);
	    
  //! Leave the last entered section
  void leave(const std::string & section_name = std::string());

  //! Print summary of information
  void summary();

private:
  
  //! Overall time
  Timer total_time;

  //! A simple data structure to hold information of one section.
  struct Section
  {
    Timer timer;
    uint n_calls;
    double time;
  };

  //! A list of all the sections and their information.
  std::map<std::string, Section> sections;

  //! A list of the sections that have been entered and not exited. 
  //! The list is kept in the order in which sections have been entered.
  std::list<std::string> active_sections;


};

#endif
