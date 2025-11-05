#include "timer.hpp"
#include <list>

Timer::Timer() : num_calls(0), total_time(0.0)
{
  // do nothing
}

void Timer::start()
{
  num_calls++;
  gettimeofday(&start_time, NULL);
}

void Timer::reset()
{
  total_time = 0.0;
  num_calls = 0;
}

void Timer::stop()
{
  gettimeofday(&stop_time, NULL);

  total_time += stop_time.tv_sec - start_time.tv_sec;
  total_time += (stop_time.tv_usec - start_time.tv_usec)/1000000.0;

  start_time = stop_time;
}

double Timer::get_total_time()
{
  return total_time;
}

int Timer::get_num_calls()
{
  return num_calls;
}

// -----------------------------------------------------------------------------

TimerSection::TimerSection()
{
  // do nothing
}

TimerSection::~TimerSection()
{
  while (active_sections.size() > 0)
    leave();

  //print_summary();
}

void TimerSection::enter(const std::string & section_name)
{
  assert(section_name.empty() == false);

  if(std::find(active_sections.begin(), active_sections.end(),
	       section_name) != active_sections.end())
  {
    std::cout << "Error in TimerSection: section " << section_name
	      << " is already active." << std::endl;
    exit(0);
  }

  if (sections.find(section_name) == sections.end())
  {
    sections[section_name].n_calls = 0;
    sections[section_name].time = 0;
  }

  sections[section_name].timer.reset();
  sections[section_name].timer.start();
  sections[section_name].n_calls++;

  active_sections.push_back(section_name);
}

double TimerSection::time(const std::string & name)
{
  // check for a section that was never created
  assert( sections.find(name) != sections.end() );
  return sections[name].time;
}
	    
void TimerSection::leave(const std::string & section_name)
{
  assert( !active_sections.empty() );

  if (section_name != "")
  {
    // check for a section that was never created
    assert( sections.find(section_name) != sections.end() );
    // check for a section that has not been entered
    assert( std::find(active_sections.begin(), active_sections.end(),
		      section_name) != active_sections.end());
  }

  // if string is empty, exit the last active
  const std::string actual_name = (section_name == "" ?
				   active_sections.back() :
				   section_name);

  sections[actual_name].timer.stop();
  sections[actual_name].time += sections[actual_name].timer.get_total_time();

  // delete from the list of active sections
  
  active_sections.erase(std::find(active_sections.begin(),
				  active_sections.end(),
				  actual_name));
}

void TimerSection::summary()
{
  double total = 0;

  for (std::map<std::string, Section>::const_iterator
	 i = sections.begin(); i!=sections.end(); ++i)
  {
    total += i->second.time;
  }

  // print a nice header
  std::cout << "\n"
       << "+----------------------------------------+------------+------------+\n"
       << "| Total time elapsed                     |            |            |\n"
       << "|                                        |            |            |\n"
       << "| Section                    | no. calls | time (sec) | % of total |\n"
       << "+----------------------------+-----------+------------+------------+\n";

  for (std::map<std::string, Section>::const_iterator
	 i = sections.begin(); i!=sections.end(); ++i)
  {
    std::string name_out = i->first;
    uint size = name_out.size();
    uint fill = 29 - size;
    std::cout << "| " << name_out;
    std::cout << std::setw(fill);
    std::cout << "| ";
    std::cout << std::setw(9);
    std::cout << i->second.n_calls << " |";
    std::cout << std::setw(10);
    std::cout << std::setprecision(3);
    std::cout << i->second.time << "s | ";
    std::cout << std::setw(4);
    std::cout << std::setprecision(2);    
    double value = i->second.time/total * 100;
    std::cout << value << " % |\n";
  }

  // print a footer
  std::cout << "+----------------------------+-----------+------------+------------+"
	    << "\n\n";
  
}
