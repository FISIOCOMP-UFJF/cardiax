#include <stdexcept>
#include <string.h>
#include "time_parameters.hpp"
#include "util/pugixml.hpp"

void TimeParameters::info()
{
  cout << "Time discretization information" << endl;
  cout << " Number of time steps: " << nsteps << endl;
  cout << fixed;
  cout << " Time step: " << time_step << endl;
  cout << " Start time: " << start << endl;
  cout << " Final time: " << stop << endl;
}

void TimeParameters::parse_xml(std::string filename, std::string nodename)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(filename.c_str());

  if (!result)
  {
    std::cout << "XML parsed with errors, attr value: ["
    << doc.child("mesh").attribute("attr").value() << "]\n";
    std::cout << "Error description: " << result.description() << "\n";
    std::cout << "Error offset: " << result.offset;
    std::cout << " (error at [..." << (result + result.offset) << "]\n\n";
    std::runtime_error("XML parsing error");
  }

  pugi::xml_node p = doc.child(nodename.c_str()).child("time_parameters");

  for (pugi::xml_node_iterator it = p.begin(); it != p.end(); ++it)
  {
    //double dt, double tend, double tr)
    if( strcmp(it->name(), "time_step") == 0  )
      time_step = it->text().as_double();

    if( strcmp(it->name(), "total_time") == 0 )
      stop = it->text().as_double();

    if( strcmp(it->name(), "time_rate" ) == 0 )
      time_rate = it->text().as_double();
  }
}