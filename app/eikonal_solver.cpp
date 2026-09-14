#include <iostream>
#include <string>
#include "util/toml.hpp"
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "pdes/eikonal.hpp"
#include "mesh/mesh.hpp"

using namespace std;

void usage()
{
  std::cout << std::endl;
  std::cout << " Usage: eikonal -f config.toml\n";
  std::cout << std::endl;
  exit(0);
}

void print_banner()
{
  cout << "========================================" << endl;
  cout << "     E I K O N A L   S O L V E R        " << endl;
  cout << "========================================" << endl;
}

int main(int argc, const char *argv[])
{ 
  CommandLineArgs::init(argc, argv);
  string inputfile = CommandLineArgs::read("-f", "");   // empty = not provided
  if (inputfile.empty())
  {
      usage();
      return 1;
  }

  // read the TOML config
  toml::table cfg;
  try
  {
    cfg = toml::parse_file(inputfile);
  }
  catch (const toml::parse_error & err)
  {
    cerr << " *** ERROR parsing '" << inputfile << "':\n" << err << endl;
    return 1;
  }

  print_banner();

  Eikonal* eikonal = new Eikonal();
  eikonal->setup(cfg);
  eikonal->solve();
  delete eikonal;

  msg("Done.");

  return 0;
}