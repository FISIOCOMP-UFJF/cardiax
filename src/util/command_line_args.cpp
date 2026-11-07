#include "command_line_args.h"

int CommandLineArgs::argc;
const char** CommandLineArgs::argv;
char CommandLineArgs::invalid[20];

std::map<std::string, CommandLineArgs::Entry> CommandLineArgs::registry;

void override_str(std::string & value, const char * flag)
{
  static const std::string SENT = "\x01__unset__\x01";
  std::string got = CommandLineArgs::read(flag, SENT);
  if (got != SENT) value = got;
}

void override_int(int & value, const char * flag)
{
  const int SENT = std::numeric_limits<int>::lowest();
  int got = CommandLineArgs::read(flag, SENT);
  if (got != SENT) value = got;
}

void override_double(double & value, const char * flag)
{
  const double SENT = std::numeric_limits<double>::lowest();
  double got = CommandLineArgs::read(flag, SENT);
  if (got != SENT) value = got;
}
