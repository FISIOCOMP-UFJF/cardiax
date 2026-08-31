#include "command_line_args.h"

int CommandLineArgs::argc;
const char** CommandLineArgs::argv;
char CommandLineArgs::invalid[20];

std::map<std::string, CommandLineArgs::Entry> CommandLineArgs::registry;
