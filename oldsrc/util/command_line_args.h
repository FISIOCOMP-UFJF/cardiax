#ifndef COMMAND_LINE_ARGS_H_
#define COMMAND_LINE_ARGS_H_

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

class CommandLineArgs
{
  static int argc;
  static const char ** argv;

  static char invalid[20];

public:

  /** Initialize command line arguments parser */
  static void init(int argc_, const char ** argv_)
  {
    strcpy(invalid, "____invalid____");
    argc = argc_;
    argv = argv_;
  }

  /** Read argument without value (switch) */
  static bool hasSwitch(const char * arg)
  {
    const char* thisArg;
    int tmp_argc = argc;
    const char ** tmp_argv = argv;

    while((thisArg = *++tmp_argv) && --tmp_argc){ // check all arguments
      if (!strcmp(thisArg, arg) && (tmp_argc > 0))
        return true; // return if match and ++tmp_argv exist
    }

    return false;
  }

  /** Read argument with double value */
  static double read(const char * arg, const double & def)
  {
    const char* pc = read(arg, invalid);

    if (strcmp(pc,invalid) == 0)
      return def;
    else
      return atof(pc);
  }

  /** Read argument with int value */
  static int read(const char * arg, const int & def)
  {
    const char* pc = read(arg, invalid);

    if (strcmp(pc, invalid) == 0)
      return def;
    else
      return atoi(pc);
  }

  /** Read argument with string value */
  static string read(const char * arg, const string & def)
  {
    const string pc = read(arg, invalid);

    if (strcmp(pc.c_str(),invalid) == 0)
      return def;
    else
      return pc;
  }

  /** Read argument with string value */
  static bool read(const char * arg, const bool & def)
  {
    const string pc = read(arg, invalid);

    if (strcmp(pc.c_str(),invalid) == 0)
      return def;
    else
    {
      if(strcmp(pc.c_str(),"false") == 0)
          return false;
      else
        return true;
    }
  }

  /** Read argument with value */
  static const char* read(const char* arg, const char * def)
  {
    const char* thisArg;
    int tmp_argc = argc;
    const char** tmp_argv = argv;

    while((thisArg = *++tmp_argv) && (--tmp_argc)) { // check all arguments
      if (!strcmp(thisArg, arg) && (tmp_argc > 0))
        return *++tmp_argv; // return if match and ++tmp_argv exist
    }

    return def; // no match, return default value
  }

};

#endif
