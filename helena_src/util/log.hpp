#ifndef LOG_HPP_
#define LOG_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <stdarg.h>

/**
 * A simple class for Logging events into a text file.
 */
class Log
{
public:

  //! Default constructor
  Log() {};

  //! Constructor
  Log(const std::string & filename);

  //! Destructor
  ~Log();

  //! Create file for logging
  void create(const std::string & filename);

  //! Operator << for string
  Log & operator<<(std::string str);

  //! Operator << for numerical values
  Log & operator<<(double value);

  //! Operator << for ostream stuff
  Log & operator<<(std::ostream& (*p) (std::ostream &));

  //! Member function that writes to the file stream.
  void write(char* logline);

  //! Variadic version of write function. 
  //! Use in order to avoid using a temporary buffer and sprintf.
  void write(const char* logline, ...);

private:

  //! File stream to write the logs.
  std::ofstream _stream;

};


#endif


