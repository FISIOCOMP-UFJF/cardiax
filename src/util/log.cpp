#include "log.hpp"

Log::Log(const std::string & filename)
{
  _stream.open(filename.c_str());
}

Log::~Log()
{
  _stream.close();
}

void Log::create(const std::string & filename)
{
  _stream.open(filename.c_str());
}

Log & Log::operator<<(std::string str)
{
  if( _stream.is_open() )
    _stream << str;
  
  return *this;
}

Log & Log::operator<<(double value)
{
  if( _stream.is_open() )
    _stream << value;
  
  return *this;
}

Log & Log::operator<<(std::ostream& (*p) (std::ostream &))
{
  if( _stream.is_open() )
    _stream << p;
  
  return *this;
}

void Log::write(char * logline)
{
  _stream << logline << std::endl;
}

void Log::write(const char * logline, ...)
{
  va_list argList;
  char cbuffer[1024];
  va_start(argList, logline);
  vsnprintf(cbuffer, 1024, logline, argList);
  va_end(argList);
  _stream << cbuffer;
}
