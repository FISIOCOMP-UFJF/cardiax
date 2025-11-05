#include <sstream>
#include "util.hpp"
#include "parameter.hpp"

// -----------------------------------------------------------------------------
// base class Parameter
// -----------------------------------------------------------------------------

Parameter::Parameter(std::string key)
  : p_is_set(false), pkey(key), pdescription("param description")
{
  // check that key name is allowed
  check_key(key);
}

Parameter::~Parameter()
{
  // do nothing ...
}

std::string Parameter::key() const
{
  return pkey;
}

std::string Parameter::info() const
{
  return pdescription;
}

bool Parameter::is_set() const
{
  return p_is_set;
}

void Parameter::reset()
{
  p_is_set = false;
}

const Parameter& Parameter::operator= (int val)
{  
  return *this;
}

const Parameter& Parameter::operator= (double val)
{
  return *this;
}

const Parameter& Parameter::operator= (std::string val)
{
  return *this;
}

const Parameter& Parameter::operator= (const char* val)
{
  return *this;
}

const Parameter& Parameter::operator= (bool val)
{
  return *this;
}

Parameter::operator int() const
{
  print_error("parameter.cpp", "int()", "base class cannot cast to int");
  return 0;
}

Parameter::operator std::size_t() const
{
  print_error("parameter.cpp", "size_t()", "base class cannot cast to size_t");
  return 0;
}

Parameter::operator double() const
{
  print_error("parameter.cpp", "double()", "base class cannot cast to double");
  return 0;
}

Parameter::operator std::string() const
{
  print_error("parameter.cpp", "string()", "base class cannot cast to string");
  return 0;
}

Parameter::operator bool() const
{
  print_error("parameter.cpp", "bool()", "base class cannot cast to bool");
  return 0;
}

void Parameter::check_key(std::string key)
{
  // Space and punctuation not allowed in key names
  for (std::size_t i = 0; i < key.size(); i++)
  {
    if (key[i] == ' ' || key[i] == '.')
      print_error("parameter.hpp", "check_key", "invalid name for key");
  }
}

// -----------------------------------------------------------------------------
// class IntParameter
// -----------------------------------------------------------------------------

IntParameter::IntParameter(std::string key)
  : Parameter(key)
{
  // Do nothing
}

IntParameter::IntParameter(std::string key, int val)
  : Parameter(key), value(val)
{
  p_is_set = true;
}

const IntParameter& IntParameter::operator= (int val)
{
  value = val;
  p_is_set = true;
  return *this;
}

IntParameter::operator int() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "int()", "parameter has not been set");

  return value;
}

std::string IntParameter::str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "str()", "parameter has not been set");

  std::stringstream s;
  s << "<parameter type=int, name=\"" << key()
    << "\", value=" << value << ">";
  return s.str();
}

std::string IntParameter::value_str() const
{
  if (!p_is_set)
      error("Parameter has not been set");

  std::stringstream s;
  s << value;
  return s.str();
};

std::string IntParameter::type_str() const
{
  return "int";
}

// -----------------------------------------------------------------------------
// class DoubleParameter
// -----------------------------------------------------------------------------

DoubleParameter::DoubleParameter(std::string key)
  : Parameter(key)
{
  // Do nothing
}

DoubleParameter::DoubleParameter(std::string key, double val)
  : Parameter(key), value(val)
{
  p_is_set = true;
}

const DoubleParameter& DoubleParameter::operator= (double val)
{
  value = val;
  p_is_set = true;
  return *this;
}

DoubleParameter::operator double() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "double()", "parameter has not been set");

  return value;
}

std::string DoubleParameter::str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "str()", "parameter has not been set");

  std::stringstream s;
  s << "<parameter type=double, name=\"" << key()
    << "\", value=" << value << ">";
  return s.str();
}

std::string DoubleParameter::value_str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "value_str()", "parameter has not been set");

  std::stringstream s;
  s << value;
  return s.str();
};

std::string DoubleParameter::type_str() const
{
  return "double";
}

// -----------------------------------------------------------------------------
// class StringParameter
// -----------------------------------------------------------------------------

StringParameter::StringParameter(std::string key) : Parameter(key)
{
  // Do nothing
}

StringParameter::StringParameter(std::string key, std::string val)
  : Parameter(key), value(val)
{
	p_is_set = true;
}

const StringParameter& StringParameter::operator= (std::string val)
{
  value = val;
  p_is_set = true;
  return *this;
}

const StringParameter& StringParameter::operator= (const char* val)
{
  std::string s(val);
  value = s;
  p_is_set = true;
  return *this;
}

StringParameter::operator std::string() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "string()", "parameter has not been set");

  return value;
}

std::string StringParameter::str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "str()", "parameter has not been set");

  std::stringstream s;
  s << "<parameter type=string, name=\"" << key()
    << "\", value=\"" << value << "\">";
  return s.str();
}

std::string StringParameter::value_str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "value_str()", "parameter has not been set");

  return value;
};

std::string StringParameter::type_str() const
{
  return "string";
}

// -----------------------------------------------------------------------------
// class BoolParameter
// -----------------------------------------------------------------------------

BoolParameter::BoolParameter(std::string key) : Parameter(key)
{
  // do nothing
}

BoolParameter::BoolParameter(std::string key, bool val)
  : Parameter(key), value(val)
{
  p_is_set = true;
}

const BoolParameter& BoolParameter::operator= (bool val)
{
  value = val;
  p_is_set = true;
  return *this;
}

BoolParameter::operator bool() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "bool()", "parameter has not been set");

  return value;
}

std::string BoolParameter::str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "str()", "parameter has not been set");

  std::stringstream s;
  s << "<parameter type=bool, name=\"" << key()
    << "\", value=" << value << ">";
  return s.str();
}

std::string BoolParameter::value_str() const
{
  if (!p_is_set)
    print_error("parameter.cpp", "value_str()", "parameter has not been set");

  if (value)
    return "true";
  else
    return "false";
}

std::string BoolParameter::type_str() const
{
  return "bool";
}