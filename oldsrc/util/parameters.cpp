#include <sstream>
#include <stdio.h>
#include <cassert>
#include "util.hpp"
#include "parameter.hpp"
#include "parameters.hpp"
#include "fem/fe.hpp"
#include "util/pugixml.hpp"

// Typedef of iterators for convenience
typedef std::map<std::string, Parameter*>::iterator parameter_iterator;
typedef std::map<std::string, Parameter*>::const_iterator const_parameter_iterator;
typedef std::map<std::string, Parameters*>::iterator parameter_set_iterator;
typedef std::map<std::string, Parameters*>::const_iterator const_parameter_set_iterator;

Parameters::Parameters(std::string key) : _key(key)
{
  Parameter::check_key(key);
}

Parameters::~Parameters()
{
  clear();
}

Parameters::Parameters(const Parameters& parameters)
{
  *this = parameters;
}

std::string Parameters::name() const
{
  return _key;
}

void Parameters::rename(std::string name)
{
  Parameter::check_key(name);
  _key = name;
}

void Parameters::clear()
{
  // clear parameters
  for (parameter_iterator it = _parameters.begin(); it != _parameters.end(); ++it)
    delete it->second;
  _parameters.clear();

  // clear (nested) parameter sets
  for (parameter_set_iterator it = _parameter_sets.begin(); it != _parameter_sets.end(); ++it)
    delete it->second;
  _parameter_sets.clear();

  // reset key
  _key = "";
}

void Parameters::add(std::string key, int value)
{
  if (find_parameter(key))
    error("parameter already defined");

  _parameters[key] = new IntParameter(key, value);
}

void Parameters::add(std::string key, double value)
{
  if (find_parameter(key))
    error("parameter already defined");

  _parameters[key] = new DoubleParameter(key, value);
}

void Parameters::add(std::string key, std::string value)
{
  if (find_parameter(key))
    error("parameter already defined");

  _parameters[key] = new StringParameter(key, value);
}

void Parameters::add(std::string key, const char* value)
{
  // This version is needed to avoid having const char* picked up by
  // the add function for bool parameters.

  if (find_parameter(key))
    error("parameter already defined");

  _parameters[key] = new StringParameter(key, value);
}

void Parameters::add(std::string key, bool value)
{
  if (find_parameter(key))
    error("parameter already defined");

  _parameters[key] = new BoolParameter(key, value);
}

void Parameters::add(const Parameters& parameters)
{
  if (find_parameter_set(parameters.name()))
    error("parameter already defined");

  Parameters* p = new Parameters("");
  *p = parameters;
  _parameter_sets[parameters.name()] = p;
}

void Parameters::remove(std::string key)
{
  if (!find_parameter(key) && !find_parameter_set(key))
    error("no parameter defined");

  // delete objects (safe to delete both even if only one is nonzero)
  delete find_parameter(key);
  delete find_parameter_set(key);

  // remove from maps (safe to remove both)
  std::size_t num_removed = 0;
  num_removed += _parameters.erase(key);
  num_removed += _parameter_sets.erase(key);
  assert(num_removed == 1);
}


Parameter& Parameters::operator[] (std::string key)
{
  Parameter* p = find_parameter(key);
  if (!p)
    error("parameter is not defined");

  return *p;
}

const Parameter& Parameters::operator[] (std::string key) const
{
  Parameter* p = find_parameter(key);
  if (!p)
    error("parameter is not defined");

  return *p;
}

Parameters& Parameters::operator() (std::string key)
{
  Parameters* p = find_parameter_set(key);
  if (!p)
    error("parameter is not defined");
	
  return *p;
}

const Parameters& Parameters::operator() (std::string key) const
{
  Parameters* p = find_parameter_set(key);
  if (!p)
    error("parameter is not defined");

  return *p;
}

const Parameters& Parameters::operator= (const Parameters& parameters)
{
  clear();

  // We rely on the default copy constructors for the Parameter subclasses
  // to do work properly since they don't use any dynamically allocated data

  // We need to copy: _key, _parameters and _parameter_sets
  // Copy key
  _key = parameters._key;

  // Copy parameters
  for (const_parameter_iterator it = parameters._parameters.begin();
       it != parameters._parameters.end(); ++it)
  {
    const Parameter& p = *it->second;
    Parameter* q = 0;

    if (p.type_str() == "int")
      q = new IntParameter(dynamic_cast<const IntParameter&>(p));
    else if (p.type_str() == "double")
      q = new DoubleParameter(dynamic_cast<const DoubleParameter&>(p));
    else if (p.type_str() == "bool")
      q = new BoolParameter(dynamic_cast<const BoolParameter&>(p));
    else if (p.type_str() == "string")
      q = new StringParameter(dynamic_cast<const StringParameter&>(p));
    else
      error("unknown parameter type");

    _parameters[p.key()] = q;
  }

  // Copy parameter sets
  for (const_parameter_set_iterator it = parameters._parameter_sets.begin();
       it != parameters._parameter_sets.end(); ++it)
  {
    const Parameters& p = *it->second;
    _parameter_sets[p.name()] = new Parameters(p);
  }

  return *this;
}

bool Parameters::has_key(std::string key) const
{
  return has_parameter(key) || has_parameter_set(key);
}

bool Parameters::has_parameter(std::string key) const
{
  return find_parameter(key) != 0;
}

bool Parameters::has_parameter_set(std::string key) const
{
  return find_parameter_set(key) != 0;
}

void Parameters::get_parameter_keys(std::vector<std::string>& keys) const
{
  keys.reserve(_parameters.size());
  for (const_parameter_iterator it = _parameters.begin(); it != _parameters.end(); ++it)
    keys.push_back(it->first);
}

void Parameters::get_parameter_set_keys(std::vector<std::string>& keys) const
{
  keys.reserve(_parameter_sets.size());
  for (const_parameter_set_iterator it = _parameter_sets.begin(); it != _parameter_sets.end(); ++it)
    keys.push_back(it->first);
 }

std::string Parameters::str() const
{
  std::stringstream s;

  s << "Parameter set \"" << name() << "\"\n"
    << " Number of parameters: " << _parameters.size() << "\n"
    << " Number of nested parameter sets: " << _parameter_sets.size() << "\n";

  if (_parameters.empty() && _parameter_sets.empty())
  {
    s << "Parameters: " << name() << " is empty\n";
    return s.str();
  }

  for (const_parameter_iterator it = _parameters.begin();
       it != _parameters.end(); ++it)
  {
    Parameter* p = it->second;
    s << " Parameter=" << p->key() << ", ";
    s << "type=" << p->type_str() << ", ";
    s << "value=" << (p->is_set() ? p->value_str() : "<unset>") << "\n";
  }

  for (const_parameter_set_iterator it = _parameter_sets.begin();
       it != _parameter_sets.end(); ++it)
  {
    s << "\n " << it->second->str();
  }

  return s.str();
}

Parameter* Parameters::find_parameter(std::string key) const
{
  const_parameter_iterator p = _parameters.find(key);
  if (p == _parameters.end())
    return 0;
  return p->second;
}

Parameters* Parameters::find_parameter_set(std::string key) const
{
  const_parameter_set_iterator p = _parameter_sets.find(key);
  if (p == _parameter_sets.end())
    return 0;
  return p->second;
}

void Parameters::parse_xml(std::string filename, std::string nodename)
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
    exit(1);
  }

  pugi::xml_node p = doc.child(nodename.c_str()).child("parameters");

  for (pugi::xml_node_iterator it = p.begin(); it != p.end(); ++it)
  {
    if(has_parameter(it->name()))
    {
      Parameter *p = find_parameter(it->name());
      if (p->type_str() == "double")
        *p = it->text().as_double();
      else if (p->type_str() == "int")
        *p = it->text().as_int();
      else if (p->type_str() == "bool")
        *p = it->text().as_bool();
      else if (p->type_str() == "string")
        *p = it->text().as_string();
    }
  }

}

/*
void Parameters::update(const Parameters& parameters)
{
  // Update the parameters
  for (const_parameter_iterator it = parameters._parameters.begin();
       it != parameters._parameters.end(); ++it)
  {
    // Get parameters
    const Parameter& other = *it->second;
    Parameter* self = find_parameter(other.key());

    // Skip parameters not in this parameter set (no new parameters added)
    if (!self)
    {
      //warning("Ignoring unknown parameter \"%s\" in parameter set \"%s\" when updating parameter set \"%s\".", other.key().c_str(), parameters.name().c_str(), name().c_str());
      continue;
    }

    // Skip unset parameters
    if (!other.is_set())
    {
      warning("Ignoring unset parameter \"%s\" in parameter set \"%s\" when updating parameter set \"%s\".",
              other.key().c_str(), parameters.name().c_str(), name().c_str());
      continue;
    }

    // Set value (will give an error if the type is wrong)
    if (other.type_str() == "int")
      *self = static_cast<int>(other);
    else if (other.type_str() == "double")
      *self = static_cast<double>(other);
    else if (other.type_str() == "bool")
      *self = static_cast<bool>(other);
    else if (other.type_str() == "string")
      *self = static_cast<std::string>(other);
    else
      error("Parameters.cpp",
			     "update parameter set",
		       "Parameter \"%s\" has unknown type: \"%s\"",
			      other.key().c_str(), other.type_str().c_str());
  }

  // Update nested parameter sets
  for (const_parameter_set_iterator it = parameters._parameter_sets.begin(); it != parameters._parameter_sets.end(); ++it)
    (*this)(it->first).update(*it->second);
}
*/


/*
void Parameters::read_vm(po::variables_map& vm, Parameters &parameters, std::string base_name)
{
  // Read values from po::variables_map
  for (parameter_iterator it = parameters._parameters.begin();
       it != parameters._parameters.end(); ++it)
  {
    Parameter& p = *it->second;
    std::string param_name(base_name + p.key());
    if (p.type_str() == "int")
    {
      const po::variable_value& v = vm[param_name];
      if (!v.empty())
        p = v.as<int>();
    }
    else if (p.type_str() == "bool")
    {
      const po::variable_value& v = vm[param_name];
      if (!v.empty())
        p = v.as<bool>();
    }
    else if (p.type_str() == "double")
    {
      const po::variable_value& v = vm[param_name];
      if (!v.empty())
        p = v.as<double>();
    }
    else if (p.type_str() == "string")
    {
      const po::variable_value& v = vm[param_name];
      if (!v.empty())
        p = v.as<std::string>();
    }
  }

  for (parameter_set_iterator it = parameters._parameter_sets.begin(); it != parameters._parameter_sets.end(); ++it)
  {
    read_vm(vm, *it->second, base_name + it->first + ".");
  }
}
*/