#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <map>
#include <vector>
#include "parameter.hpp"

/**
 This class holds parameters, where each parameter is composed of
 a key (name) and a value of some type (int, double, bool, string).
 Parameter sets can be nested. Example of usage:

   Parameters p("my_parameters");
   p.add("tolerance",  1e-15);
   p.add("maxits",       30);
   p.add("rebuild",   false);

   p["maxits"] = 50;
   int its = p["maxits"];


   Parameters q("nested_parameters");
   p.add(q);

 Nested parameters are retrieved by

   p("nested_parameters")["..."]

 TODO: need to implement something to parse from cmdline
   p.parse(argc, argv);

*/

class Parameters
{
public:

	//! Create empty parameter set
	Parameters(std::string key = "parameters");

	//! Destructor
	virtual ~Parameters();

	//! Copy constructor
	Parameters(const Parameters& parameters);

	//! Return name for parameter set
	std::string name() const;

	//! Rename parameter set
	void rename(std::string name);

	//! Clear parameter set
	void clear();

	/// Add int-valued parameter
	void add(std::string key, int value);

	//! Add double-valued parameter
	void add(std::string key, double value);

	//! Add string-valued parameter
	void add(std::string key, std::string value);

	//! Add string-valued parameter
	void add(std::string key, const char* value);

	//! Add bool-valued parameter
	void add(std::string key, bool value);

	//! Add nested parameter set
	void add(const Parameters& parameters);

  //! Parse parameters from command-line
  //virtual void parse(int argc, char* argv[]);

  //! Remove parameter or parameter set with given key
	void remove(std::string key);

	//! Update parameters with another set of parameters
	//void update(const Parameters& parameters);

	//! Return parameter for given key
	Parameter& operator[] (std::string key);

	//! Return parameter for given key (const version)
	const Parameter& operator[] (std::string key) const;

	//! Return nested parameter set for given key
	Parameters& operator() (std::string key);

	//! Return nested parameter set for given key (const)
	const Parameters& operator() (std::string key) const;

	//! Assignment operator
	const Parameters& operator= (const Parameters& parameters);

	//! Check if parameter set has key (parameter or nested parameter set)
	bool has_key(std::string key) const;

	//! Check if parameter set has given parameter
	bool has_parameter(std::string key) const;

	//! Check if parameter set has given nested parameter set
	bool has_parameter_set(std::string key) const;

	//! Return a vector of parameter keys
	void get_parameter_keys(std::vector<std::string>& keys) const;

	//! Return a vector of parameter set keys
	void get_parameter_set_keys(std::vector<std::string>& keys) const;

	//! Return informal string representation (pretty-print)
	std::string str() const;

	// Return pointer to parameter for given key and 0 if not found
	Parameter* find_parameter(std::string key) const;

	// Return pointer to parameter set for given key and 0 if not found
	Parameters* find_parameter_set(std::string key) const;

  // Read parameters from a xml file
  void parse_xml(std::string filename, std::string nodename);

private:

	// Parameter set key
	std::string _key;

	// Map from key to parameter
	std::map<std::string, Parameter*> _parameters;

	// Map from key to parameter sets
	std::map<std::string, Parameters*> _parameter_sets;

};


#endif
