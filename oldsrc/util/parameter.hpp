#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>

class Parameter
{
public:

	//! Create parameter for given key
	explicit Parameter(std::string key);

	//! Destructor
	virtual ~Parameter();

	//! Return parameter key
	std::string key() const;

	//! Return parameter description
	std::string info() const;

	//! Return true if parameter is set, return false otherwise
	bool is_set() const;

	//! Reset the parameter to empty, so that is_set() returns false.
	void reset();

	//! Assignment from int
	virtual const Parameter& operator= (int value);

	//! Assignment from double
	virtual const Parameter& operator= (double value);

	//! Assignment from string
	virtual const Parameter& operator= (std::string value);

	//! Assignment from string
	virtual const Parameter& operator= (const char* value);

	//! Assignment from bool
	virtual const Parameter& operator= (bool value);

  //! Cast parameter to int
  virtual operator int() const;

  //! Cast parameter to std::size_t
  virtual operator std::size_t() const;

  //! Cast parameter to double
  virtual operator double() const;

  //! Cast parameter to string
  virtual operator std::string() const;

  //! Cast parameter to bool
  virtual operator bool() const;

  //! Return the name of the parameter type (int, double, etc...)
  virtual std::string type_str() const = 0;

  //! Return the value as a string
  virtual std::string value_str() const = 0;

	// Check that key name is allowed
	static void check_key(std::string key);

protected:

	// Whether or not parameter has been set
	bool p_is_set;

	// Parameter key
	std::string pkey;

	// Parameter description
	std::string pdescription;

};

// -----------------------------------------------------------------------------

class IntParameter : public Parameter
{
public:

	//! Create unset int-valued
	explicit IntParameter(std::string key);

	//! Create int-valued parameter
	IntParameter(std::string key, int val);

	//! Destructor
	~IntParameter(){};

	//! Assignment
	const IntParameter& operator= (int val);

	//! Cast parameter to int
	operator int() const;

	//! Return short string description
	std::string str() const;

  //! Return int
  std::string type_str() const;

  //! Return value
  std::string value_str() const;

private:

	//! Parameter value
	int value;

};

// -----------------------------------------------------------------------------

class DoubleParameter : public Parameter
{
public:

	//! Create unset double-valued parameter
	explicit DoubleParameter(std::string key);

	//! Create double-valued parameter
	DoubleParameter(std::string key, double val);

	//! Destructor
	~DoubleParameter(){};

	//! Assignment
	const DoubleParameter& operator= (double val);

	//! Cast parameter to double
	operator double() const;

	//! Return short string description
	std::string str() const;

  //! Return double
  std::string type_str() const;

  //! Return value
  std::string value_str() const;

private:

	//! Parameter value
	double value;

};

// -----------------------------------------------------------------------------

class StringParameter : public Parameter
{
public:

	//! Create unset string-valued parameter
	explicit StringParameter(std::string key);

	//! Create string-valued parameter
	StringParameter(std::string key, std::string val);

	//! Destructor
	~StringParameter(){};

	//! Assignment
	const StringParameter& operator= (std::string val);

	//! Assignment
	const StringParameter& operator= (const char* val);

	//! Cast parameter to string
	operator std::string() const;

	//! Return short string description
	std::string str() const;

  //! Return string
  std::string type_str() const;

  //! Return value
  std::string value_str() const;

private:

	//! Parameter value
	std::string value;

};

// -----------------------------------------------------------------------------

class BoolParameter : public Parameter
{
public:

	//! Create unset bool-valued parameter
	explicit BoolParameter(std::string key);

	//! Create bool-valued parameter
	BoolParameter(std::string key, bool value);

	//! Destructor
	~BoolParameter(){};

	//! Assignment
	const BoolParameter& operator= (bool value);

	//! Cast parameter to bool
	operator bool() const;

	//! Return short string description
	std::string str() const;

  //! Return bool
  std::string type_str() const;

  //! Return value
  std::string value_str() const;

private:

	//! Parameter value
	bool value;

};



#endif
