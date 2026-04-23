#ifndef UTIL_H
#define UTIL_H

#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "input_file.h"
#include "interp.h"
#include "log.hpp"
#include "timer.hpp"
#include "command_line_args.h"
#include "time_parameters.hpp"

// -----------------------------------------------------------------------------
// Useful macros
// -----------------------------------------------------------------------------

#define error(msg) print_error(__FILE__, __func__,msg);

// -----------------------------------------------------------------------------
// General useful functions 
// -----------------------------------------------------------------------------

//! Prints a simple message to the stdout
void msg(std::string m);

//! Prints the title of the problem in boldface
void title(std::string t);

//! Prints an error message with more information and throws an error
void print_error(std::string file, std::string task, std::string reason);

//! Checks if filename exists
bool file_exists(const char * filename);

//! Checks if filename exists
bool file_exists(const std::string & filename);

//! Extracts file extension
std::string file_extension(const std::string & filename);

//! Prints a line to the stdout
std::string line();

//! Prints a nice dynamic progress bar
inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50);

//! Splits a string using a separator string 
void string_split(std::string str, 
		  std::string separator,
		  std::vector<std::string>* results);

//! Removes white spaces from string
void trim_white_spaces(std::string & str);

//! Parses a string to a arma::vec3
void parse_to_vec3(const std::string & str, arma::vec3 & v);

//! Parses a string to an std::vector of T
//void parse_to_vector(const std::string & str, std::vector<double> & v);
template<typename T>
void parse_to_vector(const std::string & str, std::vector<T> & v)
{
  std::string input(str);
  std::istringstream ss(input);
  std::string token;

  int i=0;
  T val;
  while(ss >> val)
  {
    v.push_back(val);
    if(ss.peek() == ',')
      ss.ignore();
    i++;
  }
}

// -----------------------------------------------------------------------------
// Useful functions for XML processing
// -----------------------------------------------------------------------------

void read_xml_section(const std::string & xmlname,
                      const std::string & section,
                      const std::string & subsection,
                      std::map<int,double> & pmap);

// -----------------------------------------------------------------------------
// Useful functions for using Armadillo
// -----------------------------------------------------------------------------

//! Converts a matrix to a vector in Voigt notation
void voigtvec(const int nd, const arma::mat & m, arma::mat & v);

//! Contraction of two second order tensors
double contract(const arma::mat & A, const arma::mat & B);

//! Computes the trace of a second order tensor
double trace(const arma::mat & A);

#endif
