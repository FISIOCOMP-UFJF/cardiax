#ifndef INPUT_FILE_H_INCLUDED
#define INPUT_FILE_H_INCLUDED

#include <string>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <armadillo>

#include "util.hpp"

/** Nodal data information pair for use in maps and multimaps
       int : spatial direction 0,1,2 (x,y,z)
       double: displacement, force, or something else
*/
typedef std::pair<int,double> NodalData;

/**
 * Input file class
 * */

class InputFile
{
public:

  //! Constructor - with char
  InputFile(const char * fname);

  //! Constructor - with string
  InputFile(const std::string & fname);

  //! Desctructor
  ~InputFile();

  //! Close file and finish reading
  void close();

  //! Read parameter (key = int value)
  void read(const char * key, int & value);

  //! Read parameter (key = double value)
  void read(const char * key, double & value);

  //! Read parameter (key = double value)
  void read(const char * key, std::string & value);

  //! Read all the data from an array (e.g.: material properties)
  void read_array(const char * name, std::vector<double> & array);

  /*! Given the name of a section of the type

      section_name
         inta0 intb0
         inta1 intb1
         ...
         intaN intbN
      end_section_name

      reads its content and stores in a std::map<int,int>
  */
  void read_section(const char * name, std::multimap<int,int> & parmap);

  /*! Given the name of a section of the type

      section_name
         int0 double0
         int1 double1
         ...
         intn doublen
      end_section_name

      reads its content and stores in a std::map<int,double>
  */
  void read_section(const char * name, std::map<int,double> & parmap);

  /*! Given the name of a section of the type

      section_name
         int0 fx0 fy0 fz0
         int1 fx1 fy1 fz1
         ...
         intn fxn fyn fzn
      end_section_name

      this member functions reads its content and stores
      it in a std::map<int,arma::vec3>. Useful for traction (Neumann)
      boundary conditions in elasticity problems.
  */
  void read_section(const char * name, std::map<int,arma::vec3> & parmap);

  /*! Given the name of a section of the type

      section_name
         int0 dir double0
         int1 dir double1
         ...
         intn dir doublen
      end_section_name

      where dir is:
         x = 0
         y = 1
         z = 2
      this member functions reads its content and stores
      it in a std::multimap<int,NodalData>
  */
  void read_section(const char * name, std::multimap<int,NodalData> & parmap);

private:

  //! Name of the input file to read options
  std::string filename;

  //! File reader
  std::ifstream is;

  //! Separator of key-value (=)*/
  std::string separator;

  //! Find the beginning of a section
  int find_section(const char * name);

  //! Given a string with the contents of a line, read key and value
  void read_key_value(const std::string & line, std::string & k, std::string & v);

};

#endif
