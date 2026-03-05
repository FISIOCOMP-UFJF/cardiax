#ifndef GMSHIO_HPP
#define GMSHIO_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "util/util.hpp"
#include "mesh.hpp"

class GmshIO 
{
private:

  Mesh & mesh;

  void read_mesh(std::istream & in);
  
public:

  GmshIO(Mesh & m) : mesh(m) {}

  void read(const std::string & name);

};

#endif
