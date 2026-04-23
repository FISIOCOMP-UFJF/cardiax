#ifndef WRITER_HPP
#define	WRITER_HPP

#include <string>
#include "mesh.hpp"

/**
 * Abstract base class for writing data
 * */

class Writer
{
public:

  //! Constructor
  Writer(Mesh * m);
    
  //! Virtual Destructor
  virtual ~Writer();
    
  //! Open the file to write data
  virtual void open(const std::string & file, int nsteps, double step, bool bido=false) = 0;

  //! Close the file and everything
  virtual void close() = 0;
    
  //! Write transmembrane potential
  virtual void write_vm_step(int step, const double * data) = 0;
    
  //! Write displacements
  virtual void write_displ_step(int step, const double * data) = 0;
    
protected:

  //! Reference to the computational mesh.
  Mesh * mesh;

};

#endif	/* WRITER_HPP */

