#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include <vector>
#include <fstream>
#include "mesh.hpp"
#include "writer.hpp"

/** A simple class to read and write VTK data for visualization
    using any software that supports VTK, like Paraview.
*/

class WriterVTK : public Writer
{
public:

  //! Constructor 
  WriterVTK(Mesh * m);

  //! Destructor
  ~WriterVTK();

  void open(const std::string & file, int nsteps, double step, bool bido=false);
  void close();
  void write_vm_step(int step, const double *data);
  void write_displ_step(int step, const double *data);

private:
  
  //! XML File
  std::ofstream xmlfile;
    
  //! Name of the XML file
  std::string filename;

  //! Write the XML structure
  void write_xml(const std::string & file, int nsteps, double step, bool bido);

  //! Write data to the appended section
  void write_bin(const std::string & file);

};


#endif
