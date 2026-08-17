#ifndef XDMFIO_HPP
#define	XDMFIO_HPP

#include <vector>
#include <fstream>
#include "hdf5.h"
#include "writer.hpp"

class WriterHDF5 : public Writer
{
public:

  //! Constructor
  WriterHDF5(Mesh * m);
    
  //! Destructor
  ~WriterHDF5();

  //! Add another scalar field for writing data
  void add_ve();

  //! Add fiber vector fields for writing data
  void add_fibers();

  //! Add scalar field for writing data
  void add_scalar_field(std::string & field_name);

  //! Close HDF5 and XDMF files
  void close();

  //! Open file for writing HDF5 data
  void open(const std::string & file, int nsteps, double step, bool bido=false,bool is_restart=false);

  //! Write one time step of vm scalar field
  void write_vm_step(int step, const double *data);

  void write_cell_field_step(int step, const double *data, string fieldname);

  //! Write one time step of displacement field
  void write_displ_step(int step, const double *data);

  //! Write one time step of ve scalar field
  void write_ve_step(int step, const double *data);
    
  //! Write data structure in a XDMF file
  void write_xdmf(const std::string & file, int nsteps, double step, bool bido=false);

  //! Write checkpoint state file
  // LUCAS:
  void write_checkpoint(int step, double current_time, const double *vm, const double *state_vars, int num_state_vars);
  
  void read_checkpoint_metadata(const std::string &filename, int &step, double &time, int &num_nodes, int &num_vars);

  void read_checkpoint_data(const std::string &filename, double *vm, double *state_vars);

  void write_mech_checkpoint(int step, double current_time, int load_increment, double load_factor,
                             const double *x_current, const double *fext0, int num_dofs);

  void read_mech_checkpoint_metadata(const std::string &filename, 
                                     int &step, double &time, int &load_increment, 
                                     double &load_factor, int &num_dofs);

  void read_mech_checkpoint_data(const std::string &filename, 
                                 double *x_current, double *fext0);

private:
    
  //! XMF File
  std::ofstream xmf;

  //! Name of the HDF5 dataset file
  std::string h5name;
    
  //! Write the raw data to HDF5 file
  void write_hdf5(const std::string & file, int nsteps, double step); 
    
};

#endif	/* XDMFIO_HPP */




