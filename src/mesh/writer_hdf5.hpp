#ifndef XDMFIO_HPP
#define	XDMFIO_HPP

#include <vector>
#include <string>
#include <fstream>
#include "hdf5.h"
#include "writer.hpp"

//
// Field description used to drive BOTH the HDF5 dataset creation and the
// XDMF attribute emission from a single registry. Adding a field is one
// add_field() call; it is then created, written and described automatically.
//
enum class FieldLoc  { Point, Cell };      // nodal (np rows) or cell (ne rows)
enum class FieldKind { Static, Transient };// one solution, or [nsteps, n, (ncomp)]

struct FieldSpec
{
  std::string name;        // dataset name under /point_data or /cell_data
  FieldLoc    loc;         // Point -> np, Cell -> ne
  FieldKind   kind;        // Static -> no time dim; Transient -> leading nsteps dim
  int         ncomp = 1;   // 1 = scalar, 3 = vector (e.g. displacements)
};

class WriterHDF5 : public Writer
{
public:

  //! Constructor
  WriterHDF5(Mesh * m);

  //! Destructor
  ~WriterHDF5();

  // --------------------------------------------------------------------------
  //  Field registry (NEW API)
  // --------------------------------------------------------------------------

  //! Register a field BEFORE open(). Drives dataset creation + XDMF emission.
  void add_field(const std::string & name, FieldLoc loc,
                 FieldKind kind, int ncomp = 1);

  //! Clear the registry (start a solver-specific field set from scratch).
  void clear_fields();

  //! Legacy default field set (vm, active_stress, ..., fibrosis, aha_marker,
  //! displacements). Called automatically by open() when no field was
  //! registered, so existing solvers keep producing the same file.
  //! TEMPORARY: remove the open()-time fallback once each app registers its
  //! own fields explicitly.
  void register_default_fields();

  // --------------------------------------------------------------------------
  //  Open / close / write
  // --------------------------------------------------------------------------

  //! Open file for writing HDF5 + XDMF data
  void open(const std::string & file, int nsteps, double step,
            bool bido = false, bool is_restart = false);

  //! Close HDF5 and XDMF files
  void close();

  //! Write one transient step of a node-centred scalar field
  void write_point_field_step(int step, const double * data, std::string fieldname);

  //! Write one transient step of a cell-centred scalar field
  void write_cell_field_step(int step, const double * data, std::string fieldname);

  //! Write one transient step of the displacement vector field [1, np, 3]
  void write_displ_step(int step, const double * data);

  //! Convenience wrappers (delegate to write_point_field_step)
  void write_vm_step(int step, const double * data);
  void write_ve_step(int step, const double * data);

  //! Write a STATIC field once (no time dimension). ncomp from the registry.
  void write_static_field(const double * data, std::string fieldname);

  // --------------------------------------------------------------------------
  //  Legacy field-adders kept for source compatibility (now thin wrappers)
  // --------------------------------------------------------------------------
  void add_ve();                                   //!< registers "ve" (Point, Transient)
  void add_fibers();                               //!< no-op (kept for compatibility)
  void add_scalar_field(std::string & field_name); //!< registers field_name (Point, Transient)

  //! Write data structure in a XDMF file (registry-driven)
  void write_xdmf(const std::string & file, int nsteps, double step, bool bido = false);

  //! Write static Eikonal LAT result to HDF5 and XDMF
  void write_eikonal_lat(const std::string & file, const double * lat_data);

  // --------------------------------------------------------------------------
  //  Checkpoints (unchanged; separate concern)
  // --------------------------------------------------------------------------
  void write_checkpoint(int step, double current_time, const double *vm,
                        const double *state_vars, int num_state_vars);

  void read_checkpoint_metadata(const std::string &filename, int &step,
                                double &time, int &num_nodes, int &num_vars);

  void read_checkpoint_data(const std::string &filename, double *vm, double *state_vars);

  void write_mech_checkpoint(int step, double current_time, int load_increment,
                             double load_factor, const double *x_current,
                             const double *fext0, int num_dofs);

  void read_mech_checkpoint_metadata(const std::string &filename, int &step,
                                     double &time, int &load_increment,
                                     double &load_factor, int &num_dofs);

  void read_mech_checkpoint_data(const std::string &filename,
                                 double *x_current, double *fext0);

  void write_coupled_checkpoint(int step, double current_time,
                                const double *vm, const double *state_vars,
                                int num_state_vars, int load_increment,
                                double load_factor, const double *x_current,
                                const double *fext0, int num_dofs);

private:

  //! XMF file stream
  std::ofstream xmf;

  //! Name of the HDF5 dataset file
  std::string h5name;

  //! Remembered time extent/spacing (set in open())
  std::size_t nsteps_ = 0;
  double      dt_     = 0.0;

  //! The field registry that drives HDF5 creation + XDMF emission
  std::vector<FieldSpec> fields_;

  //! Look up a registered field by name (nullptr if absent)
  const FieldSpec * find_field(const std::string & name) const;

  //! Create (truncate) the HDF5 file, write geometry/topology/time, and
  //! create every registered field's dataset.
  void write_hdf5(const std::string & file, int nsteps, double step);

  //! Emit one <Attribute> block for field `fs` at time index `step`.
  void emit_attribute(std::ostream & os, const FieldSpec & fs,
                      int step, int np, int ne, int nsteps) const;

  //! XDMF topology-type string from element node count + dimension.
  std::string xdmf_topology_type() const;
};

#endif	/* XDMFIO_HPP */


// #ifndef XDMFIO_HPP
// #define	XDMFIO_HPP

// #include <vector>
// #include <fstream>
// #include "hdf5.h"
// #include "writer.hpp"

// class WriterHDF5 : public Writer
// {
// public:

//   //! Constructor
//   WriterHDF5(Mesh * m);
    
//   //! Destructor
//   ~WriterHDF5();

//   //! Add another scalar field for writing data
//   void add_ve();

//   //! Add fiber vector fields for writing data
//   void add_fibers();

//   //! Add scalar field for writing data
//   void add_scalar_field(std::string & field_name);

//   //! Close HDF5 and XDMF files
//   void close();

//   //! Open file for writing HDF5 data
//   void open(const std::string & file, int nsteps, double step, bool bido=false,bool is_restart=false);

//   //! Write one time step of vm scalar field
//   void write_vm_step(int step, const double *data);

//   void write_cell_field_step(int step, const double *data, string fieldname);

//   //! Write one time step of a node-centred scalar field (vm, active_stress)
//   void write_point_field_step(int step, const double *data, string fieldname);

//   //! Write one time step of displacement field
//   void write_displ_step(int step, const double *data);

//   //! Write one time step of ve scalar field
//   void write_ve_step(int step, const double *data);
    
//   //! Write data structure in a XDMF file
//   void write_xdmf(const std::string & file, int nsteps, double step, bool bido=false);

//   //! Write static Eikonal LAT result to HDF5 and XDMF
//   void write_eikonal_lat(const std::string & file, const double *lat_data);
  
//   //! Write checkpoint state file
//   void write_checkpoint(int step, double current_time, const double *vm, const double *state_vars, int num_state_vars);
  
//   void read_checkpoint_metadata(const std::string &filename, int &step, double &time, int &num_nodes, int &num_vars);

//   void read_checkpoint_data(const std::string &filename, double *vm, double *state_vars);

//   void write_mech_checkpoint(int step, double current_time, int load_increment, double load_factor,
//                              const double *x_current, const double *fext0, int num_dofs);

//   void read_mech_checkpoint_metadata(const std::string &filename, 
//                                      int &step, double &time, int &load_increment, 
//                                      double &load_factor, int &num_dofs);

//   void read_mech_checkpoint_data(const std::string &filename, 
//                                  double *x_current, double *fext0);

//   void write_coupled_checkpoint(int step, double current_time, 
//                                 const double *vm, const double *state_vars, int num_state_vars,
//                                 int load_increment, double load_factor, 
//                                 const double *x_current, const double *fext0, int num_dofs);

// private:
    
//   //! XMF File
//   std::ofstream xmf;

//   //! Name of the HDF5 dataset file
//   std::string h5name;
    
//   //! Write the raw data to HDF5 file
//   void write_hdf5(const std::string & file, int nsteps, double step);
    
// };

// #endif	/* XDMFIO_HPP */




