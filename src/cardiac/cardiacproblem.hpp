#ifndef CARDIAC_PROBLEM_HPP
#define CARDIAC_PROBLEM_HPP

#include "stimulus.hpp"
#include "linalg/linalg.h"
#include "odes/odes.h"
#include "mesh/writer_hdf5.hpp"
#include "util/parameters.hpp"

typedef std::vector<arma::mat33*> ArrayMat33;

/**
   Abstract base class for transient cardiac problems
   such as Bidomain or Monodomain models. 
*/

class CardiacProblem
{
public:

  //! Constructor
  CardiacProblem();

  //! Virtual destructor
  virtual ~CardiacProblem();

  //! Interface - perform one time step
  virtual void advance() = 0;

  //! Interface - initialization  
  virtual void init() = 0;
  
  //! Interface - solution
  virtual void solve() = 0;

  //! Configure initial conditions of the cells
  virtual void initial_conditions() = 0;

  //! Return an array with the cell types
  const arma::ivec & get_cell_types() const { return cells->get_cell_types(); }

  //! Return reference to the Cells object
  const Cells & get_cells() { return *cells; }

  Stimuli & get_stimuli() { return stimuli; }

  //! Return reference to the Mesh object
  const Mesh & get_mesh() { return *mesh; }

  //! Return pointer to the Mesh object
  const Mesh * get_mesh_ptr() { return mesh; }

  //! Reference to time parameters object
  const TimeParameters & get_time_parameters() { return tip; }

  //! Return the fiber direction of a given element
  const arma::vec3 & get_fiber(int elindex)
  { return mesh->get_element(elindex).get_fiber(); }

  //! Return the sheet/transverse direction of a given element
  const arma::vec3 & get_trans(int elindex)
  { return mesh->get_element(elindex).get_trans(); }

  //! Return the sheet-normal direction of a given element
  const arma::vec3 & get_normal(int elindex)
  { return mesh->get_element(elindex).get_normal(); }

  //! Change PDE time step
  void set_timestep (double s) { timestep = s; }

  //! Change total time
  void set_totaltime(double t) { totaltime = t; }

  //! Change print rate
  void set_printrate(double r) { printrate = r; }

  //! Configure cardiac problem with given data
  void setup(std::string & b, std::string & c, std::string & m,
             double dt, double T, double pr, double pa);

  //! Configure the type of each cell
  void setup_types(std::string & f);

  //! Need to update stiffness matrix?
  void update_matrix(bool r) { re_assembly_mats = r; }
  
  //! Write solution to output file
  void write_data(const arma::vec & u, const std::string & s, int * step);

  //! Write solution and displacement vector to output file
  void write_data(const arma::vec & u, const arma::vec & displ, const std::string & s, int step);

  //! Write Vm to a text file for post-processing
  void write_data_text(const arma::vec & vm, int * step);
  
  //! To measure times
  TimerSection timer;
 
protected:

  bool re_assembly_mats;
  double timestep;
  double totaltime;
  double printrate;
  double printrate_apd;

  string cell_name;
  string odesolver;
  string mesh_filename;
  string stimuli_filename;

  Stimuli stimuli;
  Parameters parameters;
  TimeParameters tip;
  Cells * cells;
  CellModel * cellmodel;
  Mesh * mesh;
  WriterHDF5 * writer;

};


#endif /* CARDIAC_PROBLEM_HPP_ */