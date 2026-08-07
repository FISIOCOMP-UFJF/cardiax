#ifndef CELLS_H
#define CELLS_H

#include <set>

#include "linalg/linalg.h"
#include "linalg/petsc_vector.hpp"
#include "odes/cellmodel.hpp"
#include <armadillo>

using namespace std;

/** 
    Cells represents a collection of systems of ODES that describes
    the kinetics of a given cell model. 
    
    TODO: so far it uses only one cell model, we need to improve this
    and implement it to support different cell models.
*/
class Cells
{
 public:
 
  //! Default constructor
  Cells (uint n, CellModel * c); 
  
  //! Default destructor
  ~Cells();
   
  //! Set initial conditions to all system of ODEs
  void init();

  //! Advance systems of ODEs in time
  void advance(double t, double dt);

  //! Advance systems of ODEs in time with given Istim data (region)
  void advance(double t, double dt, const double istim, const std::set<uint> & snodes);

  //! Advance systems of ODEs in time with given Istim data (nodes)
  void advance(double t, double dt, const arma::vec & stim_values);

  //! Loop in time to solve the systems of ODEs
  void solve();

  //! Return an array with the cell types
  const arma::ivec & get_cell_types() const { return types; }
  
  //! Return the current time
  inline double get_time() const { return ts->time(); }

  //! Return the value in position i of the global states variable array
  inline double get_state(uint i) const { return states[i]; }

  //! Return the value of variable i in system s
  double get_state(uint s, uint i) const;

  //! Return the entire state variables array
  inline double * get_state_vars() const { return states; }

  //! Return the number of state variables of the ODE
  inline uint get_ode_size() const { return ode->get_num_state_vars(); }

  //! Return the size = number of systems of ODE * number of state vars
  inline uint get_size() const { return num_systems*ode->get_num_state_vars(); }

  //! Return an array with all the values of a given variable
  void get_var(int vindex, double * varray) const;

  //! Return an array with all the values of a given variable
  void get_var(int vindex, arma::vec &v) const;

  // Return an array with all the values of a given variable
  void get_var(int vindex, petsc::Vector &v) const;

  //! Return an array with all monitored values
  void get_monitored_values(int mindex, arma::vec & v) const;

  //! Set state variable vindex with contents of varray
  void set_var(int vindex, double * varray) const;

  //! Set state variable vindex with contents of v
  void set_var(int vindex, arma::vec & v) const;

  //! Config types
  void set_cell_types(int num, int * vtypes);

  //! Return the number of cells (systems of ODEs)
  int size() { return num_systems; }

 protected:  
  
  //! Number of cells (systems of ODEs)
  uint num_systems;

  //! The ODE system describing the cell model
  CellModel* ode;  

  //! Timing step
  TimeStepper * ts;

  //! Array to store the state variables of each ODE system
  double * states;

  //! Types for each cell
  arma::ivec types;

  //! Monitored values array
  arma::vec monitored_values;

};

#endif
