#ifndef ODEPROBLEM_H
#define ODEPROBLEM_H

#include <set>
#include <map>
#include <armadillo>
#include "timestepper.h"

using namespace std;

// Forward declarations
class ODESolver;

/**
   Enumeration to specify the correct cell type
   { 0 EPI, 1 MCELL, 2 ENDO, 3 APEX, 4 BASE }
*/
enum CellType { EPI, MCELL, ENDO, APEX, BASE };

/**
    Abstract base class for the ODE systems that describe a cell model
*/
class CellModel
{
 public:

  //! Default constructor
  CellModel (int num_states);

  //! Constructor and setup solver
  CellModel (ODESolver *s);

  //! Destructor
  virtual ~CellModel ();

  //! Advance ODE in time with ODE solver
  void advance(double * statevars, double t, double dt);

  //! Advance ODE in time with ODE solver
  void advance(double * statevars, double t, double dt, double istim);

  //! Compute finite-difference Jacobian.
  //! Ref: "Numerical Methods for Unconstrained Optimization
  //!       and Nonlinear Equations", Dennis and Schnabel
  void compute_jacobian(double * states, double t, arma::mat & jac);

  //! Factory method for the creation of CellModels
  static CellModel * create(std::string cellname);

  //! Return the cell type
  int get_celltype() const { return type; }

  inline double get_monitored_value(const int index) const { return *(monitored[index]); }

  //! Get number of state variables
  inline int get_num_state_vars() const { return num_state_vars; }

  //! Get number of monitored values
  inline int get_num_monitored() const { return monitored.size(); }

  //! Get time stepper object
  TimeStepper * get_timestepper() const { return ts; }

  //! Change the cell type
  void set_celltype(int ctype);

  //! Set the stimulus
  void set_stimulus(double stim) { i_stim = stim; }

  //! Set timestep and time to print
  void setup(string method, double timestep, double tf, double tp=1.0);

  //! Loop in time to solve ODE
  void solve();

  //! Loop in time to solve a single ODE system
  void solveTest(double stim, double sstart, double sstop,
                 const std::string & fname);

  //! Loop in time to solve a single ODE system (output to HDF5)
  void solveTestHDF5(double stim, double stime, double sdur,
                     const std::string & fname);

  //! Basic Cycle Length protocol (output to HDF5)
  void solveTestBCL(double stim, double sdur, double bcl,
                    const std::string & fname);

  //! Restitution protocol (output to HDF5)
  void solveTestRTT(double stim, double sdur, double itl, double bcl, double delta,
                    const std::string & fname);

  // Interface ----------------------------------------------------------------

  //! Set initial conditions
  virtual void init(double * values) const = 0;

  //! Compute ODE equations
  virtual void equation(const double time, const double * statevars,
                        double * values) = 0;

  //! RL variables
  std::set<int> rlvars;

 protected:

  //! Number of state variables (equations)
  const int num_state_vars;

  //! Select cell type
  CellType type;

  //! Solver
  ODESolver * ode_solver;

  //! Timing step
  TimeStepper * ts;

  //! Stimulus value
  double i_stim;

  //! Used to compute jacobian
  arma::vec f1, f2;

  //! Used to monitor some variables
  std::vector<double*> monitored;

  //! Dictionary for variable names
  std::map<int, std::string> var_names;

};

#endif
