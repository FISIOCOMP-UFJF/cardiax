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

  //! Milliseconds represented by one time unit of THIS model's equations.
  //! Ionic models (ToRORd, TenTusscher, ...) are written in ms and keep the
  //! default of 1. Phenomenological models written in seconds, such as
  //! Kerckhoffs, override this with 1000. It lets the solver run in whatever
  //! time unit it likes while each model still receives time in its own.
  virtual double native_time_unit_ms() const { return 1.0; }

  //! Index of the state variable that stores the local activation time, for
  //! phenomenological models driven directly by it (Kerckhoffs). Ionic
  //! models return -1: they are driven by a stimulus current, and their
  //! state variables are concentrations and gates that must not be
  //! overwritten with an activation time.
  virtual int lat_var_index() const { return -1; }

  //! Default constructor
  CellModel (int num_states);

  //! Constructor and setup solver
  CellModel (ODESolver *s);

  //! Destructor
  virtual ~CellModel ();

  static void help();
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

  //! Coordenada apicobasal da celula (Cobiveco ab): 0 no apice, 1 na base.
  //! O valor default 0.5 e NEUTRO -- os modelos que usam esse gradiente
  //! devem reduzir-se ao comportamento original em 0.5, para que malhas sem
  //! o bloco <ab> continuem dando exatamente o mesmo resultado de antes.
  void set_apicobasal(double a) { apicobasal = a; }

  //! Coordenada apicobasal da celula
  double get_apicobasal() const { return apicobasal; }

  //! Estiramento na direcao da fibra, lambda_f = sqrt(I4f), vindo da
  //! mecanica. O valor default 1.0 e NEUTRO: sem acoplamento mecanico o
  //! modelo celular se comporta exatamente como antes.
  void set_stretch(double l) { stretch = l; }

  //! Estiramento na direcao da fibra
  double get_stretch() const { return stretch; }

  //! Taxa de estiramento d(lambda_f)/dt, na unidade de tempo NATIVA do
  //! modelo. Default 0.0 = neutro (contracao isometrica do ponto de vista
  //! do modelo celular).
  void set_stretch_rate(double dl) { stretch_rate = dl; }

  //! Taxa de estiramento na direcao da fibra
  double get_stretch_rate() const { return stretch_rate; }

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
  
  //! Arrays for RL data - filled for rlvars indices only
  std::vector<double> rl_inf, rl_tau; 

  // is_rl[i] == true  <=>  i is a Rush-Larsen variable
  std::vector<bool> is_rl;   
  
  double dt_solver = 0.0;

 protected:

  //! Number of state variables (equations)
  const int num_state_vars;

  //! Select cell type
  CellType type;

  //! Coordenada apicobasal da celula (0 apice -> 1 base); 0.5 = neutro
  double apicobasal = 0.5;

  //! Estiramento na direcao da fibra; 1.0 = neutro (sem deformacao)
  double stretch = 1.0;

  //! Taxa de estiramento na direcao da fibra; 0.0 = neutro
  double stretch_rate = 0.0;

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
