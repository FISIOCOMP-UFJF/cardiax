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
 
  //! Milliseconds represented by one time unit of the SOLVER.
  //! 1000 when the solver runs in seconds, 1 when it runs in milliseconds.
  //! Time and time step are converted from the solver unit into each
  //! model's own unit before integration, so a model written in ms and one
  //! written in seconds can both be driven by the same solver.
  void set_solver_time_unit_ms(double ms) { solver_time_unit_ms = ms; }
  double get_solver_time_unit_ms() const { return solver_time_unit_ms; }

  //! The underlying cell model
  const CellModel & get_model() const { return *ode; }

  //! Factor converting solver time into the cell model's own time
  double time_factor() const
  { return solver_time_unit_ms / ode->native_time_unit_ms(); }

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

  //! Coordenada apicobasal (Cobiveco ab) de cada celula, 0 apice -> 1 base.
  //! Vetor vazio (o default) desliga o gradiente: o modelo celular fica no
  //! valor neutro e o resultado e identico ao de antes.
  void set_apicobasal(const arma::vec & v) { apicobasal = v; }

  //! Coordenada apicobasal de cada celula
  const arma::vec & get_apicobasal() const { return apicobasal; }

  //! Estiramento de fibra lambda_f = sqrt(I4f) de cada celula, vindo da
  //! mecanica. Vetor vazio (o default) desliga o acoplamento: o modelo
  //! celular fica em lambda = 1 e o resultado e identico ao de antes.
  void set_stretch(const arma::vec & v) { stretch = v; }

  //! Estiramento de fibra de cada celula
  const arma::vec & get_stretch() const { return stretch; }

  //! Taxa d(lambda_f)/dt de cada celula, na unidade de tempo NATIVA do
  //! modelo celular. Vetor vazio = 0 (contracao isometrica).
  void set_stretch_rate(const arma::vec & v) { stretch_rate = v; }

  //! Taxa de estiramento de cada celula
  const arma::vec & get_stretch_rate() const { return stretch_rate; }

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

  //! Coordenada apicobasal de cada celula (vazio = gradiente desligado)
  arma::vec apicobasal;

  //! Estiramento de fibra de cada celula (vazio = acoplamento desligado)
  arma::vec stretch;

  //! Taxa de estiramento de cada celula (vazio = 0)
  arma::vec stretch_rate;

  //! Monitored values array
  arma::vec monitored_values;


 private:

  double solver_time_unit_ms = 1000.0;  //!< solver runs in seconds by default
};

#endif