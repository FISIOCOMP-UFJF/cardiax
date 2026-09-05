#ifndef EIKONAL_H
#define EIKONAL_H

#include "cardiacproblem.hpp"
#include "fem/fem.h"
#include "util/pugixml.hpp"

/** 
    This class implements the Monodomain model.
    Conductivities are calculated as follows:
    sigma_k = (sigma_k_i * sigma_k_e)/(sigma_k_i + sigma_k_e)
    where k = l, t or n, denoting longitudinal, transverse and normal
*/
struct EikonalNode {
    double cost;
    int id;
    bool operator>(const EikonalNode& other) const {
        return cost > other.cost; 
    }
};

class Eikonal : public CardiacProblem
{
public:

  enum CondTensorType { 
      S_ISOTROPIC,    // spatial isotropy
			S_TRANSVERSE,   // spatial transverse isotropy
			S_ORTHOTROPIC,  // spatial orthotropy
			M_ISOTROPIC,    // material isotropy
			M_TRANSVERSE,   // material transverse isotropy
			M_ORTHOTROPIC   // material orthotropy
  };

  //! Constructor
  Eikonal();

  //! Destructor
  virtual ~Eikonal();

  //! Perform one time step
  void advance() override;
  
  //! Change conductivity tensor
  void set_conductivity(int cond);

  //! Initialize solver (mats, vecs, ics, etc)
  void init();

  //! Set initial conditions on cells
  void initial_conditions();

  //! Set cell state variable value
  void set_stimulus_value(int index, double val);

  //! Stimulate every node whose local activation time has been reached,
  //! for `duration` after it. Needed by ionic cell models (ToRORdLand),
  //! which require a depolarising current rather than reading the LAT.
  //! Both lat and duration are in the solver time unit.
  //!
  //! `period` > 0 makes the stimulus repeat once per beat: the activation
  //! window is then tested against t modulo period, so a multi-beat run
  //! re-excites the tissue instead of stimulating only the first beat.
  //! It must be given in the SOLVER time unit, like lat and duration.
  void apply_lat_stimulus(double amplitude, double duration,
                          double period = 0.0);

  //! Local activation times (solver time unit)
  const arma::vec & get_lat() const { return lat; }

  //! Solve the problem
  void solve();

  //! Solve the problem with meshfile given
  void solve(const string &mshfile);

  //! Dijkstra
  void solve_dijkstra(const std::vector<int>& root_nodes, 
                              const std::vector<double>& root_times,
                              const std::vector<std::vector<std::pair<int, double>>>& adj_cost);


   //! Setup
   void setup(std::string & b, std::string & c, std::string & m,
             double dt, double T, double pr, double pa);

protected:
  
  //! Conductivity type
  CondTensorType condtype;

  //! Number of degrees of freedom 
  uint ndofs;

  //! Stimulus value
  double stim_val;

  //! Stimulus control
  bool stim_apply;

  //! Nodal stimulus control
  bool stim_apply_nodes;

  //! Nodes to apply stimulus
  std::set<uint> stim_nodes;

  //! Finite element space of the solution
  H1FESpace fespace;

  //! Vector local activation time
  arma::vec lat;

  //! Nodal stimuli values
  arma::vec stim_values;
   
   //! PETSc linear solver handler
  petsc::LinearSolver solver;

  //! Advance systems of ODEs in time
  void solve_odes();

};

#endif

