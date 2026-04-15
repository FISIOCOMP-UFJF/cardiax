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

  //! Solve the problem
  void solve();

  //! Solve the problem with meshfile given
  void solve(const string &mshfile);

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


