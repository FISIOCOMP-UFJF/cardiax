#ifndef MONOPURKINJE_H
#define MONOPURKINJE_H

#include <mesh/purkinje_mesh.hpp>
#include "cardiacproblem.hpp"
#include "fem/fem.h"

/** 
    This class implements the Monodomain model for the Purkinje system
		using the Cubic Hermite Finite Element.
		
    Conductivities are calculated as follows:
      sigma_k = (sigma_k_i * sigma_k_e)/(sigma_k_i + sigma_k_e)
    where k = l, t or n, denoting longitudinal, transverse and normal
*/

class MonodomainPurkinje : public CardiacProblem
{
public:

  //! Constructor
  MonodomainPurkinje();

  //! Destructor
  virtual ~MonodomainPurkinje();

  //! Perform one time step
  void advance();

  //! Initialize solver (mats, vecs, ics, etc)
  void init();

  //! Set initial conditions on cells
  void initial_conditions();

  //! Solve the problem
  void solve();

protected:

  //! Number of degrees of freedom 
  uint ndofs;

  //! Intracellular conductivity in longitudinal direction (mS/um)
  double sigma_l;

  //! Purkinje fiber radius (um: micrometer)
  double radius;

  //! Nodes to apply stimulus
  std::set<uint> stim_nodes;

  //! Bifurcation nodes
  std::vector<int> bif_dofs;

  //! Finite element space of the Purkinje System solution
	CubicHermiteFESpace pk_fespace;

  //! Vector of the transmembrane potential at current time
  arma::vec vm;
   
  //! Assembled global mass and stiffness matrices
  petsc::Matrix Mi, Ai;
  
    //! Vector of solutions
  petsc::Vector v0, v1, f;

  //! PETSc linear solver handler
  petsc::LinearSolver solver;

  //! Assemble mass and stiffness matrix for solving monodomain
  virtual void assemble_matrices();

  //! Compute element stiffness matrix
  virtual void calc_elmat_stiff(const int eindex, const FiniteElement & fe,
                                arma::mat & elmat);

  //! Compute element mass matrix
  void calc_elmat_mass(const int eindex, const FiniteElement & fe, 
											 arma::mat & elmat);

  //! Advance systems of ODEs in time
  void solve_odes();

  //! Solve parabolic problem 
  virtual void solve_parabolic();
};

#endif


