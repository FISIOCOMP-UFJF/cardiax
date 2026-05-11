#ifndef BIDOMAIN_HPP_
#define BIDOMAIN_HPP_

#include <unordered_set> 
#include "cardiacproblem.hpp"
#include "fem/fem.h"

/** Concrete class that implements the (uncoupled) Bidomain problem
    The conductivities were obtained from literature.
    References (implementation): 
       1. Joakim Sundnes
          "Computing the electrical activity in the heart"

    References (conductivities):
       1. Bradley J. Roth 
          "Electrical Conductivity Values Used with the Bidomain
	  Model of Cardiac Tissue"

       2. L. S. Graham
	  "Estimation of the Bidomain Conductivity Parameters of
	  Cardiac Tissue From Extracellular Potential Distributions
	  Initiated by Point Stimulation"

   Values for conductivities:

      Conductivities from L. Clerc (units in mS/um)
         sigma_il(0.000170) sigma_it(0.000019)
         sigma_el(0.000620) sigma_et(0.000240)

*/

class Bidomain : public CardiacProblem
{
public:
  
  //! Constructor
  Bidomain();
  
  //! Destructor
  ~Bidomain();

  //! Perform one time step
  virtual void advance(){};

  //! Initialization of the problem
  virtual void init();

  //! Setup up initial conditions of the systems of ODEs. 
  void initial_conditions();
  
  //! Solve the problem in time
  virtual void solve();
  
  //! Reconfigure cell types for bidomain problem where dofs are different
  void setup_types(std::string & f);
  
  //! TODO: document 
  std::map<uint,uint> bidomap;

protected:

  //! Number of state variables in cell model 
  uint neq;

  //! Total number of degrees of freedom (tissue+bath) 
  uint ndof;

  //! Number of Vm degrees of freedom (tissue only) 
  uint ndof_tissue;

  //! Number of Ve degrees of freedom (bath only) 
  uint ndof_bath;

  //! Number of elements in bath region 
  uint num_bath_elems;

  //! Number of elements in tissue region 
  uint num_tissue_elems;

  //! Value for intracellular stimulus. 
  double stimval;

  //! Apply intracellular stimulus or not ? 
  bool intra_apply;

  //! Bidomain intra-/extra-cellular conductivities (mS/cm)
  double sigma_il, sigma_it, sigma_in;
  double sigma_el, sigma_et, sigma_en;

  //! TODO: document 
  std::unordered_set<uint> tissueset, bathset;

  //! Connectivity for Vm+Ve in tissue+bath, reordered but not renumbered 
  std::vector< std::vector<uint> > connec;

  //! Connectivity for Vm using renumbered DoFs 
  std::vector< std::vector<uint> > connec_vm;

  //! Vector of index nodes to apply intracellular stimulus 
  std::set<uint> si_nodes;

  //! Temporary vectors for manipulating data
  arma::vec tmp_vm, tmp_ve;

  //! Wrappers to PETSc linear system solvers and matrices. 
  petsc::LinearSolver par_solver;
  petsc::LinearSolver ell_solver;
  petsc::Matrix Ai, Aie, Mi, Ki;
  petsc::Vector vm0, vm1;
  petsc::Vector ve, b0, b1, be, aux_vec;

  //! Finite element H1 space 
  H1FESpace fespace;

  // Methods -------------------------------------------------------------------

  //! Assemble FEM matrices to solve the problem. 
  void assemble_matrices();

  //! Compute global conductivity tensors for intracellular and extracellular.
  void calc_cond_tensors(const int index, const int ndim,
			 arma::mat & sigma_i, arma::mat & sigma_e);

  //! Compute element stiffness matrix for parabolic and elliptic problems. 
  void calc_elmat_stiff_ie(const int eindex, const FiniteElement & fe,
			       arma::mat & elmat_i, arma::mat & elmat_ie);

  //! Compute element (bath) stiffness matrix of elliptic problem only. 
  void calc_elmat_bath(const int eindex, const FiniteElement & fe,
		       arma::mat & elmat);

  //! Compute element mass matrix for parabolic problem. 
  void calc_elmat_mass(const int eindex, const FiniteElement & fe,
		       arma::mat & elmat);

  //! Fix boundary values to solve singularity of the elliptic problem. 
  void fix_phie_values(std::vector<int> & nodes);

  //! Get the Bidomain numbering of DoFs 
  void get_dof_nums(int i, std::vector<int> & dnums);

  /*! Renumbering of Bidomain degrees of freedom for parabolic 
      and elliptic problems. Creates the tissue_map, tissue_connec, 
      bath_map and bath_connec structures to assemble linear 
      systems and output solution (vm,ve).
  */
  void setup_nodes();
 
  //! Advance systems of ODEs in time 
  void solve_odes();

  //! Solve elliptic problem defined by:
  //!    div((sigma_i + sigma_e) grad(Vm)) = - div(sigma_i grad(vm))
  void solve_elliptic();

  /*! Solve parabolic problem
      b0 = M * vm0  (matrix vector multiplication)
      b  = Ki * ve0 (matrix vector multiplication)
      A * vm1 = b   (linear system solve)
  */
  void solve_parabolic();

  //! Write VTU (VTK) files for tissue+bath region 
  //void write_data_ve(const arma::vec & u, const std::string & s, int * step);

  //! Write VTU (VTK) files for tissue region only 
  //void write_data_vm(const arma::vec & u, const std::string & s, int * step);

};

#endif

