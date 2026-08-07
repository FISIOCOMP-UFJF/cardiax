#ifndef ELASTICITY_HPP_
#define ELASTICITY_HPP_

#include <map>
#include "petscksp.h"
#include "fem/fem.h"
#include "util/timer.hpp"
#include "util/input_file.h"
#include "util/parameters.hpp"
#include "linalg/linalg.h"
#include "linalg/petsc_matrix.hpp"
#include "linalg/petsc_vector.hpp"
#include "linalg/petsc_linear_solver.hpp"

enum ElasticityType
{ 
  PLANE_STRAIN, 
  PLANE_STRESS, 
  THREE_DIM
};

/** 
    Class to solve Linear Elasticity Problem in terms of
    the Lame parameters. The partial differential equations
    are given by:
    \[
      (2 mu e(u),e(v)) + (lami div u , div v) = (f,v) + (gN, v)
    \]
*/
class Elasticity
{
public:

  //! Default constructor
  Elasticity();

  //! Default destructor
  ~Elasticity();

  //! Configure problem with parameters from file
  void config(std::string & configfile);

  //! Run the solver with the given mesh
  void run(const string & meshfilename);

  //! Set the type of the elasticity problem
  void set_type(ElasticityType type);

  //! Set Neumann boundary conditions
  void set_neumann(const std::map<int,arma::vec3> & nm) { neumann_map = nm; }

  //! Set homogeneous Dirichlet boundary conditions
  void set_dirichlet(const std::multimap<int,int> & dm) { dirichlet_map = dm; }

  //! Set nodes with prescribed displacements
  void set_fixed_nodes(const std::multimap<int,NodalData> & fm) { fixed_nodes_map = fm; }

  //! Output data to file
  void write_data(const string & filename);

protected:

  //! Handle the Elasticity problem parameters (Young modulus, ...)
  Parameters parameters;
  ElasticityType elastype;
  Mesh msh;
  VecH1FESpace fespace;
  WriterHDF5 writer;

  std::string filename;
  std::map<int,arma::vec3> neumann_map;
  std::multimap<int,int> dirichlet_map;
  std::multimap<int,NodalData> fixed_nodes_map;

  arma::vec solution;
  petsc::Matrix K;
  petsc::Vector u,f;

  //! Assemble system of linear equations Ku=f
  void assemble_system();

  //! Compute element stiffness matrix
  void calc_elmat_stiffness (const int eindex, const MixedFiniteElement * fe,
                             arma::mat & elmat);

  //! Compute element load vector
  void calc_elvec_source (const int eindex, const MixedFiniteElement * fe,
                          arma::vec & elvec);

  //! Compute element Robin matrix (boundary condition)
  void calc_robin_elmat (const int eindex, const MixedFiniteElement * fe,
                         arma::mat & elmat);

  //! Compute element Robin load vector (boundary condition)
  void calc_robin_elvec (const int eindex, const MixedFiniteElement * fe,
                         arma::vec & elvec);

  //! Compute stress for a given element
  void calc_element_stress(Mapping &em, const MixedFiniteElement *fe,
                           const arma::vec &u, const arma::vec3 &X,
                           arma::vec &sigma);

  //! Compute FE matrix (gradients)
  void calc_B_matrix (const arma::mat & gradn, arma::mat & B);

  //! Compute FE elasticity matrix (6x6 matrix in the general case)
  void calc_D_matrix (const arma::vec3 & X, int index, arma::mat & D);

  //! Compute FE matrix associated with Robin matrix
  void calc_H_matrix (const arma::vec & shape, arma::mat & H);

  //! Compute FE matrix associated with Robin matrix (coefficient)
  void calc_C_matrix (const arma::vec3 & X, arma::mat & C);

  //! Compute values to impose boundary conditions (matrix contribution - Dirichlet)
  arma::vec3 coeff_robin_mat(int index);

  //! Compute values to impose boundary conditions (load vector contribution - Neumann)
  arma::vec3 coeff_robin_vec(int index);

  //! Compute Mu given Young modulus and Poisson ration
  double eval_mu(double e, double v) const { return e/(2.0*(1+v)); }

  //! Compute Lambda given Young modulus and Poisson ration
  double eval_lambda(double e, double v) const { return (v*e)/((1+v)*(1-2.*v)); }

  //! Get the solution data from a specfic element (with dof numbers - vdnums)
  void get_elvec(const std::vector<int> & vdnums, const arma::vec & u,
                 arma::vec & elu);

  //! Initilization (read mesh)
  void init();

  //! Solve the linear elasticity problem
  void solve();

  //! Calculate Von Mises stress of mixed FE Variable
  void von_mises(const arma::vec & u, arma::vec & sigma);

};

#endif /* ELASTICITY_HPP_ */
