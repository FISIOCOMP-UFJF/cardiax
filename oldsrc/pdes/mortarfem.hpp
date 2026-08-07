#ifndef MORTARFEM_HPP
#define MORTARFEM_HPP

#include <map>
#include "fem/fem.h"
#include "util/timer.hpp"
#include "util/input_file.h"
#include "linalg/linalg.h"

/*!
 * A simple class to solve Poisson problem given by:
 *
 *   -\Delta u = f, in \Omega
 *
*/
class Mortarfem
{
public:

  //! Default constructor
  Mortarfem();

  //! Destructor
  virtual ~Mortarfem();

  //! Compute mass matrix in order to compute L2 error norm
  void calc_mass_matrix(arma::sp_mat & mass);

  //! Read parameters (boundary conditions, coefficients) from file
  void config(std::string optFile);

  //! Return reference to mesh
  const Mesh & get_mesh() const { return mesh; }

  //! Return reference to solution vector
  const arma::vec & get_solution() const { return solution; }

  //Return the solution of the primal variable problem
  const arma::vec & get_primal_solution() const { return primal_solution; }

  //Return the map of multiplier to global DGOFs
  const std::map<int, vector<int> > & get_map_multiplier() const { return map_multiplier; }

  //! Set Neumann boundary conditions map
  void set_neumann(const std::map<int,double> & nm);

  //! Working with multiplier position and values to
  //Properly compute action potential
  // void computing_potential_diff(const std::map< int, vector<int> > map_multiplier, arma::vec & lambda_sol);
  void computing_potential_diff(const std::map< int, vector<int> > map_multiplier, int t_step);

  //! Set Dirichlet boundary conditions map
  void set_dirichlet(const std::map<int,double> & dm);

  //! Set fixed nodes (Dirichlet) boundary conditions map
  void set_fixed_nodes(const std::map<int,double> & dm);

  //! Output VTU file for visualization in Paraview.
  void write_data(const string & filename);

  void write_data(const string & filename, arma::vec & u_exact);
  //Fill Block A and primary vector of Global Mat
  void fill_primary_matrix_vector();

  //Fill B block from non-matching mesh interface contributions
  void fill_interfaces_contributions_non_matching(int num_interface_eleme, int ndofs, int t_step);

  //Fill B block from matching mesh interface contributions
  void fill_interfaces_contributions_matching(int num_interface_eleme, int ndofs, int t_step);
  //! Init, assemble and solve
  void run(const string & filename);

protected:

  //! The finite element mesh
  Mesh mesh;

  //! H1 finite element space
  H1FESpace fespace;
  //CubicHermiteFESpace fespace;

  //! Neumann boundary conditions description
  std::map<int,double> neumann_map;

  //! Dirichlet boundary (edge or surface) conditions description
  std::map<int,double> dirichlet_map;

  //! Dirichlet boundary (at the nodes) conditions description
  std::map<int,double> fixed_nodes_map;

  //! Relations between multiplier and global DGOFs
  std::map<int, vector<int> > map_multiplier;

  //! Solution of EDO in average by element for the new iteration
  std::map<int, double > interface_potential;

  //! Map each lagrange mult to the state var
  std::map<int, vector<double> > interface_state_var;

  //! Auxiliary vector to hold solution
  arma::vec solution;
  //! Auxiliary vector to hold primal var solutions
  arma::vec primal_solution;
  //! For solving with direct methods of Armadillo
  arma::vec primary_sol;
  //! Auxiliary vector to hold multiplier solution
  arma::vec lambda_solution;
  ////Used to hold solution and compute transmembrane current
  arma::vec lambda_sol;
  /////to save the source of the jump in interfaces
  arma::vec Vm_vec;
  //This matrix are needed to compute
  //the primal solution as post-processing
  arma::mat A_inv;
  arma::mat global_A;
  arma::mat B;
  arma::mat Bt;
  arma::vec primary_vec;
  //! PETSc system matrix
  petsc::Matrix K, L;
  petsc::Matrix Schur;

  //! PETSc system vectors. u: solution, f: right-hand side
  petsc::Vector u,f,f_lambda,sol_lambda;
  //! PETSc system vectors. u: solution, f: right-hand side
  petsc::Vector lambda, load_Vm;

  //! VTK output
  //WriterVTK vtkout;
  WriterHDF5 writer;

  //! Input filename
  string filename;

  //! Conductivity
  bool has_cond;

  //! Read mesh and setup FE space
  void init();

  //! Assemble K u = f
  virtual void assemble_system(int time_step);

  //! Get Robin boundary condition matrix
  double coeff_robin_mat(int index);

  //! Get Robin boundary condition vector
  double coeff_robin_vec(int index);

  //! Computes the element mass matrix
  void calc_elmat_mass (const int iel, const FiniteElement & fe,
												arma::mat & elmat);

  //! Computes the element stiffness matrix
  void calc_elmat_poisson (const int iel, const FiniteElement & fe,
													 arma::mat & elmat);

    //! Computes the element 'force' vector
  void calc_elvec_source (const int iel, const FiniteElement & fe,
													const ScalarFunction<double> & rhs,
													arma::vec & elvec);

  //! Computes the matrix associated with boundary conditions
  void calc_robin_elmat (const int iel, const FiniteElement & fe,
												 arma::mat & elmat);
  //! Computes the matrix associated with interface conditions
  void calc_interface_elmat (const int iel, const FiniteElement & fe,
                        arma::mat & elmat);
  //! Computes the vector associated with boundary conditions
  void calc_interface_elvec (const int iel, const FiniteElement & fe,
                         arma::vec & elvec, double Vm);

  //! Computes the vector associated with boundary conditions
  void calc_robin_elvec (const int iel, const FiniteElement & fe,
												 const ScalarFunction<double> & gD,
												 const ScalarFunction<double> & gN,
												 arma::vec & elvec);

  int count_interface_elements();

  //! Solve linear system of equations using PETSc
  void solve(int t_step);

  void write_interface_potential(string file_name);

  void write_gradient(const string & filename);

};

// ------------------- inline and template functions --------------------------

inline Mortarfem::Mortarfem()
  : writer(&mesh), has_cond(false)
{
  //vtkout
  // do nothing
}

inline Mortarfem::~Mortarfem()
{
  // do nothing
}

// ------------------- Right hand side functions ------------------------------

template <typename T>
class RHSExample0 : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const {
    return 0.0;
  }
};


template <typename T>
class RHSExampleOne : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const {
    return 1.0;
  }
};


template <typename T>
class RHSExample1 : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const {
    const T x = pt(0);
    const T y = pt(1);
    return (sin(M_PI*x)*sin(M_PI*y));
  }
};

template <typename T>
class RHSExample2 : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const {
    const T x = pt(0);
    const T y = pt(1);
    return -2.0*x*(x-1.0)-2.0*y*(y-1.0);
  }
};

template <typename T>
class RHSExample3D : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const {
    const T x = pt(0);
    const T y = pt(1);
    const T z = pt(2);
    return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
  }
};

template <typename T>
class RHSExample3D2 : public ScalarFunction<T>
{
public:
  T operator()(const arma::vec3 & pt) const {
    const T x = pt(0);
    const T y = pt(1);
    const T z = pt(2);
    return -2.0*x*(x-1.0)-2.0*y*(y-1.0)-2.0*z*(z-1.0);
  }
};

#endif