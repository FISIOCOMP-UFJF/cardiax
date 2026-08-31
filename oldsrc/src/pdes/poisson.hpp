#ifndef POISSON_HPP
#define POISSON_HPP

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
class Poisson
{
public:

  //! Default constructor
  Poisson();

  //! Destructor
  virtual ~Poisson();

  //! Compute mass matrix in order to compute L2 error norm
  void calc_mass_matrix(arma::sp_mat & mass);

  //! Compute L2 error
  double calc_l2_error(ScalarFunction<double> & exact);  
  
  //! Read parameters (boundary conditions, coefficients) from file
  void config(std::string optFile);

  //! Return reference to mesh
  const Mesh & get_mesh() const { return mesh; }

  //! Return reference to solution vector
  const arma::vec & get_solution() const { return solution; }

  //! Set Neumann boundary conditions map
  void set_neumann(const std::map<int,double> & nm);

  //! Set Dirichlet boundary conditions map
  void set_dirichlet(const std::map<int,double> & dm);

  //! Set fixed nodes (Dirichlet) boundary conditions map
  void set_fixed_nodes(const std::map<int,double> & dm);

  //! Output VTU file for visualization in Paraview.
  void write_data(const string & filename);

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

  //! Auxiliary vector to hold solution
  arma::vec solution;

  //! PETSc system matrix
  petsc::Matrix K;

  //! PETSc system vectors. u: solution, f: right-hand side
  petsc::Vector u,f;

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
  virtual void assemble_system();

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
  
  //! Computes the vector associated with boundary conditions
  void calc_robin_elvec (const int iel, const FiniteElement & fe,
                         const ScalarFunction<double> & gD,
                         const ScalarFunction<double> & gN,
                         arma::vec & elvec);

  //! Solve linear system of equations using PETSc
  void solve();
  
};

// ------------------- inline and template functions --------------------------

inline Poisson::Poisson()
  : writer(&mesh), has_cond(false)
{
  //vtkout
  // do nothing
}

inline Poisson::~Poisson()
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
    //return 2.0;
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
    //return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    return 1.0;
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
