//
// Created by rocha on 1/11/16.
//

#ifndef CARDIAX_PROJECTION_HPP
#define CARDIAX_PROJECTION_HPP

#include <iostream>
#include <map>
#include "fem/fem.h"
#include "util/timer.hpp"
#include "util/input_file.h"
#include "linalg/linalg.h"

using namespace std;

class Projection
{
public:

  //! Default constructor
  Projection() : writer(&msh) {};

  //! Destructor
  virtual ~Projection() {};

  //! Compute mass matrix in order to compute L2 error norm
  void calc_mass_matrix(arma::sp_mat & mass);

  //! Read parameters (boundary conditions, coefficients) from file
  void config(std::string optFile);

  //! Return reference to mesh
  const Mesh & get_mesh() const { return msh; }

  //! Return reference to solution vector
  const arma::vec & get_solution() const { return solution; }

  //! Output VTU file for visualization in Paraview.
  void write_data(const string & filename, bool as_text);

  //! Init, assemble and solve
  void run(const string & filename);

protected:

  //! The finite element mesh
  Mesh msh;

  //! H1 finite element space
  H1FESpace fespace;

  //! Auxiliary vector to hold solution
  arma::vec solution;

  //! PETSc system matrix
  petsc::Matrix K;

  //! PETSc system vectors. u: solution, f: right-hand side
  petsc::Vector u,f;

  //! Output
  WriterHDF5 writer;

  //! Input filename
  string filename;

  //! Read mesh and setup FE space
  void init();

  //! Assemble K u = f
  void assemble_system();


  //! Computes the element mass matrix
  void calc_elmat_mass (const int iel, const FiniteElement & fe,
                          arma::mat & elmat);

  void calc_elmat_mixed_mass(const int ielF, const int ielC,
                            const FiniteElement & feF,
                        const FiniteElement & feC,
                        arma::mat & elmat);

  //! Computes the element 'force' vector
  void calc_elvec_source (const int iel, const FiniteElement & fe,
                            const ScalarFunction<double> & rhs,
                            arma::vec & elvec);

  //! Solve linear system of equations using PETSc
  void solve();

  //! Transfer the FE solution from one mesh to the target mesh
  void transfer(Mesh & target_mesh, arma::vec & target_solution);

};


template <typename T>
class RHSExample1 : public ScalarFunction<T> {
public:
  T operator()(const arma::vec3 &pt) const {
    const T x = pt(0);
    const T y = pt(1);
    return (sin(M_PI * x) * sin(M_PI * y));
    //return 2.0;
  }
};

#endif //CARDIAX_PROJECTION_HPP
