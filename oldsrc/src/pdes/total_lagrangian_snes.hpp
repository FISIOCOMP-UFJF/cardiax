#ifndef TOTAL_LAGRANGIAN_SNES_HPP
#define TOTAL_LAGRANGIAN_SNES_HPP

#include "petscsnes.h"
#include "total_lagrangian.hpp"

/** This class implements the Total Lagrangian formulation
    for the solution of nonlinear elasticity problems using
    PETSc SNES (Scalable Nonlinear Equations Solvers)
*/

class TotalLagrangianSNES : public TotalLagrangian
{
public:
  TotalLagrangianSNES(){};
  ~TotalLagrangianSNES(){};
  void init_resid();
  void init_stiff();
  void func();
  void jacob();
  void solve();
  void update_coordinates(PetscScalar *xx);

  petsc::Vector & residual() { return r; }
  petsc::Matrix & jacobian() { return K; }
  
  petsc::Vector x_aux;
};

#endif
