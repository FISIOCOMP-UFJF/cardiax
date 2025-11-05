#ifndef _BONDARENKO_H_
#define _BONDARENKO__H_

#include <cmath>

#define NUMBER_EQUATIONS_CELL_MODEL 40

// Precistion to be used for the calculations
typedef float Real;

void
setIC_ode_cpu (Real *sv, int cellID, int NEQ);

void
solve_ode_cpu (Real time, Real dt, Real *sv, int cellID, int NEQ, int col);

void
solve_Forward_Euler_cpu (Real time, Real *sv, Real dt, int NEQ, int cellID, int col);

void
RHS_Bondarenko_cpu (Real time, Real *rY, Real *rDY, int col);

#endif

