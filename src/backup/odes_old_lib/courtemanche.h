#ifndef _COURTEMANCHE_1998_H_
#define _COURTEMANCHE_1998_H_

#include <cmath>

#define NUMBER_EQUATIONS_CELL_MODEL 20

typedef float Real;

void setIC_ode_cpu (Real *sv, int cellID, int NEQ);

void solve_ode_cpu (Real time, Real dt, Real *sv, int cellID, int NEQ, int col);

void solve_Forward_Euler_cpu (Real time, Real *sv, Real dt, int NEQ, int cellID, int col);

void RHS_Courtemanche_cpu(Real time, Real *rY, Real *rDY, int col);

#endif

