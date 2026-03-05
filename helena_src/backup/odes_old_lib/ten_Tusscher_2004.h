#ifndef _TEN_TUSSCHER_2004_H_
#define _TEN_TUSSCHER_2004_H_

#define NUMBER_EQUATIONS_CELL_MODEL 17

#include <cmath>

// Precision to be used for the calculations
typedef float Real;

void setIC_ode_cpu(Real *sv, int cellID, int NEQ);
void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, int NEQ, int col);
void solve_Rush_Larsen_cpu(Real time, Real *sv, Real dt, int NEQ, int cellID, int col);

#endif
