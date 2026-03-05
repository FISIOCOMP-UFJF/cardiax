#ifndef _TEN_TUSSCHER_2006_H_
#define _TEN_TUSSCHER_2006_H_

#define NUMBER_EQUATIONS_CELL_MODEL 19

#include <cmath>
#include "datatype.h"

// Precision to be used for the calculations
//typedef float real;

void setIC_ode_cpu(Real *sv, int cellID, int NEQ);
void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, int NEQ, int col);
void solve_Forward_Euler_cpu(Real time, Real *sv, Real dt, int NEQ, int cellID, int col);
void RHS_ten_Tusscher_2006_cpu(Real time, Real *rY, Real *rDY, int col);

#endif
