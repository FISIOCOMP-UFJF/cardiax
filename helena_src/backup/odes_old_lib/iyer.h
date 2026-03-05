#ifndef _IYER_H_
#define _IYER_H_

#define NUMBER_EQUATIONS_CELL_MODEL 67

#include <cmath>

// Precistion to be used for the calculations
typedef float Real;

void
setIC_ode_cpu (Real *sv, int cellID, int NEQ);

void
solve_ode_cpu (Real time, Real dt, Real *sv, int cellID, int NEQ, int col);

void
solve_Forward_Euler_cpu (Real time, Real *sv, Real dt, int NEQ, 
		  	 int cellID, int col);

void
RHS_Iyer_cpu (Real time, Real *rY, Real *rDY, int col);

#endif
