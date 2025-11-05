#ifndef _MAHAJAN_2008_H_
#define _MAHAJAN_2008_H_
#include <cmath>

// Precision to be used for the calculations
typedef float Real;

#define NUMBER_EQUATIONS_CELL_MODEL 26

void 
setIC_ode_cpu(Real *sv, int cellID, int NEQ);

void 
solve_ode_cpu(Real time, Real dt, Real *sv, 
	      int cellID, int NEQ, int col);

void
solve_Forward_Euler_cpu(Real time, Real *sv, Real dt, int NEQ, 
			int cellID, int col);

void 
RHS_mahajan_cpu(Real time_, Real *sv_, Real *rDY_, int col);

#endif

