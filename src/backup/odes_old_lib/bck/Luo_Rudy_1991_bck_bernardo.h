#ifndef LUO_RUDY_1991_H
#define LUO_RUDY_1991_H

#include <math.h>
#include "../datatype.h"

#define NUMBER_EQUATIONS_CELL_MODEL 8

void setIC_ode_cpu (Real *sv, int cellID, int NEQ);
void solve_ode_cpu (Real time, Real dt, Real *sv,
			      int cellID, int NEQ, int col);

#endif
