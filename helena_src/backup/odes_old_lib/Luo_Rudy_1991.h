#ifndef _LUO_RUDY_1991_H_
#define _LUO_RUDY_1991_H_

#include <cmath>
#include "datatype.h"
#include <iostream>

#define NUMBER_EQUATIONS_CELL_MODEL 8

// Precistion to be used for the calculations
//typedef float real;

void setIC_ode_cpu(Real *sv, int cellID, int NEQ);
void solve_ode_cpu(Real time, Real dt, Real *sv, int cellID, int NEQ, int col);
void solve_Forward_Euler_cpu(Real time, Real *sv, Real dt, int NEQ, int cellID, int col);
void RHS_Luo_Rudy_1991_cpu(Real time, Real *rY, Real *rDY, int col);

#endif

