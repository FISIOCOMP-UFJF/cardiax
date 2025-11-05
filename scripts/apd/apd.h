#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>

using namespace std;

//! Calculates the average of the APD in one given direction (ix)
double apd_avg(int n, int ix, double **apd);

//! Calculates the action potential duration at X% of the repolarization
void calc_apd(double dt, int n, double *vm, double x, double  *apd, 
              double *act, double *rep);

//! Calculates the repolarization time
double calc_rep_time(double dt, int n, double *vm, double x);

//! Maximum value of the transmembrane potential
void max(int n, double *v, int *imax, double *vmax);

//! Minimum value of the transmembrane potential
void min(int n, double *v, int *imin, double *vmin);
