#include "apd.h"

double apd_avg(int n, int ix, double **apd)
{
  double avg = 0.0;
  for(int j=0; j<n; j++)
    avg += apd[j][ix];
  return avg/(double)n;
}

void calc_apd(double dt, int n, double *vm, double x, double  * apd, 
              double *act, double * rep)
{
  double ti = 0.0;
  double tf = (n-1)*dt;
  double *t = new double[n];
  for(int i=0; i<n; i++) t[i] = i*dt;

  int idmax = 0;
  double dv, dvmax = 0.0;
  for(int i=0; i<n-1; i++)
  {
    dv = vm[i+1] - vm[i];
    if(dv > dvmax)
    {
      dvmax = dv;
      idmax = i;
    }
  }

  double act_time = t[idmax];
  act_time = t[idmax];

  int imax, imin;
  double vmax, vmin;
  min(n,vm,&imin,&vmin);
  max(n,vm,&imax,&vmax);

  double amp = vmax - vmin;
  double vmx = (1.0 - x/100.0)*amp;

  // Look for vmx in vm
  int iapdx = imax;
  for(int i=imax; i<n; i++)
  {
    if(vm[i] <= vmin + vmx)
    {
      iapdx = i;
      break;
    }
  }

  // calculating the APD at x% of repolarization
  double apdx  = t[iapdx] - t[idmax];
  double vapdx = vm[iapdx];
  double tapdx = t[iapdx];
  
  // output
  *apd = apdx;
  *act = t[idmax];
  *rep = t[iapdx];

#ifdef VERBOSE
  printf("\n");
  printf("min   = %.2f mV\n" , vmin);
  printf("max   = %.2f mV\n" , vmax);
  printf("APA   = %.2f mV\n" , vmax - vmin);
  printf("MUV   = %.2f V/s\n", dvmax);
  printf("APD%d = %.2f ms\n" , (int)x, apdx);
  printf("LAT   = %.2f ms\n" , idmax*dt);
#endif
  delete [] t;
  //return apdx;
}

double calc_rep_time(double dt, int n, double *vm, double x)
{
  double ti = 0.0;
  double tf = (n-1)*dt;
  double *t = new double[n];
  for(int i=0; i<n; i++) t[i] = i*dt;

  int idmax = 0;
  double dv, dvmax = 0.0;
  for(int i=0;i<n-1;i++)
  {
    dv = vm[i+1]-vm[i];
    if(dv > dvmax)
    {
      dvmax = dv;
      idmax = i;
    }
  }

  double act_time = t[idmax];
  act_time = t[idmax];

  int imax, imin;
  double vmax, vmin;
  min(n,vm,&imin,&vmin);
  max(n,vm,&imax,&vmax);

  double amp = vmax - vmin;
  double vmx = (1.0 - x/100.0)*amp;

  // Look for vmx in vm
  int iapdx = imax;
  for(int i=imax; i<n; i++)
  {
    if(vm[i] <= vmin + vmx)
    {
      iapdx = i;
      break;
    }
  }

  // Calculating the APD at x% of repolarization
  double reptime = t[iapdx];
  double vapdx = vm[iapdx];
  double tapdx = t[iapdx];

#ifdef VERBOSE
  cout << endl << "REPtime = " << reptime << " ms" << endl;
#endif
  return reptime;
}

void max(int n, double *v, int *imax, double *vmax)
{
  *imax = 0;
  *vmax = v[0];
  for(int i= 1; i<n; i++)
    if(v[i] > *vmax)
    {
      *vmax = v[i];
      *imax = i;
    }
}

void min(int n, double *v, int *imin, double *vmin)
{
  *imin = 0;
  *vmin = v[0];
  for(int i= 1; i<n; i++)
    if(v[i] < *vmin)
    {
      *vmin = v[i];
      *imin = i;
    }
}
