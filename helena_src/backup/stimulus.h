// Author Bernardo M. Rocha.
//
// First added:  2009-12-01
// Last changed: 2009-09-15
//
// Header file for stimulus interface
//

#ifndef STIMULUS_H
#define STIMULUS_H

#include <iostream>
#include <fstream>
#include "util/util.h"
#include "util/malloc.h"
#include "util/datatype.h"
#include "fem/fem.h"

struct Stimulus
{
  int it;
  Real time;
  Real value;
  Real x0, x1, y0, y1, z0, z1; 
  
  Stimulus(int iter, Real v, Real xb, Real xe, Real yb, Real ye, Real zb, Real ze)
     : it(iter), value(v), x0(xb), x1(xe), y0(yb), y1(ye), z0(zb), z1(ze) {}
};

class Stimuli
{
 private:
  int num_stimuli;
  int tstim;
  int astim;
  vector<Stimulus*> svec;

 public:
  Stimuli(): num_stimuli(0), astim(0), tstim(0) {}
  ~Stimuli(){}

  void read(char * filename);
  void apply(int neq, Real *sv, Stimulus *s, int ndof, Real **nodes);
  void check(int it, Real time, int neq, Real *sv, const Mesh &msh);

};

/** Stimulus functions prototypes */

Real **
stimulus_read (char * filename, int *nstim);


void
stimulus_check (int *nstim, int *astim, int *tstim, 
		Real **stimv, int neq, Real *sv, 
		int it, Real time, const Mesh & msh);

#endif
