/**
 * Stimulus (.stim)
 * Format:
 *    num_stimuli
 *    it val x0 x1 y0 y1 z0 z1
 *    it val x0 x1 y0 y1 z0 z1
 *    ...
 *    it val x0 x1 y0 y1 z0 z1
 *
 * Note: use it ordered in time
 */

#include "stimulus.h"
#include <iostream>
#include <fstream>

using namespace std;

void Stimuli::read(char * filename)
{
  int it;
  Real time;  
  Real x0, x1, y0, y1, z0, z1, strength;
  char sfile[256];

  sprintf(sfile, "%s.stim", filename);

  ifstream in(sfile);
  in >> num_stimuli;
  svec.resize(num_stimuli);

  tstim = num_stimuli;

  for(int i=0; i<num_stimuli; i++)
  {
    in >> it >> strength >> x0 >> y0 >> z0 >> x1 >> y1 >> z1;
    svec[i] = new Stimulus(it, strength, x0, x1, y0, y1, z0, z1);
  }

  in.close();
}

void Stimuli::check(int it, Real time, int neq, Real *sv, const Mesh &msh)
{
  int i, ndof;
  Real **nodes;
  ndof  = msh.n_points;
  nodes = msh.nodes;

  if(num_stimuli != 0)
  {
    for(i=astim; i<tstim; i++)
    {
      if (it == (svec[i]->it)+1)
      {
	cout << " Applying stimulus at time " << time << endl;
	apply(neq, sv, svec[i], ndof, nodes);	
	astim++;
	num_stimuli--;
      }
    }
  }
}

void Stimuli::apply(int neq, Real *sv, Stimulus *s, int ndof, Real **nodes)
{
  int i;
  int index, sindex;
  Real x,y,z;

  for(i=0; i<ndof; i++)
  {
    x = nodes[0][i];
    y = nodes[1][i];
    z = 0.0;  // nodes[2][i];

    if( ((x>=s->x0)&&(x<=s->x1)) && ((y>=s->y0)&&(y<=s->y1)) && 
        ((z>=s->z0)&&(z<=s->z1)) )
      index = i; 
    else
      index = -1;

    if(index != -1)
    {
      sindex = neq * index;
      sv[sindex] = s->value;
    }
  }
  cout << endl;
}





Real **
stimulus_read (char *filename, int *nstim)
{
  int i, it;
  Real x0, x1, y0, y1, z0, z1, strength;
  Real **stim_array;
  char sfile[256];

  sprintf(sfile, "%s.stim", filename);

  ifstream in(sfile);
  in >> *nstim;

  stim_array = new_2D_array<Real> (*nstim, 8);

  for(i=0; i<*nstim; i++){
 
    in >> it >> strength >> x0 >> y0 >> z0 >> x1 >> y1 >> z1;

    stim_array[i][0] = it;
    stim_array[i][1] = strength;
    stim_array[i][2] = x0;
    stim_array[i][3] = y0;
    stim_array[i][4] = z0;
    stim_array[i][5] = x1;
    stim_array[i][6] = y1;
    stim_array[i][7] = z1;
  }
  
  in.close();

  return stim_array;
}

void
stimulus_apply (int neq, Real *sv, Real *stim_data,
		int ndof, Real **nodes)
{
  int i;
  int index, sindex;
  Real x, y, z;

  // fetch information about stimulus from stim_data
  Real val = stim_data[1];
  Real x0 = stim_data[2];
  Real y0 = stim_data[3];
  Real z0 = stim_data[4];
  Real x1 = stim_data[5];
  Real y1 = stim_data[6];
  Real z1 = stim_data[7];

  for(i=0; i<ndof; i++)
  {
    x = nodes[0][i];
    y = nodes[1][i];
    z = 0.0;  // nodes[2][i];

    if(((x>=x0)&&(x<=x1)) && ((y>=y0)&&(y<=y1)) && ((z>=z0)&&(z<=z1)))
    {
      index = i; 
    } else {
      index = -1;
    }

    if(index != -1)
    {
      sindex = neq * index;
      sv[sindex] = val;
    }
  }
  cout << endl;
}

void stimulus_check (int *nstim, int *astim, int *tstim, 
		     Real **stimv, int neq, Real *sv, 
		     int it, Real time, const Mesh &msh)
{
  int i, ndof;
  Real **nodes;
  ndof  = msh.n_points;
  nodes = msh.nodes;

  if(*nstim != 0){  
    
    for(i=*astim; i<*tstim; i++){

      if (it == stimv[i][0]+1){
	fprintf(stdout," Applying stimulus at time = %f ...\n",time);

	stimulus_apply (neq, sv, stimv[i], ndof, nodes);

        *nstim--; 
        *astim++;
      }
    }
  }
}

