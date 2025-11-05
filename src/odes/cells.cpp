#include "cells.hpp"

Cells::Cells(uint n, CellModel *c) 
  : num_systems(n), ode(c), types()
{
  // allocate memory
  uint mem = c->get_num_state_vars() * num_systems;
  states = new double[mem];

  // initialize timestepper with the one from the given cellmodel 
  ts = ode->get_timestepper();

  // initialize monitored values array
  monitored_values.resize(ode->get_num_monitored() * num_systems);
}

Cells::~Cells()
{
  delete [] states;
}

void Cells::advance(double t, double dt)
{       
  for (uint system=0; system<num_systems; system++)    
  {
    // Compute offset for ODE
    const uint offset = system*ode->get_num_state_vars();

    if (types.size() != 0) ode->set_celltype( types(system) );

    // Time-stepping
    ode->advance(states+offset, t, dt);

    // Update monitored values
    if (ode->get_num_monitored() > 0)
    {
      for(int j=0; j<ode->get_num_monitored(); j++)
      {
        double value = ode->get_monitored_value(j);
        // BUG
        //const uint moffset = system * ode->get_num_monitored() + j;
        const uint moffset = system * ode->get_num_monitored();
        monitored_values(moffset) = value;
      }
    }
  }
}

void Cells::advance(double t, double dt, const double istim,
		                const std::set<uint> & snodes)
{
  std::set<uint>::iterator it;
    
  for (uint system=0; system<num_systems; system++)
  {
    // compute offset for ODE
    const uint offset = system * ode->get_num_state_vars();

    // searching for node I in snodes system
    it = snodes.find(system);

    if (types.size() != 0) ode->set_celltype( types(system) );

    // time-stepping
    if (it != snodes.end())
      ode->advance(states+offset, t, dt, istim);
    else
      ode->advance(states+offset, t, dt);

    // update monitored values
    if (ode->get_num_monitored() > 0)
    {
      for(int j=0; j<ode->get_num_monitored(); j++)
      {
        double value = ode->get_monitored_value(j);
        const uint moffset = system * ode->get_num_monitored();
        monitored_values(moffset) = value;
      }
    }
  }
}

void Cells::advance(double t, double dt, const arma::vec & stim_values)
{
  for (uint system=0; system<num_systems; system++)
  {
    // compute offset for ODE
    const uint offset = system * ode->get_num_state_vars();

    if (types.size() != 0) ode->set_celltype( types(system) );

    const double istim = stim_values(system);
    ode->advance(states+offset, t, dt, istim);

    // update monitored values
    if (ode->get_num_monitored() > 0)
    {
      for(int j=0; j<ode->get_num_monitored(); j++)
      {
        double value = ode->get_monitored_value(j);
        const uint moffset = system * ode->get_num_monitored();
        monitored_values(moffset) = value;
      }
    }
  }
}

void Cells::get_monitored_values(int mindex, arma::vec &v) const
{
  uint mstart = mindex * num_systems;
  uint mend = mstart + num_systems;
  for(uint i=mstart; i<mend; i++)
    v(i) = monitored_values(i);
}

double Cells::get_state(uint s, uint i) const
{
  assert(ode);
  const uint offset = s * ode->get_num_state_vars();
  return states[offset + i];
}

void Cells::get_var(int vindex, double *varray) const
{
  uint odesize = ode->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    varray[i] = states[vindex+(i*odesize)];
}

void Cells::get_var(int vindex, arma::vec &v) const
{
  uint odesize = ode->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    v(i) = states[vindex+(i*odesize)];
}

void Cells::get_var(int vindex, petsc::Vector &v) const
{
  uint odesize = ode->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    v.set(i , states[vindex+(i*odesize)] );
}

void Cells::init()
{
  cout << "Setting initial conditions" << endl;
  
  // Setup initial conditions for each system of ODEs (CellModel)
  for (uint system=0; system<num_systems; system++)
  {        
    // compute offset for ODE
    const uint offset = system*ode->get_num_state_vars();
    
    // change type of cell model (endo, mid, epi)
    if (types.size() != 0) ode->set_celltype( types(system) );

    // set initial conditions of this system
    ode->init(states+offset);
  }
}

void Cells::set_var(int vindex, double *varray) const
{
  uint odesize = ode->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    states[vindex+(i*odesize)] = varray[i];
}

void Cells::set_var(int vindex, arma::vec &v) const
{
  uint odesize = ode->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    states[vindex+(i*odesize)] = v(i);
}

void Cells::set_cell_types(int num, int * vec)
{
  types.set_size(num);

  for(int i=0; i<num; i++)
    types(i) = vec[i];
}

void Cells::solve()
{
  double t;
  double dt = ts->timestep();

  // Set initial conditions to all cells
  init();
   
  cout << " Solving " << num_systems << " cell models" << endl;

  // Start solving
  while(!(ts->finished()))
  {
    ts->increase_time();
    t = ts->time();
    
    // Write to screen
    if(ts->time_to_print())
    {
      cout << "  at time " << t << endl;
      cout << states[0] << "\t" << states[8] << "\t" << states[16] << endl;
    }

    advance(t,dt);
  }

  cout << "\nDone.\n" << endl;

}


