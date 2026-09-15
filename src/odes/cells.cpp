#include "cells.hpp"

Cells::Cells(uint n, CellModel *c) 
  : num_systems(n), types()
{
  ode = new CellModel*[1];
  ode[0] = c;

  // allocate memory
  uint mem = c->get_num_state_vars() * num_systems;
  states = new double[mem];

  // initialize timestepper with the one from the given cellmodel 
  ts = ode[0]->get_timestepper();

  // initialize monitored values array
  monitored_values.resize(ode[0]->get_num_monitored() * num_systems);

}

Cells::Cells(uint n, string cell_model_name, string ode_solver, double dt, double totaltime, double tp)
    : num_systems(n), types()
{
  // Isso ficou um pouco confuso, tem a funçaõ omp_get_num_threads, mas parece que ela tem o risco de 
  // retornar um valor inferior ao número de threads (se o numero de threads não for setado, tenho que testar isso). 
  // Mas deixando max_threads, com certeza não vai faltar o objeto ode para nenhuma thread, porém tem o risco de criar 
  // mais objetos do que o necessário!! Tenho que testar isso!! 

    int num_threads = omp_get_max_threads(); 

    ode = new CellModel*[num_threads];

    for (int i = 0; i < num_threads; i++)
    {
        ode[i] = CellModel::create(cell_model_name);

        ode[i]->setup(
            ode_solver,
            dt,
            totaltime,
            tp
        );
    }

    uint mem = ode[0]->get_num_state_vars() * num_systems;
    states = new double[mem];

    ts = ode[0]->get_timestepper();
    monitored_values.resize(
        ode[0]->get_num_monitored() * num_systems
    );
}

Cells::~Cells()
{
  delete [] states;
}

void Cells::advance(double t, double dt)
{
  // convert solver time into the model's own time unit
  // const double tf = time_factor();
  // t  *= tf;
  // dt *= tf;

  // for (uint system=0; system<num_systems; system++)    
  // {
  //   // Compute offset for ODE
  //   const uint offset = system*ode[0]->get_num_state_vars();

  //   if (types.size() != 0) ode[0]->set_celltype( types(system) );
  //   if (apicobasal.n_elem != 0) ode[0]->set_apicobasal( apicobasal(system) );
  //   if (stretch.n_elem != 0)      ode[0]->set_stretch( stretch(system) );
  //   if (stretch_rate.n_elem != 0) ode[0]->set_stretch_rate( stretch_rate(system) );

  //   // Time-stepping
  //   ode[0]->advance(states+offset, t, dt);

  //   // Update monitored values
  //   if (ode[0]->get_num_monitored() > 0)
  //   {
  //     for(int j=0; j<ode[0]->get_num_monitored(); j++)
  //     {
  //       double value = ode[0]->get_monitored_value(j);
  //       //cout << "el: "<<system<<" - " << value << endl;
  //       // BUG
  //       //const uint moffset = system * ode[0]->get_num_monitored() + j;
  //       const uint moffset = system * ode[0]->get_num_monitored();
  //       monitored_values(moffset) = value;
  //     }
  //   }
  // }
  cout<<"This function should be legacy, and if it's been called, there is something wrong! : advance(double, double)" <<endl; 
  exit(10);
}

void Cells::advance(double t, double dt, const double istim,
		                const std::set<uint> & snodes)
{
  const double tf = time_factor();
  t  *= tf;
  dt *= tf;

  #pragma omp parallel for
  for (uint system=0; system<num_systems; system++)
  {
    int tid = omp_get_thread_num(); //thread_id

    // compute offset for ODE
    const uint offset = system * ode[tid]->get_num_state_vars();

    // searching for node I in snodes system
    bool apply_stimulus = snodes.find(system) != snodes.end();

    if (types.size() != 0) ode[tid]->set_celltype( types(system) );
    if (apicobasal.n_elem != 0) ode[tid]->set_apicobasal( apicobasal(system) );
    if (stretch.n_elem != 0)      ode[tid]->set_stretch( stretch(system) );
    if (stretch_rate.n_elem != 0) ode[tid]->set_stretch_rate( stretch_rate(system) );

    // time-stepping
    if (apply_stimulus)
      ode[tid]->advance(states+offset, t, dt, istim);
    else
      ode[tid]->advance(states+offset, t, dt);

    // update monitored values
    if (ode[tid]->get_num_monitored() > 0)
    {
      for(int j=0; j<ode[tid]->get_num_monitored(); j++)
      {
        double value = ode[tid]->get_monitored_value(j);
        const uint moffset = system * ode[tid]->get_num_monitored() + j;
        monitored_values(moffset) = value;
      }
    }
  }
}

void Cells::advance(double t, double dt, const arma::vec & stim_values)
{
  const double tf = time_factor();
  t  *= tf;
  dt *= tf;

  for (uint system=0; system<num_systems; system++)
  {
    // compute offset for ODE
    const uint offset = system * ode[0]->get_num_state_vars();

    if (types.size() != 0) ode[0]->set_celltype( types(system) );
    if (apicobasal.n_elem != 0) ode[0]->set_apicobasal( apicobasal(system) );
    if (stretch.n_elem != 0)      ode[0]->set_stretch( stretch(system) );
    if (stretch_rate.n_elem != 0) ode[0]->set_stretch_rate( stretch_rate(system) );

    const double istim = stim_values(system);
    ode[0]->advance(states+offset, t, dt, istim);

    // update monitored values
    if (ode[0]->get_num_monitored() > 0)
    {
      for(int j=0; j<ode[0]->get_num_monitored(); j++)
      {
        double value = ode[0]->get_monitored_value(j);
        const uint moffset = system * ode[0]->get_num_monitored();
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
  const uint offset = s * ode[0]->get_num_state_vars();
  return states[offset + i];
}

void Cells::get_var(int vindex, double *varray) const
{
  uint odesize = ode[0]->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    varray[i] = states[vindex+(i*odesize)];
}

void Cells::get_var(int vindex, arma::vec &v) const
{
  uint odesize = ode[0]->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    v(i) = states[vindex+(i*odesize)];
}

void Cells::get_var(int vindex, petsc::Vector &v) const
{
  uint odesize = ode[0]->get_num_state_vars();
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
    const uint offset = system*ode[0]->get_num_state_vars();
    
    // change type of cell model (endo, mid, epi)
    if (types.size() != 0) ode[0]->set_celltype( types(system) );
    if (apicobasal.n_elem != 0) ode[0]->set_apicobasal( apicobasal(system) );
    if (stretch.n_elem != 0)      ode[0]->set_stretch( stretch(system) );
    if (stretch_rate.n_elem != 0) ode[0]->set_stretch_rate( stretch_rate(system) );

    // set initial conditions of this system
    ode[0]->init(states+offset);
  }
}

void Cells::set_var(int vindex, double *varray) const
{
  uint odesize = ode[0]->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    states[vindex+(i*odesize)] = varray[i];
}

void Cells::set_var(int vindex, arma::vec &v) const
{
  uint odesize = ode[0]->get_num_state_vars();
  for(uint i=0; i<num_systems; i++)
    states[vindex+(i*odesize)] = v(i);
}

void Cells::set_system_state(uint s, const double * y)
{
  assert(ode);
  assert(y);
  assert(s < num_systems);

  const uint n = ode[0]->get_num_state_vars();
  double * dst = states + s * n;
  for (uint i = 0; i < n; i++) dst[i] = y[i];
}

void Cells::get_system_state(uint s, double * y) const
{
  assert(ode);
  assert(y);
  assert(s < num_systems);

  const uint n = ode[0]->get_num_state_vars();
  const double * src = states + s * n;
  for (uint i = 0; i < n; i++) y[i] = src[i];
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


