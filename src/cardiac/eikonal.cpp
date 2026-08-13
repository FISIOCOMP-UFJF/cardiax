#include "eikonal.hpp"
#include "monodomain.hpp"

#include <cmath>

/*!
 * Conductivities should be given in mS/um
 * 
 * My old setting: sigma_l(0.00012),   sigma_t(0.00006226), sigma_n(0.00002556)
 * Rodrigo's code: sigma_l(0.0001543), sigma_t(0.0000324),  sigma_n(0.0000120)
 * Benchmark     : sigma_l(0.0001334), sigma_t(0.0000176),  sigma_n(0.0000176)
 * Rabbit        : sigma_l(0.000204),  sigma_t(0.0000102),  sigma_n(0.000037)
 * Sundnes       : sigma_l(0.000300), sigma_t(0.000100), sigma_n(0.000031525)
*/

//#define MONO_ISOTROPIC

Eikonal::Eikonal()
  : CardiacProblem(),
    //sigma_l(0.0001334), sigma_t(0.0000176),  sigma_n(0.0000176),
    stim_apply_nodes(false)
{
  cout << "Eikonal" << endl; 
  mesh = new Mesh();
  writer = new WriterHDF5(mesh);

  //? Precisa disso? 
  parameters.rename("Eikonal_parameters");
  parameters.add("sigma_l", 0.0001334);
  parameters.add("sigma_t", 0.0000176);
  parameters.add("sigma_n", 0.0000176);
}

Eikonal::~Eikonal()
{
  delete cells;
  delete cellmodel;
}

void Eikonal::advance()
{
  if( !tip.finished() )
  {
    tip.increase_time();
    timer.enter("ODEs");
    solve_odes();
    timer.leave();
  }
}

void Eikonal::init()
{
  tip = TimeParameters(timestep, totaltime, printrate);

  mesh->read_xml(mesh_filename);
  stimuli.read_xml(stimuli_filename);

  fespace.set_mesh(mesh);
  ndofs = mesh->get_n_points(); //Before: get_n_points() 

  // Per-node stimulus buffer. Monodomain sizes it in its own init(); Eikonal
  // never did, so it stayed empty and any per-node stimulus indexed out of
  // bounds. It is needed here for LAT-driven stimulation of ionic models.
  stim_values.zeros(ndofs);

  // setup data writer to write at every 1 ms
  // potential field and displacements
  std::size_t pos  = mesh_filename.find(".xml");
  std::string output = mesh_filename.substr(0,pos) + "_output";
  int nsteps = tip.get_size(); 
  writer->open(output, nsteps+1, timestep);

  // setup model and cells
  cellmodel = CellModel::create(cell_name);
  cellmodel->setup(odesolver, timestep, totaltime, 1.0);
  cells = new Cells(ndofs, cellmodel);

  cells->init();
}

void Eikonal::set_conductivity(int cond)
{
    CondTensorType tcond = static_cast<CondTensorType>(cond);
    condtype = tcond;
  
    cout << "Conductivity type: ";
    switch(condtype)
    {
        case S_ISOTROPIC:   cout << "S_ISOTROPIC"   << endl; break;
        case S_TRANSVERSE:  cout << "S_TRANSVERSE"  << endl; break;
        case S_ORTHOTROPIC: cout << "S_ORTHOTROPIC" << endl; break;
        case M_ISOTROPIC:   cout << "M_ISOTROPIC"   << endl; break;
        case M_TRANSVERSE:  cout << "M_TRANSVERSE"  << endl; break;
        case M_ORTHOTROPIC: cout << "M_ORTHOTROPIC" << endl; break;
    }
}

void Eikonal::initial_conditions()
{
  tip.reset();
  
  cells->init();

  // Pre-condicionamento do modelo celular 0D: substitui as condicoes
  // iniciais publicadas pelo estado do ciclo limite no BCL da simulacao.
  // Desligado por default (-warmup 1 liga), e portanto sem efeito nenhum
  // sobre as rodadas existentes.
  //
  // Tem de vir DEPOIS de cells->init() -- que reescreve todo o vetor de
  // estado -- e ANTES da escrita do LAT abaixo, para nao apagar o tempo de
  // ativacao dos modelos fenomenologicos.
  warmup_cells(lat);

  // Only phenomenological models store the activation time as a state
  // variable. Writing the LAT into variable 1 unconditionally destroyed the
  // initial conditions of ionic models -- in ToRORd variable 1 is the
  // intracellular sodium (11.9 mM), so it was being zeroed and the cell
  // could never produce an action potential.
  int lat_idx = cells->get_model().lat_var_index();
  if (lat_idx >= 0)
    cells->set_var(lat_idx, lat);

  // loop in time
  int step=0;

  cells->advance(tip.time(), timestep, stim_val, stim_nodes);
  
}

void Eikonal::set_stimulus_value(int index, double val)
{
  stim_apply_nodes = true;
  stim_values(index) = val;
}

void Eikonal::apply_lat_stimulus(double amplitude, double duration,
                                 double period)
{
  // Stimulate every node whose local activation time has been reached, for
  // `duration` after it. This is what drives ionic cell models such as
  // ToRORdLand, which need an actual depolarising current: unlike the
  // phenomenological Kerckhoffs model, they do not read the activation time
  // directly.
  //
  // lat and duration are in the solver time unit (seconds here, so a 2 ms
  // stimulus is duration = 0.002).
  if (stim_values.n_elem != lat.n_elem)
  {
    std::cout << " Warning: stimulus buffer (" << stim_values.n_elem
              << ") and LAT vector (" << lat.n_elem
              << ") have different sizes; resizing." << std::endl;
    stim_values.zeros(lat.n_elem);
  }

  double t = tip.time();

  // lat is fixed and t grows monotonically, so without folding the time back
  // into a single beat the activation window [lat, lat+duration) is crossed
  // only during the first beat and every later beat runs unstimulated.
  if (period > 0.0) t = std::fmod(t, period);

  bool any = false;

  for (arma::uword i = 0; i < lat.n_elem; i++)
  {
    if (t >= lat(i) && t < lat(i) + duration)
    {
      stim_values(i) = amplitude;
      any = true;
    }
  }

  if (any)
    stim_apply_nodes = true;
}

void Eikonal::solve(const string &mshfile)
{
  cout << " -- Reading local activation time from mesh file --" << endl; 
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());

  // zeros(), not set_size(): nodes missing from <eikonal> would otherwise
  // keep uninitialised memory.
  lat.zeros(ndofs);

  pugi::xml_node pvloop_data = doc.child("mesh").child("pvloop");
  double begin_active_stress = 0.0;
  // passive_time is given in MILLISECONDS in the mesh file, while the solver
  // may work in seconds; ms_to_solver_time() does the conversion.
  if(pvloop_data)
    begin_active_stress = pvloop_data.attribute("passive_time").as_double() * ms_to_solver_time();

  pugi::xml_node eikonal_data = doc.child("mesh").child("eikonal");
  bool has_eikonal = false;
  int n_read = 0;

  if(eikonal_data)
  {
    has_eikonal = true;
    for(pugi::xml_node node = eikonal_data.child("node"); node; node = node.next_sibling("node"))
    {
      int index = node.attribute("id").as_int();
      if(index < 0 || index >= (int) ndofs)
      {
        std::cout << " *** WARNING: <eikonal> node id " << index
                  << " is outside [0," << ndofs << "); ignored." << std::endl;
        continue;
      }
      // The per-node LAT is given in MILLISECONDS.
      lat(index) = node.attribute("lat").as_double();
      n_read++;
    }
  }

  if(has_eikonal)
  {
    // The per-node LAT from the mesh is used AS IT IS: the only operation is
    // the conversion from ms to the solver time unit. There is no shift and
    // no rescaling -- the activation sequence stored in the mesh is the
    // activation sequence the cells see. passive_time is only a fallback for
    // meshes without an <eikonal> section.
    lat *= ms_to_solver_time();

    std::cout << " Local activation time (per node, from the mesh)" << std::endl;
    std::cout << " Loaded LATs: " << n_read << " of " << ndofs << " nodes" << std::endl;
    if(n_read != (int) ndofs)
      std::cout << " *** WARNING: " << (ndofs - n_read) << " node(s) without LAT"
                << " were left at 0 and will activate at the start of the beat."
                << std::endl;
    std::cout << " Earliest activation: " << lat.min()
              << "  Latest activation: " << lat.max()
              << "  (spread " << (lat.max() - lat.min())
              << ", solver time units)" << std::endl;
  }
  else
  {
    lat.fill(begin_active_stress); //if there isn't lat in the mesh file, we use the passive_time for all elements. 
    std::cout << " No per-node LAT in the mesh: uniform activation at "
              << begin_active_stress << " (solver time units, from passive_time)"
              << std::endl;
  }
}


void Eikonal::solve()
{
  //TODO: solve eikonal and set lat into the cellmodel
  //For now, this function is been used only for debuging
  cout << "\nSimulating" << endl;

  initial_conditions();
    
  // loop in time
  int step=0;
  //int step_apd=0;

  //stimuli.check(tip.time(), *mesh, stim_nodes, &stim_val, &stim_apply);
  cells->advance(tip.time(), timestep, stim_val, stim_nodes);
  //cells->get_var(0, v0);
  arma::vec ta(ndofs);

  while( !tip.finished() )
  {
    tip.increase_time();
    tip.show_time();
    
    timer.enter("ODEs");
    solve_odes();
    cells->get_monitored_values(0, ta);

    arma::uvec non_zero_indices = arma::find(ta != 0.0);
    int count_non_zero = non_zero_indices.n_elem;
    int count_zero = ta.n_elem - count_non_zero;

    timer.leave();
  }  

  timer.summary();
}

void Eikonal::solve_odes()
{
  // todo:  simplify this function
  stimuli.check(tip.time(), *mesh, stim_nodes, &stim_val, &stim_apply);
  
  if (stim_apply)
  {
    cells->advance(tip.time(), timestep, stim_val, stim_nodes);
    stim_nodes.clear();
  }
  else if(stim_apply_nodes)
  {
    //cout << "Aplicando estimulos " << tip.time() << endl;
    cells->advance(tip.time(), timestep, stim_values);
    stim_values.fill(0);
    stim_apply_nodes = false;
  }
  else
  {
    cells->advance(tip.time(), timestep);
  }

  // NOTE: there used to be a second, unconditional cells->advance() here.
  // It made the cell model integrate TWO steps per time step, so the cell
  // clock ran at twice the simulation clock and activation happened at half
  // the intended time.
}