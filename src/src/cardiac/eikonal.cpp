#include "eikonal.hpp"
#include "monodomain.hpp"

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
  ndofs = mesh->get_n_elements(); //Before: get_n_points()

  // setup data writer to write at every 1 ms
  // potential field and displacements
  std::size_t pos  = mesh_filename.find(".xml");
  std::string output = mesh_filename.substr(0,pos) + "_output";
  int nsteps = tip.get_size(); //tip.get_nsteps();// /(1.0/timestep); //troquei para conseguir rodar casos maiores (juvs)
  writer->open(output, nsteps+1, timestep);

  // setup model and cells
  cellmodel = CellModel::create(cell_name);
  cellmodel->setup(odesolver, timestep, totaltime, 1.0);
  cells = new Cells(ndofs,cellmodel);

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
  cells->set_var(1, lat); //TODO: Será que eu devo fazer isso?

  // loop in time
  int step=0;

  cells->advance(tip.time(), timestep, stim_val, stim_nodes);
  
}

void Eikonal::set_stimulus_value(int index, double val)
{
  stim_apply_nodes = true;
  stim_values(index) = val;
}

void Eikonal::solve(const string &mshfile)
{
  cout << " -- Reading local activation time from mesh file --" << endl; 
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());
  
  lat.set_size(ndofs);
  pugi::xml_node element_data = doc.child("mesh").child("element_data");
  bool has_eikonal = false; 
  if(element_data)
  {
    for (pugi::xml_node elem = element_data.child("element"); elem; elem = elem.next_sibling("element"))
    {
        pugi::xml_node eikonal_p_elem = elem.child("eikonal");
        if(eikonal_p_elem)
        {
          has_eikonal = true; 
          int index; 
          index = elem.attribute("id").as_int();
          lat(index) = std::stod(eikonal_p_elem.text().as_string());
        }
        else
        {
          if(has_eikonal)
          {
            std::cerr << "ERROR: Some elements have a local activation time, and others do not.\n";
            exit(6);
          }
        }
    }
  }

  if(has_eikonal)
  {
    double min_val = lat.min(); 
    double max_val = lat.max(); 

    std::cout<<"Local activation time"<<std::endl;

    //TODO: read this information from file? 
    double begin_active_stress = 0.136;
    double latest_lat = 0.136 + 0.146;

    lat = begin_active_stress + (lat - min_val) * (latest_lat - begin_active_stress)/ (max_val-min_val);
    std::cout << " Earliest activation: " << lat.min() << "  Latest activation: " << lat.max()  << std::endl; 
    std::cout << " Loaded LATs: " << lat.n_elem << " values.\n";
  }
  else
  {
    lat.fill(0.136); //if there isn't lat in the mesh file, we use 0.136 for all elements.
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

    cout << "Zeros exatos: " << count_zero << endl;
    cout << "Diferentes de zero: " << count_non_zero << endl;

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

  cells->advance(tip.time(), timestep);
}



