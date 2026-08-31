// TO-DO: corrigir lance do ODE_solver



#include "eikonal.hpp"
#include "monodomain.hpp"
#include <queue>
#include <vector>
#include <utility>
#include "mesh/writer_hdf5.hpp" 

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

  parameters.rename("Eikonal_parameters");
  parameters.add("vel_f", 0.006);
  parameters.add("vel_s", 0.0002);
  parameters.add("vel_n", 0.0002);

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


void Eikonal::setup(std::string & b, std::string & c, std::string & m,
                           double dt, double T, double pr, double pa)
{
  timestep  = dt;
  totaltime = T;
  printrate = pr;
  printrate_apd = pa;
  
  mesh_filename = b;
  stimuli_filename = b;

  cell_name = c;
  odesolver = m;
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
  std::cout << " -- Initializing local activation time --" << std::endl; 
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

  bool loaded_lat = false;
  bool has_root_nodes = false;
  int n_read = 0;
  
  std::vector<int> root_nodes;
  std::vector<double> root_times;

  if(eikonal_data)
  {
    // Verify if there is already lat for all nodes in the mesh
    if (eikonal_data.child("node")) 
    {
      loaded_lat = true;
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
      if(n_read != (int) ndofs)
      {
        std::cout << " *** WARNING: prescribed LAT incomplete (" << n_read << " of " << ndofs << " nodes)." << std::endl;
        std::cout << "     Falling back to Eikonal solving or uniform activation." << std::endl;
        loaded_lat = false; 
      }
    }
    // Verify if there are root nodes for solving eikonal
    if (eikonal_data.child("root_node") && !loaded_lat) 
    {
      has_root_nodes = true;
      for(pugi::xml_node node = eikonal_data.child("root_node"); node; node = node.next_sibling("root_node"))
      {
        root_nodes.push_back(node.attribute("id").as_int());
        root_times.push_back(node.attribute("time").as_double());
      }
    }
  }

  //If the lat is already precribed, don't solve eikonal just load the local activation time
  if(loaded_lat)
  {
    // The per-node LAT from the mesh is used AS IT IS: the only operation is
    // the conversion from ms to the solver time unit. There is no shift and
    // no rescaling -- the activation sequence stored in the mesh is the
    // activation sequence the cells see. passive_time is only a fallback for
    // meshes without an <eikonal> section.
    lat *= ms_to_solver_time();

    std::cout << " Local activation time (per node, from the mesh)" << std::endl;
    std::cout << " Loaded LATs: " << n_read << " of " << ndofs << " nodes" << std::endl;
    std::cout << " Earliest activation: " << lat.min()
              << "  Latest activation: " << lat.max()
              << "  (spread " << (lat.max() - lat.min())
              << ", solver time units)" << std::endl;
  }
  else if (has_root_nodes)
  {
    std::cout << " -- Computing local activation time via Eikonal Solver --" << std::endl;
    
    double vf = parameters["vel_f"]; 
    double vs = parameters["vel_s"];
    double vn = parameters["vel_n"];
    
    if (eikonal_data.attribute("vel_f")) vf = eikonal_data.attribute("vel_f").as_double();
    if (eikonal_data.attribute("vel_s")) vs = eikonal_data.attribute("vel_s").as_double();
    if (eikonal_data.attribute("vel_n")) vn = eikonal_data.attribute("vel_n").as_double();

    // conductivity tensor squared
    arma::mat33 g(arma::fill::zeros);
    g(0,0) = vf * vf;
    g(1,1) = vs * vs;
    g(2,2) = vn * vn;

    std::vector<std::map<int, double>> edge_costs(ndofs);

    int num_elements = mesh->get_n_elements();
    for (int i = 0; i < num_elements; ++i) 
    {
        const Element& el = mesh->get_element(i);
        
        std::vector<int> pnums;
        mesh->get_element_pt_nums(i, pnums);

        arma::vec3 f = el.get_fiber();
        arma::vec3 s = el.get_trans();
        arma::vec3 n = el.get_normal();

        arma::mat33 fsn;
        fsn.col(0) = f;
        fsn.col(1) = s;
        fsn.col(2) = n;

        // conductivity tensor projected in the global system
        arma::mat33 aux2 = fsn * g * fsn.t(); 
        
        // Em vez de resolver o sistema linear para cada aresta repetidamente, 
        // inverte-se a matriz 3x3 uma única vez por elemento para otimização de CPU.
        arma::mat33 aux2_inv;
        bool is_invertible = arma::inv(aux2_inv, aux2); 

        // Creating edges by permutating all nodes
        for (size_t a = 0; a < pnums.size(); ++a) 
        {
            for (size_t b = a + 1; b < pnums.size(); ++b) 
            {
                int u = pnums[a];
                int v = pnums[b];

                arma::vec3 edge_vec = mesh->get_point(v) - mesh->get_point(u);
                double cost = 0.0;

                if (is_invertible) 
                {
                    // aux2*x = edge_vec is the same as x = aux2_inv * edge_vec
                    arma::vec3 x = aux2_inv * edge_vec;
                    double val = arma::dot(x, edge_vec);
                    cost = (val > 0.0) ? std::sqrt(val) : 0.0;
                }

                if (edge_costs[u].find(v) == edge_costs[u].end() || cost < edge_costs[u][v]) 
                {
                    edge_costs[u][v] = cost;
                    edge_costs[v][u] = cost;
                }
            }
        }
    }

    // adjacency list by node
    std::vector<std::vector<std::pair<int, double>>> adj_cost(ndofs);
    for (uint u = 0; u < ndofs; ++u) 
    {
        for (auto const& edge : edge_costs[u]) 
        {
            adj_cost[u].push_back({edge.first, edge.second});
        }
    }
    
    solve_dijkstra(root_nodes, root_times, adj_cost);
    
    std::cout << " Computed Earliest activation: " << lat.min() << "  Latest activation: " << lat.max() << std::endl;

    WriterHDF5 writer(mesh);
    writer.write_eikonal_lat(mshfile, lat.memptr());

    std::ofstream arquivo("eikonal.txt");
    for (uint u = 0; u<ndofs; u++)
    {
      arquivo << "<node id=\"" <<u<<"\" lat=\"" <<lat[u]<<"\" />" <<endl; 
    }
    std::cout << " -- LAT saved successfully to HDF5/XDMF format. --" << std::endl;
  }
  else
  {
    lat.fill(begin_active_stress); // if there isn't lat in the mesh file, we use the passive_time for all elements. 
    std::cout << " No valid/complete per-node LAT or root nodes found in the mesh." << std::endl;
    std::cout << " Uniform activation at "
              << begin_active_stress << " (solver time units, from passive_time)"
              << std::endl;
  }
}


void Eikonal::solve()
{

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
}

void Eikonal::solve_dijkstra(const std::vector<int>& root_nodes, 
                             const std::vector<double>& root_times,
                             const std::vector<std::vector<std::pair<int, double>>>& adj_cost) 
{
    lat.set_size(ndofs);
    lat.fill(0.0);
    
    std::vector<bool> visited(ndofs, false);
    
    std::vector<double> temp_times(ndofs, std::numeric_limits<double>::infinity());
    
    std::priority_queue<EikonalNode, std::vector<EikonalNode>, std::greater<EikonalNode>> min_heap;

    for (size_t i = 0; i < root_nodes.size(); ++i) {
        int root = root_nodes[i];
        double time = root_times[i];
        
        temp_times[root] = time;
        min_heap.push({time, root});
    }

    while (!min_heap.empty()) {
        EikonalNode current = min_heap.top();
        min_heap.pop();

        int u = current.id;
        double current_cost = current.cost;

        if (visited[u]) continue;

        visited[u] = true;
        lat(u) = current_cost; 

        for (const auto& edge : adj_cost[u]) {
            int v = edge.first;
            double edge_weight = edge.second;

            if (!visited[v]) {
                double new_cost = current_cost + edge_weight;
                if (new_cost < temp_times[v]) {
                    temp_times[v] = new_cost;
                    min_heap.push({new_cost, v});
                }
            }
        }
    }
}