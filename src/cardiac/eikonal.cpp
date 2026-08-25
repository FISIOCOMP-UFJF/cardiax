// TO-DO: corrigir lance do ODE_solver



#include "eikonal.hpp"
#include "monodomain.hpp"
#include <queue>
#include <vector>
#include <utility>
#include "mesh/writer_hdf5.hpp" 

Eikonal::Eikonal()
  : CardiacProblem(),
    //sigma_l(0.0001334), sigma_t(0.0000176),  sigma_n(0.0000176),
    stim_apply_nodes(false)
{
  cout << "Eikonal" << endl; 
  mesh = new Mesh();

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

  // setup data writer to write at every 1 ms
  // potential field and displacements
  std::size_t pos  = mesh_filename.find(".xml");
  int nsteps = tip.get_size(); 

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
  std::cout << " -- Initializing local activation time --" << std::endl; 
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(mshfile.c_str());
  
  lat.set_size(ndofs);
  
  pugi::xml_node eikonal_data = doc.child("mesh").child("eikonal");
  
  bool loaded_lat = false;
  bool has_root_nodes = false;
  
  std::vector<int> root_nodes;
  std::vector<double> root_times;

  if(eikonal_data)
  {
    // Verify if there is already lat for all nodes in the mesh
    if (eikonal_data.child("node")) 
    {
      int lat_node_values = 0; 
      loaded_lat = true;
      for(pugi::xml_node node = eikonal_data.child("node"); node; node = node.next_sibling("node"))
      {
        int index = node.attribute("id").as_int(); 
        lat(index) = std::stod(node.attribute("lat").as_string()); 
        lat_node_values++;
      }
      if(lat_node_values != ndofs)
      {
        cout<<"Warning: precribed LAT incomplete (falling back to eikonal solving)" <<endl; 
        loaded_lat = false; 
      }
    }
    // Verify if there is root nodes for solving eikonal
    if (eikonal_data.child("root_node")  && !loaded_lat) 
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
    double min_val = lat.min(); 
    double max_val = lat.max(); 

    std::cout << " -- Reading local activation time from mesh file --" << std::endl;

    pugi::xml_node pvloop_data = doc.child("mesh").child("pvloop");
    double begin_active_stress = 0.0; 
    
    if(pvloop_data) begin_active_stress = std::stod(pvloop_data.attribute("passive_time").as_string()) / 1000.0; 
    double latest_lat = begin_active_stress + 0.146; 

    if(max_val - min_val == 0.0) 
      lat.fill(0.0);
    else 
      lat = begin_active_stress + (lat - min_val) * (latest_lat - begin_active_stress) / (max_val - min_val);
      
    std::cout << " Earliest activation: " << lat.min() << "  Latest activation: " << lat.max()  << std::endl; 
    std::cout << " Loaded LATs: " << lat.n_elem << " values.\n";
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
    std::cout << " -- No LAT or root nodes found. Using passive_time for all nodes --" << std::endl;
    
    pugi::xml_node pvloop_data = doc.child("mesh").child("pvloop");
    double begin_active_stress = 0.0; 
    
    if(pvloop_data) begin_active_stress = std::stod(pvloop_data.attribute("passive_time").as_string()); 
    else std::cout<<" -- No passive_time, using lat = 0.0 for all nodes -- " <<endl; 
    lat.fill(begin_active_stress); 
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

  cells->advance(tip.time(), timestep);
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