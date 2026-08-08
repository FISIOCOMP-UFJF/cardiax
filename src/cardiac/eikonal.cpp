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
  // writer = new WriterHDF5(mesh);

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

  // setup data writer to write at every 1 ms
  // potential field and displacements
  std::size_t pos  = mesh_filename.find(".xml");
  // std::string output = mesh_filename.substr(0,pos) + "_output";
  int nsteps = tip.get_size(); 
  // writer->open(output, nsteps+1, timestep);

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
    // Prioridade 1: Verifica se existem tags <node> especificando o LAT mapeado[cite: 2]
    if (eikonal_data.child("node")) 
    {
      loaded_lat = true;
      for(pugi::xml_node node = eikonal_data.child("node"); node; node = node.next_sibling("node"))
      {
        int index = node.attribute("id").as_int(); 
        lat(index) = std::stod(node.attribute("lat").as_string()); 
      }
    }
    // Prioridade 2: Caso não existam <node>, verifica tags <root_node> para invocar o Eikonal
    else if (eikonal_data.child("root_node")) 
    {
      has_root_nodes = true;
      for(pugi::xml_node node = eikonal_data.child("root_node"); node; node = node.next_sibling("root_node"))
      {
        root_nodes.push_back(node.attribute("id").as_int());
        root_times.push_back(node.attribute("time").as_double());
      }
    }
  }

  // Bloco de Execução baseado nas marcações lidas
  if(loaded_lat)
  {
    // Aplica o reescalonamento/normalização clássico caso os dados venham do arquivo[cite: 2]
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
    
    // 1. Calcule as penalidades (custos) de navegação para a malha atual
    // std::vector<std::vector<std::pair<int, double>>> adj_cost = compute_navigation_costs();
    
    // 2. Chame a função Dijkstra recém-criada
    // solve_dijkstra(root_nodes, root_times, adj_cost);
    
    // (Opcional) Aplique a mesma normalização matemática do bloco acima ao 'lat' recém computado,
    // garantindo que os tempos da simulação fiquem limitados a 'latest_lat'.
  }
  else
  {
    // Prioridade 3: Fallback padrão[cite: 2]
    std::cout << " -- No LAT or root nodes found. Using passive_time for all elements --" << std::endl;
    
    pugi::xml_node pvloop_data = doc.child("mesh").child("pvloop");
    double begin_active_stress = 0.0; 
    
    if(pvloop_data) begin_active_stress = std::stod(pvloop_data.attribute("passive_time").as_string()); 
    
    lat.fill(begin_active_stress); 
  }
}
void Eikonal::solve()
{
  //TODO: solve eikonal and set lat into the cellmodel
  //For now, this function is been used only for debuging
  
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
    // Inicializa os vetores de estado
    lat.set_size(ndofs);
    lat.fill(0.0);
    
    std::vector<bool> visited(ndofs, false);
    
    // Substitui o valor mágico 1e6 por infinity real do limite numérico
    std::vector<double> temp_times(ndofs, std::numeric_limits<double>::infinity());
    
    // Declaração do min-heap
    std::priority_queue<EikonalNode, std::vector<EikonalNode>, std::greater<EikonalNode>> min_heap;

    // Inicialização dos root-nodes
    for (size_t i = 0; i < root_nodes.size(); ++i) {
        int root = root_nodes[i];
        double time = root_times[i];
        
        temp_times[root] = time;
        min_heap.push({time, root});
    }

    // Processamento dos caminhos mínimos
    while (!min_heap.empty()) {
        EikonalNode current = min_heap.top();
        min_heap.pop();

        int u = current.id;
        double current_cost = current.cost;

        // Pula nós que receberam atualizações mais rápidas após a inserção na fila
        if (visited[u]) continue;

        visited[u] = true;
        lat(u) = current_cost; // Consolida o tempo de ativação (LAT)

        // Expansão geométrica (Relaxamento das arestas)
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