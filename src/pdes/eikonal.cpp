#include "eikonal.hpp"
#include "mesh/writer_hdf5.hpp"
#include <queue>
#include <iostream>
#include <fstream>
#include <cmath>

static int conductivity_from_string(const std::string & s)
{
  if (s == "S_ISOTROPIC")   return S_ISOTROPIC;
  if (s == "S_TRANSVERSE")  return S_TRANSVERSE;
  if (s == "S_ORTHOTROPIC") return S_ORTHOTROPIC;
  if (s == "M_ISOTROPIC")   return M_ISOTROPIC;
  if (s == "M_TRANSVERSE")  return M_TRANSVERSE;
  if (s == "M_ORTHOTROPIC") return M_ORTHOTROPIC;

  cout << " *** WARNING: unknown conductivity_type '" << s
       << "', using M_TRANSVERSE." << endl;
  return M_TRANSVERSE;
}

Eikonal::Eikonal() : condtype(M_TRANSVERSE), 
                     mesh(nullptr),
                     default_vf(1.0), 
                     default_vs(1.0), 
                     default_vn(1.0) {}

Eikonal::~Eikonal() { if (owns_mesh) delete mesh; }                

void Eikonal::set_velocities(double vf, double vs, double vn) {
    default_vf = vf;
    default_vs = vs;
    default_vn = vn;
}

void Eikonal::setup(const toml::table & cfg)
{
    auto meshcfg = cfg["problem"]["mesh"].value<std::string>();
    if (!meshcfg)
        throw std::runtime_error("Missing mandatory config: problem.mesh");

    mesh = new Mesh();
    mesh->read_xml(*meshcfg);

    double vf = cfg["physical"]["vel_f"].value_or(default_vf);
    double vs = cfg["physical"]["vel_s"].value_or(default_vs);
    double vn = cfg["physical"]["vel_n"].value_or(default_vn);
    set_velocities(vf, vs, vn);

    std::cout << " Conduction velocities [f, s, n]: "
              << vf << ", " << vs << ", " << vn << std::endl;

    if (auto c = cfg["physical"]["conductivity_type"].value<std::string>())
        set_conductivity(conductivity_from_string(*c));

    root_nodes.clear();
    root_times.clear();
    if (auto arr = cfg["activation"]["root"].as_array()) {
        for (auto & e : *arr) {
            if (auto t = e.as_table()) {
                root_nodes.push_back((*t)["id"].value_or(-1));
                root_times.push_back((*t)["time"].value_or(0.0));
            }
        }
    }

    output_filename = cfg["output"]["filename"].value_or(std::string("eikonal_output"));
}

void Eikonal::solve() {

    std::cout << "Computing activation time via Eikonal Solver" << std::endl;
    uint ndofs = mesh->get_n_points();
    if (mesh == nullptr) {
        std::cerr << "ERROR: setup() was not called before solve()." << std::endl;
        return;
    }
    if (root_nodes.empty()) {
        std::cout << "No root nodes provided; skipping Eikonal solve." << std::endl;
        return;
    }

    lat.zeros(ndofs);

    arma::mat33 g(arma::fill::zeros);
    g(0,0) = default_vf*default_vf;
    g(1,1) = default_vs*default_vs;
    g(2,2) = default_vn*default_vn;

    std::vector<std::map<int, double>> edge_costs(ndofs);
    int num_elements = mesh->get_n_elements();

    for (int i = 0; i < num_elements; ++i) {
        const Element& el = mesh->get_element(i);
        std::vector<int> pnums;
        mesh->get_element_pt_nums(i, pnums);

        arma::mat33 fsn;
        fsn.col(0) = el.get_fiber();
        fsn.col(1) = el.get_trans();
        fsn.col(2) = el.get_normal();

        arma::mat33 aux2 = fsn * g * fsn.t(); 
        arma::mat33 aux2_inv;
        bool is_invertible = arma::inv(aux2_inv, aux2); 

        for (size_t a = 0; a < pnums.size(); ++a) {
            for (size_t b = a + 1; b < pnums.size(); ++b) {
                int u = pnums[a], v = pnums[b];
                arma::vec3 edge_vec = mesh->get_point(v) - mesh->get_point(u);
                double cost = 0.0;

                if (is_invertible) {
                    double val = arma::dot(aux2_inv * edge_vec, edge_vec);
                    cost = (val > 0.0) ? std::sqrt(val) : 0.0;
                }

                if (edge_costs[u].find(v) == edge_costs[u].end() || cost < edge_costs[u][v]) {
                    edge_costs[u][v] = cost;
                    edge_costs[v][u] = cost;
                }
            }
        }
    }

    std::vector<std::vector<std::pair<int, double>>> adj_cost(ndofs);
    for (uint u = 0; u < ndofs; ++u) {
        for (auto const& edge : edge_costs[u]) {
            adj_cost[u].push_back({edge.first, edge.second});
        }
    }
    
    solve_dijkstra(ndofs, root_nodes, root_times, adj_cost);
    
    std::cout << " Earliest activation (min): " << lat.min() << std::endl;
    std::cout << " Latest activation (max): " << lat.max() << std::endl;

    WriterHDF5 writer(mesh);
    writer.write_eikonal_lat(output_filename, lat.memptr());

    std::ofstream arquivo(output_filename + ".txt");
    for (uint u = 0; u < ndofs; u++) {
        arquivo << "<node id=\"" << u << "\" lat=\"" << lat[u] << "\" />\n";
    }
    
    std::cout << "Output saved successfully." << std::endl;
}

void Eikonal::solve_dijkstra(int ndofs, const std::vector<int>& root_nodes, const std::vector<double>& root_times, const std::vector<std::vector<std::pair<int, double>>>& adj_cost) {
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

void Eikonal::set_conductivity(int cond) {
    CondTensorType tcond = static_cast<CondTensorType>(cond);
    condtype = tcond;
  
    std::cout << "Conductivity type: ";
    switch(condtype) {
        case S_ISOTROPIC:   std::cout << "S_ISOTROPIC"   << std::endl; break;
        case S_TRANSVERSE:  std::cout << "S_TRANSVERSE"  << std::endl; break;
        case S_ORTHOTROPIC: std::cout << "S_ORTHOTROPIC" << std::endl; break;
        case M_ISOTROPIC:   std::cout << "M_ISOTROPIC"   << std::endl; break;
        case M_TRANSVERSE:  std::cout << "M_TRANSVERSE"  << std::endl; break;
        case M_ORTHOTROPIC: std::cout << "M_ORTHOTROPIC" << std::endl; break;
    }
}