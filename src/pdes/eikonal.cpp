#include "eikonal.hpp"
#include "mesh/writer_hdf5.hpp"
#include "util/pugixml.hpp"
#include <queue>
#include <iostream>
#include <fstream>
#include <cmath>

Eikonal::Eikonal() : default_vf(0.006), default_vs(0.0002), default_vn(0.0002) {}

void Eikonal::set_velocities(double vf, double vs, double vn) {
    default_vf = vf;
    default_vs = vs;
    default_vn = vn;
}

void Eikonal::solve(Mesh* mesh, const std::string &mshfile) {
    std::cout << " -- Computing local activation time via Eikonal Solver --" << std::endl;
    uint ndofs = mesh->get_n_points();
    pugi::xml_document doc;
    doc.load_file(mshfile.c_str());

    lat.zeros(ndofs);
    pugi::xml_node eikonal_data = doc.child("mesh").child("eikonal");

    std::vector<int> root_nodes;
    std::vector<double> root_times;

    if(eikonal_data && eikonal_data.child("root_node")) {
        for(pugi::xml_node node = eikonal_data.child("root_node"); node; node = node.next_sibling("root_node")) {
            root_nodes.push_back(node.attribute("id").as_int());
            root_times.push_back(node.attribute("time").as_double());
        }
    }

    if (root_nodes.empty()) {
        std::cerr << "ERROR: invalid or missing root nodes information for eikonal" << std::endl;
        exit(1);
    }

    double vf = default_vf, vs = default_vs, vn = default_vn;
    if (eikonal_data.attribute("vel_f")) vf = eikonal_data.attribute("vel_f").as_double();
    if (eikonal_data.attribute("vel_s")) vs = eikonal_data.attribute("vel_s").as_double();
    if (eikonal_data.attribute("vel_n")) vn = eikonal_data.attribute("vel_n").as_double();

    arma::mat33 g(arma::fill::zeros);
    g(0,0) = vf * vf; g(1,1) = vs * vs; g(2,2) = vn * vn;

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
    
    std::cout << " Computed Earliest activation: " << lat.min() << "  Latest activation: " << lat.max() << std::endl;

    WriterHDF5 writer(mesh);
    writer.write_eikonal_lat(mshfile, lat.memptr());

    std::ofstream arquivo("eikonal.txt");
    for (uint u = 0; u<ndofs; u++) {
        arquivo << "<node id=\"" << u << "\" lat=\"" << lat[u] << "\" />\n"; 
    }
    std::cout << " -- LAT saved successfully. --" << std::endl;
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