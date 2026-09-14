#ifndef EIKONAL_HPP
#define EIKONAL_HPP

#include <vector>
#include <string>
#include <utility>
#include <map>
#include <armadillo>
#include "mesh/mesh.hpp"
#include "util/toml.hpp"

// Struct auxiliar para a fila de prioridade do Dijkstra
struct EikonalNode {
    double cost;
    int id;
    bool operator>(const EikonalNode& other) const {
        return cost > other.cost;
    }
};

enum CondTensorType {
    S_ISOTROPIC, S_TRANSVERSE, S_ORTHOTROPIC,
    M_ISOTROPIC, M_TRANSVERSE, M_ORTHOTROPIC
};

class Eikonal {
private:
    CondTensorType condtype;

    arma::vec lat;

    double default_vf;
    double default_vs;
    double default_vn;

    Mesh* mesh;
    std::vector<int>    root_nodes;
    std::vector<double> root_times;

    std::string output_filename = "eikonal_output";

    void solve_dijkstra(int ndofs,
                        const std::vector<int>& root_nodes,
                        const std::vector<double>& root_times,
                        const std::vector<std::vector<std::pair<int, double>>>& adj_cost);

    bool owns_mesh = true;                    

public:
    Eikonal();
    ~Eikonal();

    const arma::vec& get_lat() const { return lat; }

    void set_conductivity(int cond);

    void set_mesh(Mesh* m) { mesh = m; owns_mesh = false; }

    void set_velocities(double vf, double vs, double vn);

    void set_root_nodes(const std::vector<int>& nodes,
                        const std::vector<double>& times) {
        root_nodes = nodes;
        root_times = times;
    }
    const std::vector<int>&    get_root_nodes() const { return root_nodes; }
    const std::vector<double>& get_root_times() const { return root_times; }

    void setup(const toml::table & cfg);

    void solve();
};

#endif