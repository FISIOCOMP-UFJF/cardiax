#ifndef EIKONAL_HPP
#define EIKONAL_HPP

#include <vector>
#include <string>
#include <utility>
#include <map>
#include <armadillo>
#include "mesh/mesh.hpp"

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

    void solve_dijkstra(int ndofs,
                        const std::vector<int>& root_nodes, 
                        const std::vector<double>& root_times,
                        const std::vector<std::vector<std::pair<int, double>>>& adj_cost);

public:
    Eikonal();
    ~Eikonal() = default;

    void set_velocities(double vf, double vs, double vn);

    void solve(Mesh* mesh, const std::string &mshfile);

    const arma::vec& get_lat() const { return lat; }

    void set_conductivity(int cond);
};

#endif