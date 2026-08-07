#ifndef GRAPH_COLORING_HPP
#define GRAPH_COLORING_HPP

#include <iostream>
#include <armadillo>
#include <armadillo>
#include "mesh.hpp"

using namespace std;

// Forward declaration
class Mesh;

arma::umat greedy_coloring(const arma::umat & L, const int threads);

arma::umat comm_matrix(const Mesh & mesh);


#endif
