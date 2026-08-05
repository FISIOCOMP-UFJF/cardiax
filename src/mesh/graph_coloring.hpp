#ifndef GRAPH_COLORING_HPP
#define GRAPH_COLORING_HPP

#include <vector>
#include "mesh.hpp"


class Mesh;
/**
 * Computa a coloração gulosa (greedy coloring) para os elementos da malha.
 * 
 * @param mesh Objeto Mesh contendo a conectividade.
 * @return Um std::vector onde o índice externo é a ID da cor (0, 1, 2...), 
 *         e o std::vector interno contém as IDs dos elementos daquela cor.
 */
std::vector<std::vector<int>> compute_element_colors(const Mesh & mesh);

#endif