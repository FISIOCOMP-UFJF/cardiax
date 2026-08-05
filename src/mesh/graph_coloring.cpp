#include "graph_coloring.hpp"
#include <iostream>
#include <algorithm>
#include "mesh.hpp"


std::vector<std::vector<int>> compute_element_colors(const Mesh & mesh)
{
    // Recupera informações básicas da malha usando os métodos já existentes[cite: 5]
    int num_nodes = mesh.get_n_points();
    int num_elements = mesh.get_n_elements();

    // -------------------------------------------------------------------------
    // PASSO 1: Mapeamento Nó -> Elemento (Lista de Adjacência Invertida)
    // Evita a varredura O(N^2) na busca de vizinhos[cite: 5].
    // -------------------------------------------------------------------------
    std::vector<std::vector<int>> node_to_elems(num_nodes);
    
    for (int e = 0; e < num_elements; ++e) 
    {
        std::vector<int> edofs;
        mesh.get_element_pt_nums(e, edofs); //[cite: 5]
        
        for (int node : edofs) 
        {
            node_to_elems[node].push_back(e);
        }
    }

    // -------------------------------------------------------------------------
    // PASSO 2: Coloração Gulosa (Greedy Coloring)
    // -------------------------------------------------------------------------
    std::vector<int> element_color(num_elements, -1);
    int max_color_used = -1;

    for (int e = 0; e < num_elements; ++e) 
    {
        // 2.a Encontrar os elementos vizinhos (que compartilham nós)
        std::vector<int> neighbors;
        std::vector<int> edofs;
        mesh.get_element_pt_nums(e, edofs); //[cite: 5]
        
        for (int node : edofs) 
        {
            for (int neighbor_elem : node_to_elems[node]) 
            {
                if (neighbor_elem != e) 
                {
                    neighbors.push_back(neighbor_elem);
                }
            }
        }

        // 2.b Mapear as cores que já estão sendo usadas pelos vizinhos
        // Alocamos o array de disponibilidade com tamanho 'max_color_used + 2'
        // para garantir espaço caso precisemos de uma cor nova.
        std::vector<bool> color_available(max_color_used + 2, true);
        
        for (int neighbor_elem : neighbors) 
        {
            int c = element_color[neighbor_elem];
            if (c != -1) 
            {
                color_available[c] = false; // Cor bloqueada
            }
        }

        // 2.c Encontrar a menor cor livre
        int chosen_color = 0;
        while (!color_available[chosen_color]) 
        {
            chosen_color++;
        }

        // 2.d Atribuir a cor e atualizar o máximo
        element_color[e] = chosen_color;
        if (chosen_color > max_color_used) 
        {
            max_color_used = chosen_color;
        }
    }

    // -------------------------------------------------------------------------
    // PASSO 3: Agrupar elementos por cor
    // -------------------------------------------------------------------------
    int total_colors = max_color_used + 1;
    std::vector<std::vector<int>> color_groups(total_colors);

    // Pré-alocar memória para evitar realocações durante o push_back (Otimização)
    // Estimativa: elementos divididos uniformemente pelas cores
    int avg_elements_per_color = num_elements / total_colors;
    for (int c = 0; c < total_colors; ++c) {
        color_groups[c].reserve(avg_elements_per_color * 1.5); 
    }

    // Distribuir os elementos
    for (int e = 0; e < num_elements; ++e) 
    {
        color_groups[element_color[e]].push_back(e);
    }

    std::cout << "Malha agrupada utilizando " << total_colors << " cores independentes.\n";

    return color_groups;
}