#ifndef CELLS_GPU_HPP
#define CELLS_GPU_HPP

#include <vector>
#include <set>                       // <-- ADICIONAR
#include <armadillo>                
#include "linalg/petsc_vector.hpp"   // <- ADICIONAR


typedef unsigned int uint;           // <-- ADICIONAR

// Forward declaration do wrapper CUDA
void launch_euler_kernel(double* d_states, const double* d_istim, const int* d_celltypes, const double dt, const int num_systems);

class CellsGpu {
private:
    int num_systems;
    int num_states;

    // Ponteiros no Device (VRAM)
    double* d_states;
    double* d_istim;
    int* d_celltypes;

    // Ponteiros no Host (RAM) para I/O
    double* h_states;

public:
    CellsGpu(int n_systems);
    ~CellsGpu();

    void init_default_conditions(const std::vector<int>& host_celltypes);
    void set_stimuli(const std::vector<double>& host_istim);
    void retrieve_states();
    
    // Em cells_gpu.hpp
    void send_to_device_vector(int vindex, double* dest_d_ptr) const;
    void retrieve_from_device_vector(int vindex, const double* src_d_ptr);
    
    // <- ADICIONAR: Função para enviar os estados da RAM de volta pra VRAM (usado no Checkpoint)
    void push_states(); 

    int get_num_state_vars() const { return num_states; }
    double* get_state_vars_ptr();

    void get_var_from_device(int var_idx, double* host_array) const;
    void get_var(int vindex, arma::vec &v) const;
    void get_var(int vindex, petsc::Vector &v) const;
    void set_var(int vindex, const arma::vec &v) const;

    
    // NOVO: Métodos antigos de leitura pontual
    double get_state(int cell_idx, int var_idx) const;
    void get_monitored_values(int mindex, arma::vec &v) const;

    // NOVO: Sobrecargas do advance imitando a API antiga
    void advance(double dt);
    void advance(double t, double dt, double stim_val, const std::set<uint>& stim_nodes);
    void advance(double t, double dt, const arma::vec& stim_values);
    void advance(double t, double dt);

    double* get_var_device_ptr(int vindex) const;

};

#endif