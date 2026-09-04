#include "cells_gpu.hpp"
#include <cuda_runtime.h>
#include <iostream>

CellsGpu::CellsGpu(int n_systems) : num_systems(n_systems), num_states(50) {
    size_t total_bytes_states = num_systems * num_states * sizeof(double);
    size_t total_bytes_stim = num_systems * sizeof(double);
    size_t total_bytes_types = num_systems * sizeof(int);

    // Alocação no Host
    h_states = new double[num_systems * num_states];

    // Alocação no Device
    cudaMalloc((void**)&d_states, total_bytes_states);
    cudaMalloc((void**)&d_istim, total_bytes_stim);
    cudaMalloc((void**)&d_celltypes, total_bytes_types);

    // Zera os estímulos inicialmente na GPU
    cudaMemset(d_istim, 0, total_bytes_stim);
}

CellsGpu::~CellsGpu() {
    delete[] h_states;
    cudaFree(d_states);
    cudaFree(d_istim);
    cudaFree(d_celltypes);
}

void CellsGpu::init_default_conditions(const std::vector<int>& host_celltypes) {
    // Configura os estados iniciais no Host respeitando o layout SoA
    for (int i = 0; i < num_systems; i++) {
        h_states[0* num_systems + i] = -8.863699e+01;
        h_states[1* num_systems + i] = 1.189734e+01;
        h_states[2* num_systems + i] = 1.189766e+01;
        h_states[3* num_systems + i] = 1.412345e+02;
        h_states[4* num_systems + i] = 1.412344e+02;
        h_states[5* num_systems + i] = 7.267473e-05;
        h_states[6* num_systems + i] = 6.337870e-05;
        h_states[7* num_systems + i] = 1.532653e+00;
        h_states[8* num_systems + i] = 1.533946e+00;
        h_states[9* num_systems + i] = 8.280078e-04;
        h_states[10* num_systems + i] = 6.665272e-01;
        h_states[11* num_systems + i] = 8.260208e-01;
        h_states[12* num_systems + i] = 8.260560e-01;
        h_states[13* num_systems + i] = 8.258509e-01;
        h_states[14* num_systems + i] = 1.668686e-04;
        h_states[15* num_systems + i] = 5.228306e-01;
        h_states[16* num_systems + i] = 2.859696e-01;
        h_states[17* num_systems + i] = 9.591370e-04;
        h_states[18* num_systems + i] = 9.996012e-01;
        h_states[19* num_systems + i] = 5.934016e-01;
        h_states[20* num_systems + i] = 4.886961e-04;
        h_states[21* num_systems + i] = 9.996011e-01;
        h_states[22* num_systems + i] = 6.546687e-01;
        h_states[23* num_systems + i] = 9.500075e-32;
        h_states[24* num_systems + i] = 1.000000e+00;
        h_states[25* num_systems + i] = 9.392580e-01;
        h_states[26* num_systems + i] = 1.000000e+00;
        h_states[27* num_systems + i] = 9.998984e-01;
        h_states[28* num_systems + i] = 9.999783e-01;
        h_states[29* num_systems + i] = 4.448162e-04;
        h_states[30* num_systems + i] = 7.550725e-04;
        h_states[31* num_systems + i] = 1.000000e+00;
        h_states[32* num_systems + i] = 1.000000e+00;
        h_states[33* num_systems + i] = 2.424047e-01;
        h_states[34* num_systems + i] = 1.795377e-04;
        h_states[35* num_systems + i] = -6.883086e-25;
        h_states[36* num_systems + i] = 1.117498e-02;
        h_states[37* num_systems + i] = 9.980366e-01;
        h_states[38* num_systems + i] = 8.588018e-04;
        h_states[39* num_systems + i] = 7.097447e-04;
        h_states[40* num_systems + i] = 3.812617e-04;
        h_states[41* num_systems + i] = 1.357116e-05;
        h_states[42* num_systems + i] = 2.302525e-23;
        h_states[43* num_systems + i] = 1.561941e-04;
        h_states[44* num_systems + i] = 2.351289e-04;
        h_states[45* num_systems + i] = 8.077631e-03;
        h_states[46* num_systems + i] = 9.993734e-01;
        h_states[47* num_systems + i] = 0.000000e+00;
        h_states[48* num_systems + i] = 0.000000e+00;
        h_states[49* num_systems + i] = 0.000000e+00;
    }

    // Copia estados e tipos celulares para a GPU (HtoD uma única vez)
    cudaMemcpy(d_states, h_states, num_systems * num_states * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_celltypes, host_celltypes.data(), num_systems * sizeof(int), cudaMemcpyHostToDevice);
}

void CellsGpu::set_stimuli(const std::vector<double>& host_istim) {
    // Esse HtoD deve ocorrer apenas quando o estímulo mudar, não em todo passo de dt
    cudaMemcpy(d_istim, host_istim.data(), num_systems * sizeof(double), cudaMemcpyHostToDevice);
}

void CellsGpu::advance(double dt) {
    // A integração acontece puramente na placa. Nenhuma transferência de memória aqui.
    launch_euler_kernel(d_states, d_istim, d_celltypes, dt, num_systems);
}
void CellsGpu::get_var_from_device(int var_idx, double* host_array) const {
    double* d_var_ptr = d_states + (var_idx * num_systems);
    cudaMemcpy(host_array, d_var_ptr, num_systems * sizeof(double), cudaMemcpyDeviceToHost);
}

void CellsGpu::retrieve_states() {
    cudaMemcpy(h_states, d_states, num_systems * num_states * sizeof(double), cudaMemcpyDeviceToHost);
}

double* CellsGpu::get_state_vars_ptr() {
    // Puxa os dados da VRAM para a RAM antes de gravar no arquivo HDF5
    retrieve_states(); 
    return h_states;
}

// Função necessária para restaurar o checkpoint HDF5
void CellsGpu::push_states() {
    cudaMemcpy(d_states, h_states, num_systems * num_states * sizeof(double), cudaMemcpyHostToDevice);
}

void CellsGpu::get_var(int vindex, arma::vec &v) const {
    // Como a memória na GPU é SoA, a variável inteira é um bloco contíguo
    double* d_var_ptr = d_states + (vindex * num_systems);
    cudaMemcpy(v.memptr(), d_var_ptr, num_systems * sizeof(double), cudaMemcpyDeviceToHost);
}

void CellsGpu::get_var(int vindex, petsc::Vector &v) const {
    // PETSc requer que puxemos para um array nativo primeiro
    double* temp_array = new double[num_systems];
    double* d_var_ptr = d_states + (vindex * num_systems);
    
    cudaMemcpy(temp_array, d_var_ptr, num_systems * sizeof(double), cudaMemcpyDeviceToHost);
    
    for(int i = 0; i < num_systems; i++) {
        v.set(i, temp_array[i]);
    }
    delete[] temp_array;
}

void CellsGpu::set_var(int vindex, const arma::vec &v) const {
    double* d_var_ptr = d_states + (vindex * num_systems);
    cudaMemcpy(d_var_ptr, v.memptr(), num_systems * sizeof(double), cudaMemcpyHostToDevice);
}

// Emula o loop de EDOs com set de nós estimulados
void CellsGpu::advance(double t, double dt, double stim_val, const std::set<uint>& stim_nodes) {
    std::vector<double> host_istim(num_systems, 0.0);
    for (uint node : stim_nodes) {
        host_istim[node] = stim_val;
    }
    set_stimuli(host_istim); // HtoD
    advance(dt);             // Executa na GPU
}

// Emula o loop de EDOs com vetor denso de estímulos
void CellsGpu::advance(double t, double dt, const arma::vec& stim_values) {
    cudaMemcpy(d_istim, stim_values.memptr(), num_systems * sizeof(double), cudaMemcpyHostToDevice);
    advance(dt);
}

// Emula o loop de EDOs sem estímulo
void CellsGpu::advance(double t, double dt) {
    std::vector<double> host_istim(num_systems, 0.0);
    set_stimuli(host_istim);
    advance(dt);
}

void CellsGpu::get_monitored_values(int mindex, arma::vec &v) const {
    // No TorordLand original, a Tensão Ativa (Ta) era monitorada. 
    // Na nossa estrutura CUDA, o Ta é a variável de estado 49.
    get_var(49, v);
}

double CellsGpu::get_state(int cell_idx, int var_idx) const {
    double val;
    double* d_var_ptr = d_states + (var_idx * num_systems) + cell_idx;
    cudaMemcpy(&val, d_var_ptr, sizeof(double), cudaMemcpyDeviceToHost);
    return val;
}

// Em cells_gpu.cpp
void CellsGpu::send_to_device_vector(int vindex, double* dest_d_ptr) const {
    double* src_d_ptr = d_states + (vindex * num_systems);
    // Cópia D2D: Velocidade da ordem de Terabytes por segundo (ex: 360 GB/s na RTX 3060)
    cudaMemcpy(dest_d_ptr, src_d_ptr, num_systems * sizeof(double), cudaMemcpyDeviceToDevice);
}

void CellsGpu::retrieve_from_device_vector(int vindex, const double* src_d_ptr) {
    double* dest_d_ptr = d_states + (vindex * num_systems);
    cudaMemcpy(dest_d_ptr, src_d_ptr, num_systems * sizeof(double), cudaMemcpyDeviceToDevice);
}

double* CellsGpu::get_var_device_ptr(int vindex) const {
    // Retorna o ponteiro físico da GPU para o bloco contíguo da variável vindex
    return d_states + (vindex * num_systems);
}