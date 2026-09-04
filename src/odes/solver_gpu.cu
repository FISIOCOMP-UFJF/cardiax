#include "odes/torord_land_gpu.cuh"
#include <stdio.h>

// Launch bounds limitam o número máximo de registradores por thread.
// Se o limite for excedido, o compilador alerta em tempo de compilação, 
// impedindo que a performance caia silenciosamente por "register spilling" para a memória L1/VRAM.
__global__ 
__launch_bounds__(256, 4) 
void kernel_advance_cells(
    double* d_states, 
    const double* d_istim, 
    const int* d_celltypes,
    const double dt, 
    const int num_systems) 
{
    // Identificação global da thread
    int cell_id = blockIdx.x * blockDim.x + threadIdx.x;

    if (cell_id < num_systems) {
        // Puxa o tipo da célula e o estímulo atual (pode ser 0.0 se não for nó de estímulo)
        int celltype = d_celltypes ? d_celltypes[cell_id] : 0;
        double istim = d_istim ? d_istim[cell_id] : 0.0;

        // O parâmetro pitch é o próprio número de sistemas (num_systems), 
        // mantendo os dados contíguos por variável de estado (SoA).
        compute_and_advance_torord_land(d_states, dt, cell_id, num_systems, istim, celltype);
    }
}

// Wrapper em C++ para ser chamado pelo Host (sua classe CellsGpu)
void launch_euler_kernel(
    double* d_states, 
    const double* d_istim, 
    const int* d_celltypes,
    const double dt, 
    const int num_systems) 
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_systems + threadsPerBlock - 1) / threadsPerBlock;

    kernel_advance_cells<<<blocksPerGrid, threadsPerBlock>>>(
        d_states, d_istim, d_celltypes, dt, num_systems
    );
    
    // Opcional em release, mas essencial para debug inicial
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Erro no Kernel CUDA: %s\n", cudaGetErrorString(err));
    }
}