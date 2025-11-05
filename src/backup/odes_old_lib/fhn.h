#ifndef _IONIC_DEVICE_H_
#define _IONIC_DEVICE_H_

#include <cstdio>
#include <cuda.h>
#include <cutil.h>
#include <cublas.h>
#include "util.h"

// -----------------------------------------------------------------------------
// FitzHugh Nagumo constants
// -----------------------------------------------------------------------------

#define V_star -1.19940803524
#define W_star -0.624260044055

#define EPS 0.2
#define GAMMA 0.8
#define BETA 0.7

// -----------------------------------------------------------------------------
// Modified FitzHugh Nagumo constants (Rogers and McCulloch)
// -----------------------------------------------------------------------------
#define m_V_star 0.0
#define m_W_star 0.0

#define p_a  0.130
#define p_b  0.013
#define p_d  1.000
#define p_c1 0.26
#define p_c2 0.1

// ----------------------------------------------------------------------------
// (Scaled) Modified FitzHugh Nagumo constants (Rogers and McCulloch)
// ----------------------------------------------------------------------------

#define d_G    1.5
#define d_vth  13.0
#define d_vp   100.0
#define d_eta1 4.4
#define d_eta2 0.012
#define d_eta3 1.0

// -----------------------------------------------------------------------------
// Functions
// -----------------------------------------------------------------------------

// declare texture reference for 1D real texture
texture<float> tex_v;
texture<float> tex_w;

// wrappers
void cudaFHN_FE_step      (int it, int n, Real dt, Real *sv_v, Real *sv_w, int useTex);
void cudaFHN_RMcC_FE_step (int it, int n, Real dt, Real *sv_v, Real *sv_w, int useTex);

// FHN time steppers
__global__ void cudaFHNFEstep         (int it, int nodes, Real dt, Real *sv_v, Real *sv_w);
__global__ void cudaFHNFEstepTex1D    (int it, int nodes, Real dt, Real *sv_v, Real *sv_w);

__global__ void cudaFHN_RMcC_FEStep (int it, int nodes, int neq, Real dt, Real *sv);

__global__ void cudaFHN_RMcC_FEStepOld   (int it, int nodes, Real dt, Real *sv_v, Real *sv_w);

__global__ void cudaFHN_ScRMcC_FEStep (int it, int nodes, Real dt, Real *sv_v, Real *sv_w);




#endif
