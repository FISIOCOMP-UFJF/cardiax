#ifndef IONIC_H
#define IONIC_H

#include <cmath>
#include "../util/datatype.h"

typedef enum {LR1, FHN, TT2} ionic_model;

// -------------------------------------------------------------------
// FitzHugh Nagumo constants
// -------------------------------------------------------------------

#define V_star -1.19940803524
#define W_star -0.624260044055

#define EPS 0.2
#define GAMMA 0.8
#define BETA 0.7

// -------------------------------------------------------------------
// Modified FitzHugh Nagumo constants (Rogers and McCulloch)
// -------------------------------------------------------------------
#define m_V_star 0.0
#define m_W_star 0.0

#define p_a  0.130
#define p_b  0.013
#define p_d  1.000
#define p_c1 0.26
#define p_c2 0.1

// -------------------------------------------------------------------
// (Scaled) Modified FitzHugh Nagumo constants (Rogers and McCulloch)
// -------------------------------------------------------------------
#define m_v_star 0.0
#define m_w_star 0.0

#define d_G    1.5
#define d_vth  13.0
#define d_vp   100.0
#define d_eta1 4.4
#define d_eta2 0.012
#define d_eta3 1.0

// -------------------------------------------------------------------
// Functions prototypes
// -------------------------------------------------------------------

void FHN_SetIC  (int nodes, Real **sv_v, Real **sv_w);
void FHN_FEStep (int it, int nodes, Real dt, Real **sv_v, Real **sv_w);

void FHN_RMcC_SetIC  (int nodes, Real **sv_v, Real **sv_w);
void FHN_RMcC_FEStep (int it, int nodes, Real dt, Real **sv_v, Real **sv_w);

void FHN_ScRMcC_FEStep (int it, int nodes, Real dt, Real **sv_v, Real **sv_w);

// New interface
void FHN_RMcC_SetICNew (int nodes, int neq, Real *sv);
void FHN_RMcC_FEStepNew (int it, int nodes, int neq, Real dt, Real *sv);

#endif
