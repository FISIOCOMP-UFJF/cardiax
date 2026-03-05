#include <stdio.h>
#include "ionic.h"

// ===================================================================
// Standard FitzHugh-Nagumo
// ===================================================================

/*!
 Set the initial conditions for the standard FitzHugh-Nagumo model
 @param[in]  nodes   Number of nodes
 @param[out] sv_v    State variable v
 @param[out] sv_w    State variable w
 */

void FHN_SetIC (int nodes, Real **sv_v, Real **sv_w)
{
  for(int i=0; i<nodes; i++){
    sv_v[0][i] = V_star;
    sv_w[0][i] = W_star;
  }
}

/*!
 Compute an Euler time step for the standard FitzHugh-Nagumo model
 @param[in]  nodes   Number of nodes
 @param[in]  dt      Time step
 @param[out] sv_v    State variable v
 @param[out] sv_w    State variable w
 */

void FHN_FEStep (int it, int nodes, Real dt, Real **sv_v, Real **sv_w)
{
  Real dv_dt;
  Real dw_dt;
  Real v3;
  for (int i=0; i<nodes; i++){
    v3 = pow(sv_v[it-1][i],3);
    dv_dt = ((1.0/EPS) * (sv_v[it-1][i]-((1.0/3.0)*v3)-sv_w[it-1][i]));
    dw_dt = EPS * (sv_v[it-1][i] - GAMMA * sv_w[it-1][i] + BETA);
    sv_v[it][i] = sv_v[it-1][i] + dt * dv_dt;
    sv_w[it][i] = sv_w[it-1][i] + dt * dw_dt;
  }
}

// ===================================================================
// FitzHugh-Nagumo (modified Rogers and McCulloch)
// ===================================================================

/*!
 Set the initial conditions for the modified FitzHugh-Nagumo model
 @param[in]  nodes   Number of nodes
 @param[out] sv_v    State variable v
 @param[out] sv_w    State variable w
 */

void FHN_RMcC_SetIC (int nodes, Real **sv_v, Real **sv_w)
{
  for(int i=0; i<nodes; i++){
    sv_v[0][i] = m_V_star;
    sv_w[0][i] = m_W_star;
  }
}

/*!
 Compute an Euler time step for the modified FitzHugh-Nagumo model
 @param[in]  nodes   Number of nodes
 @param[in]  dt      Time step
 @param[out] sv_v    State variable v
 @param[out] sv_w    State variable w
 */

void FHN_RMcC_FEStep (int it, int nodes, Real dt, Real **sv_v, Real **sv_w)
{
  Real dv_dt;
  Real dw_dt;

  for(int i=0; i<nodes; i++){
      dv_dt = p_c1*sv_v[it-1][i]*(sv_v[it-1][i]-p_a)*(1.0-sv_v[it-1][i])-p_c2*sv_v[it-1][i]*sv_w[it-1][i];
      dw_dt = p_b*(sv_v[it-1][i]-p_d*sv_w[it-1][i]);
      sv_v[it][i] = sv_v[it-1][i] + dt * dv_dt;
      sv_w[it][i] = sv_w[it-1][i] + dt * dw_dt;
  }
}

// ===================================================================
// Scaled FitzHugh-Nagumo (modified Rogers and McCulloch)
// ===================================================================

/*!
 Compute an Euler time step for the modified-scaled FitzHugh-Nagumo model
 @param[in]  nodes   Number of nodes
 @param[in]  dt      Time step
 @param[out] sv_v    State variable v
 @param[out] sv_w    State variable w
 */

void FHN_ScRMcC_FEStep (int it, int nodes, Real dt, Real **sv_v, Real **sv_w)
{
  Real dv_dt;
  Real dw_dt;

  for(int i=0; i<nodes; i++){
      dv_dt = d_G    *  sv_v[it-1][i]*(1.0-(sv_v[it-1][i]/d_vth))*(1.0-(sv_v[it-1][i]/d_vp)) - d_eta1 * sv_v[it-1][i] * sv_w[it-1][i];
      dw_dt = d_eta2 * (sv_v[it-1][i]/d_vp - d_eta3*sv_w[it-1][i]);
      sv_v[it][i] = sv_v[it-1][i] + dt * dv_dt;
      sv_w[it][i] = sv_w[it-1][i] + dt * dw_dt;
  }
}
