================================================================================
                                   EXAMPLES
================================================================================

CardiaX provides multiple executables, each serving a specific function in 
cardiac simulation. This document serves as a quick guide on how to use them.

TIP: We strongly recommend opening and inspecting the structure of the .xml 
files located in this directory. They contain crucial simulation parameters 
that you will likely need to modify for your specific use cases.

--------------------------------------------------------------------------------
1. MONODOMAIN
--------------------------------------------------------------------------------
The 'monodomain' executable is an electrophysiology simulator that solves the 
reaction-diffusion equation (monodomain model).

USAGE:
  ../build/app/monodomain -f ./monodomain_torord.xml -d 0.1 -t 100 -c NP -m ExplicitEuler -ep monodomain

FLAGS:
  -f   : Specifies the input configuration file (e.g., mesh and parameters).
  -d   : Sets the time integration step (dt).
  -t   : Sets the final simulation time.
  -c   : Specifies the cell model (e.g., NP for Nash-Panfilov, TT for Ten Tusscher).
  -m   : Specifies the ODE solver method (e.g., ExplicitEuler, RungeKutta).
  -ep  : Selects the physics model (use 'monodomain' or 'bidomain').

NOTE: To view all available cell models, check the source code in the 'src/odes' directory.

--------------------------------------------------------------------------------
2. NONLINEARELAS
--------------------------------------------------------------------------------
The 'nonlinearelas' executable simulates the passive mechanical behavior of 
cardiac tissue (elasticity only, without electrical coupling).

USAGE:
  ../build/app/nonlinearelas -m mesh.xml -s ul -amgx ../configs/CG_DILU.json

FLAGS:
  -m   : Specifies the mesh/configuration file.
  -s   : Specifies the formulation strategy (e.g., 'ul' for Updated Lagrangian).
  -amgx: Specifies the path to the AMGX configuration JSON file (linear solver settings).

--------------------------------------------------------------------------------
3. ELECTROMECH
--------------------------------------------------------------------------------


USAGE:
  ../build/app/electromech -f ./examples/pvloop.xml -s ul -amgx ../configs/CG_DILU.json

FLAGS:
  -f   : Specifies the input configuration file.
  -s   : Specifies the mechanical formulation (e.g., 'ul').
  -amgx: Specifies the path to the AMGX configuration JSON file (linear solver settings).

[USING EIKONAL ACTIVATION]
To apply different local activation times for each cell, ensure your mesh contains the 'eikonal' 
field inside 'element_data', as demonstrated in this example:

USAGE:
  ../build/app/electromech -f ./examples/pvloop_eikonal.xml -s ul -amgx ../configs/CG_DILU.json 

--------------------------------------------------------------------------------
4. COUPLED ELECTROMECH
--------------------------------------------------------------------------------
* WORK IN PROGRESS *


--------------------------------------------------------------------------------
5. PETSC TIPS
--------------------------------------------------------------------------------
Below are command line options for 3D elasticity preconditioners. You can pass
these directly to the CardiaX executables.

[BOOMERAMG]
  -ksp_type gmres -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_max_iter 1 -ksp_view -ksp_monitor

[GAMG]
Option 1:
  -ksp_type cg -pc_type gamg -pc_gamg_type agg -log_summary -ksp_monitor -ksp_view -options_left -mg_levels_ksp_max_it 1

OR

Option 2:
  -ksp_type cg -pc_type gamg -pc_gamg_type agg -log_summary -ksp_monitor -ksp_view -options_left -mg_levels_ksp_type richardson -mg_levels_pc_type sor