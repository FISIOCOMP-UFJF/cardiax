#include "petsc_linear_solver.hpp"
#include "util/command_line_args.h"
#include <chrono>

//#define PETSC_SOLVER
#define AMGX_SOLVER

namespace petsc
{

void callback(const char *msg, int length)
{
        const char *aux = msg;
}

void LinearSolver::init()
{
  ierr = KSPCreate(PETSC_COMM_WORLD, &_ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  set_solver_type();

  // Set runtime options, e.g.,
  // -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  //
  // These options will override those specified above as long as
  // KSPSetFromOptions() is called _after_ any other customization
  // routines.

  ierr = KSPSetFromOptions (_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

#ifdef AMGX_SOLVER
  AMGX_SAFE_CALL(AMGX_register_print_callback(&callback));
#endif
}


void LinearSolver::converged_reason()
{
  KSPConvergedReason kr;
  ierr = KSPGetConvergedReason(_ksp, &kr);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  cout << endl;
  cout << "PETSc KSP Converged Reason = " << kr << endl;
  cout << endl;
}

void LinearSolver::set_solver_type()
{
//  ierr = KSPSetType(_ksp, KSPCG);
//  ierr = KSPSetType(_ksp, KSPGMRES);

  ierr = KSPGetPC(_ksp, &_pc);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

//  ierr = PCSetType(_pc, PCJACOBI);
//  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void LinearSolver::set_solver_type(const char * type)
{
  string ksptype(type);
  set_solver_type(ksptype);
}

void LinearSolver::set_solver_type(std::string & type)
{
  if(type == "cg")
  {
    ierr = KSPSetType(_ksp, KSPCG);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "bicg")
  {
    ierr = KSPSetType(_ksp, KSPBICG);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "gmres")
  {
    ierr = KSPSetType(_ksp, KSPGMRES);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "bcgs")
  {
    ierr = KSPSetType(_ksp, KSPBCGS);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "tfqmr")
  {
    ierr = KSPSetType(_ksp, KSPTFQMR);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "richardson")
  {
    ierr = KSPSetType(_ksp, KSPRICHARDSON);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }

}
  
void LinearSolver::set_preconditioner(const char * type)
{
  string pctype(type);
  set_preconditioner(type);
}

void LinearSolver::set_preconditioner(std::string & type)
{
  ierr = KSPGetPC(_ksp, &_pc);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if(type == "none")
  {
    ierr = PCSetType(_pc, PCNONE);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "jacobi")
  {
    ierr = PCSetType(_pc, PCJACOBI);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "bjacobi")
  {
    ierr = PCSetType(_pc, PCBJACOBI);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "sor")
  {
    ierr = PCSetType(_pc, PCSOR);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "eisenstat")
  {
    ierr = PCSetType(_pc, PCEISENSTAT);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "icc")
  {
    ierr = PCSetType(_pc, PCICC);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  } 
  else if(type == "ilu")
  {
    ierr = PCSetType(_pc, PCILU);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  } 
  else if(type == "asm")
  {
    ierr = PCSetType(_pc, PCASM);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }   
  else if(type == "lu")      // direct solver
  {
    ierr = PCSetType(_pc, PCLU);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if(type == "cholesky") // direct solver
  {
    ierr = PCSetType(_pc, PCCHOLESKY);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
}

void LinearSolver::set_ordering(const char * otype)
{
  string type(otype);

  if (type == "rcm")
  
{    ierr = PCFactorSetMatOrderingType(_pc,"rcm");
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
  else if (type == "amd")
  {
    ierr = PCFactorSetMatOrderingType(_pc,"amd");
    CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }
}


void LinearSolver::MatlabView(petsc::Matrix & A)
{
      /*
      PetscViewer vm;

      PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Amat.m", &vm);

      PetscViewerSetFormat(vm, PETSC_VIEWER_ASCII_MATLAB);

      MatView(A.mat(),vm);
      PetscViewerDestroy(&vm);
      */
}

std::pair<PetscInt, PetscReal> LinearSolver::solveFieldSplit (petsc::Matrix & A,
                                                    petsc::Vector & x,
                                                    petsc::Vector & b,
                                                    const double tol)
{
  PetscInt its;
  PetscReal rnorm;
  KSPConvergedReason reason;
  int N = A.size();

  std::cout<<"Solve using Field Split\n";

  // Create Index Sets
  PetscInt T_indices[N/2];
  PetscInt P_indices[N/2];

  for (int i = 0; i < N/2; ++i)
  {
    P_indices[i] = i;
    T_indices[i] = i + N/2;
  }

  IS P_IS, T_IS;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/2, P_indices, PETSC_COPY_VALUES, &P_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/2, T_indices, PETSC_COPY_VALUES, &T_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);

// Solve the system
  //ierr = KSPCreate(PETSC_COMM_WORLD, &_ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetType(_ksp, KSPGMRES);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetOperators(_ksp, A.mat(), A.mat());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetFromOptions(_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPGetPC(_ksp, &_pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);

// Define the fieldsplit for the global matrix
  ierr = PCFieldSplitSetIS(_pc, "P", P_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PCFieldSplitSetIS(_pc, "T", T_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetUp(_ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSolve(_ksp,b.vec(),x.vec());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetIterationNumber(_ksp,&its);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetResidualNorm(_ksp, &rnorm);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetConvergedReason(_ksp, &reason);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if (reason < 0)
  {
    cout << "PETSC Converged Reason = " << reason << endl;
    print_error("petsc_linear_solver.cpp", "solve", "KSP not converged");
  }

  std::pair<PetscInt,PetscReal> ir(its,rnorm);
  return ir;
}

std::pair<PetscInt, PetscReal> LinearSolver::solve_3_FieldSplit (petsc::Matrix & A,
                                                              petsc::Vector & x,
                                                              petsc::Vector & b,
                                                              const double tol)
{
  PetscInt its;
  PetscReal rnorm;
  KSPConvergedReason reason;
  int N = A.size();

  //std::cout<<"Solve using 3 Field Split\n";

  // Create Index Sets
  PetscInt X_indices[N/3];
  PetscInt Y_indices[N/3];
  PetscInt Z_indices[N/3];

  for (int i = 0; i < N/3; ++i)
  {
    X_indices[i] = i;
    Y_indices[i] = i + N/3;
    Z_indices[i] = i + 2*N/3;
  }

  IS X_IS, Y_IS, Z_IS;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/3, X_indices, PETSC_COPY_VALUES, &X_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/3, Y_indices, PETSC_COPY_VALUES, &Y_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/3, Z_indices, PETSC_COPY_VALUES, &Z_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

// Solve the system
  ierr = KSPCreate(PETSC_COMM_WORLD, &_ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetType(_ksp, KSPGMRES);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetOperators(_ksp, A.mat(), A.mat());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetFromOptions(_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPGetPC(_ksp, &_pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);

// Define the fieldsplit for the global matrix
  ierr = PCFieldSplitSetIS(_pc, "X", X_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PCFieldSplitSetIS(_pc, "Y", Y_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PCFieldSplitSetIS(_pc, "Z", Z_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPSetUp(_ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);

  //ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  //CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSolve(_ksp,b.vec(),x.vec());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetIterationNumber(_ksp,&its);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetResidualNorm(_ksp, &rnorm);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetConvergedReason(_ksp, &reason);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = ISDestroy(&X_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISDestroy(&Y_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISDestroy(&Z_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPDestroy(&_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if (reason < 0)
  {
    cout << "PETSC Converged Reason = " << reason << endl;
    print_error("petsc_linear_solver.cpp", "solve", "KSP not converged");
  }

  std::pair<PetscInt,PetscReal> ir(its,rnorm);
  return ir;
}



const char* LinearSolver::getConfig(const std::string &configName)
{
    const char *cfg;
    string cfgString;

    if (configName == "AMG_CLASSICAL_AGGRESSIVE_CHEB_L1_TRUNC"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"scope\": \"main\", \"solver\": \"PCG\", \"store_res_history\": 1, \"print_solve_stats\": 0, \"obtain_timings\": 1, \"preconditioner\": {\"print_grid_stats\": 0, \"scope\": \"amg_solver\", \"interpolator\": \"D2\", \"solver\": \"AMG\", \"max_levels\": 24, \"selector\": \"PMIS\", \"cycle\": \"V\", \"presweeps\": 0, \"postsweeps\": 3, \"coarsest_sweeps\": 2, \"min_coarse_rows\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, \"error_scaling\": 3, \"print_grid_stats\": 0, \"aggressive_levels\": 1, \"interp_max_elements\": 4, \"smoother\": {\"relaxation_factor\": 0.91, \"scope\": \"jacobi\", \"solver\": \"CHEBYSHEV\", \"max_iters\": 1, \"preconditioner\": {\"solver\": \"JACOBI_L1\", \"max_iters\": 1}, \"chebyshev_polynomial_order\": 2, \"chebyshev_lambda_estimate_mode\": 2}}, \"max_iters\": 100, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AGGREGATION_MULTI_PAIRWISE"){
    
        cfgString = "{\"config_version\": 2, \"solver\": {\"use_scalar_norm\": 1, \"min_coarse_rows\": 32, \"matrix_coloring_scheme\": \"PARALLEL_GREEDY\", \"smoother\": \"BLOCK_JACOBI\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"relaxation_factor\": 0.5, \"weight_formula\": 1, \"tolerance\": 0.001, \"norm\": \"L2\", \"postsweeps\": 3, \"merge_singletons\": 2, \"presweeps\": 0, \"selector\": \"MULTI_PAIRWISE\", \"scope\": \"main\", \"convergence\": \"RELATIVE_MAX_CORE\", \"cycle\": \"F\", \"print_grid_stats\": 1, \"algorithm\": \"AGGREGATION\", \"solver\": \"AMG\", \"aggregation_passes\": 1, \"max_uncolored_percentage\": 0.05, \"coarse_solver\": \"DENSE_LU_SOLVER\", \"max_iters\": 1, \"monitor_residual\": 1, \"max_levels\": 20}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AGGREGATION_THRUST_BJ"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"algorithm\": \"AGGREGATION\", \"coarseAgenerator\": \"THRUST\", \"solver\": \"AMG\", \"smoother\": \"BLOCK_JACOBI\", \"print_solve_stats\": 1, \"presweeps\": 1, \"selector\": \"SIZE_2\", \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"scope\": \"main\", \"postsweeps\": 1, \"tolerance\": 0.1, \"cycle\": \"V\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AGGREGATION_THRUST_DILU"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"obtain_timings\": 1, \"print_grid_stats\": 1, \"max_uncolored_percentage\": 0.02, \"algorithm\": \"AGGREGATION\", \"coarseAgenerator\": \"THRUST\", \"solver\": \"AMG\", \"smoother\": \"MULTICOLOR_DILU\", \"print_solve_stats\": 1, \"presweeps\": 1, \"selector\": \"SIZE_2\", \"coarsest_sweeps\": 2, \"max_iters\": 2000, \"monitor_residual\": 1, \"postsweeps\": 1, \"scope\": \"main\", \"max_levels\": 1000, \"matrix_coloring_scheme\": \"MIN_MAX\", \"tolerance\": 0.1, \"norm\": \"L1\", \"cycle\": \"V\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AGGREGATION_THRUST_GS"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"coarsest_sweeps\": 2, \"matrix_coloring_scheme\": \"MIN_MAX\", \"smoother\": \"MULTICOLOR_GS\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"relaxation_factor\": 0.75, \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L1\", \"postsweeps\": 2, \"max_uncolored_percentage\": 0.15, \"presweeps\": 2, \"selector\": \"SIZE_2\", \"convergence\": \"RELATIVE_INI\", \"cycle\": \"V\", \"print_grid_stats\": 1, \"algorithm\": \"AGGREGATION\", \"coarseAgenerator\": \"THRUST\", \"solver\": \"AMG\", \"max_iters\": 50, \"monitor_residual\": 1, \"max_levels\": 1000}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "agg_cheb4"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"print_grid_stats\": 1, \"algorithm\": \"AGGREGATION\", \"obtain_timings\": 1, \"error_scaling\": 3, \"solver\": \"AMG\", \"smoother\": {\"solver\": \"CHEBYSHEV\", \"preconditioner\": {\"solver\": \"JACOBI_L1\", \"max_iters\": 1}, \"max_iters\": 1, \"chebyshev_polynomial_order\": 4, \"chebyshev_lambda_estimate_mode\": 2}, \"presweeps\": 0, \"postsweeps\": 1, \"print_solve_stats\": 1, \"selector\": \"SIZE_8\", \"coarsest_sweeps\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"min_coarse_rows\": 2, \"scope\": \"main\", \"max_levels\": 1000, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"cycle\": \"V\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_AGGRREGATION_CG"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"algorithm\": \"AGGREGATION\", \"selector\": \"SIZE_4\", \"print_solve_stats\": 1, \"smoother\": \"JACOBI_L1\", \"presweeps\": 0, \"postsweeps\": 3, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"CG\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_CLASSICAL_AGGRESSIVE_L1"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"solver\": \"AMG\", \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 1, \"selector\": \"PMIS\", \"coarsest_sweeps\": 1, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 1}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 10, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_CLASSICAL_AGGRESSIVE_L1_TRUNC"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"solver\": \"AMG\", \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"interp_max_elements\": 4, \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 2, \"selector\": \"PMIS\", \"coarsest_sweeps\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 2}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 10, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_CLASSICAL_CG"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"interpolator\": \"D2\", \"presweeps\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"CG\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_CLASSICAL_L1_AGGRESSIVE_HMIS"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"solver\": \"AMG\", \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 1, \"selector\": \"HMIS\", \"coarsest_sweeps\": 1, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 1}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 10, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_CLASSICAL_L1_TRUNC"){
        
        cfgString = "{\"config_version\": 2, \"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": {\"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"print_grid_stats\": 1, \"solver\": \"AMG\", \"interp_max_elements\": 4, \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 1, \"coarsest_sweeps\": 1, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"scope\": \"amg_solver\", \"max_levels\": 25, \"postsweeps\": 1, \"cycle\": \"V\"}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 10, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "AMG_CLASSICAL_PMIS"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"obtain_timings\": 1, \"solver\": \"GMRES\", \"print_solve_stats\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"solver\": \"AMG\", \"cycle\": \"V\", \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 2, \"postsweeps\": 2, \"selector\": \"PMIS\", \"coarsest_sweeps\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"interp_max_elements\": 4}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 10, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "CG_DILU"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"precond\", \"solver\": \"MULTICOLOR_DILU\"}, \"solver\": \"CG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 20, \"monitor_residual\": 1, \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "CHEB_SOLVER_NOPREC"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"solver\": \"NOSOLVER\", \"scope\": \"precond\"}, \"use_scalar_norm\": 1, \"solver\": \"CHEBYSHEV\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"max_iters\": 2000, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "CLASSICAL_CG_CYCLE"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"interpolator\": \"D2\", \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"CG\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "CLASSICAL_F_CYCLE"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"obtain_timings\": 1, \"interpolator\": \"D2\", \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"F\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "CLASSICAL_V_CYCLE"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"interpolator\": \"D2\", \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"V\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "CLASSICAL_W_CYCLE"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"interpolator\": \"D2\", \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"W\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "F"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"interpolator\": \"D2\", \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"F\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "FGMRES"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"print_solve_stats\": 0, \"interpolator\": \"D2\", \"presweeps\": 1, \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 20, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 0.0001, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "FGMRES_AGGREGATION_JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"error_scaling\": 0, \"print_grid_stats\": 1, \"algorithm\": \"AGGREGATION\", \"solver\": \"AMG\", \"smoother\": \"BLOCK_JACOBI\", \"presweeps\": 0, \"selector\": \"SIZE_2\", \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"min_coarse_rows\": 32, \"relaxation_factor\": 0.75, \"scope\": \"amg\", \"max_levels\": 100, \"postsweeps\": 3, \"cycle\": \"V\"}, \"use_scalar_norm\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 32, \"convergence\": \"RELATIVE_INI\", \"scope\": \"main\", \"tolerance\": 0.0001, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "FGMRES_CLASSICAL_AGGRESSIVE_HMIS"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"solver\": \"AMG\", \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"interp_max_elements\": 4, \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 2, \"selector\": \"HMIS\", \"coarsest_sweeps\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 2}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 100, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "FGMRES_CLASSICAL_AGGRESSIVE_PMIS"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"solver\": \"AMG\", \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"interp_max_elements\": 4, \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 2, \"selector\": \"PMIS\", \"coarsest_sweeps\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 2}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 100, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "FGMRES_NOPREC"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"amg\", \"solver\": \"NOSOLVER\"}, \"use_scalar_norm\": 1, \"solver\": \"FGMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 0.1, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "GMRES"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"amg\", \"solver\": \"NOSOLVER\"}, \"use_scalar_norm\": 1, \"solver\": \"GMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "GMRES_AMG_D2"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"exception_handling\": 1, \"solver\": {\"print_grid_stats\": 1, \"store_res_history\": 1, \"solver\": \"GMRES\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"preconditioner\": {\"interpolator\": \"D2\", \"print_grid_stats\": 1, \"solver\": \"AMG\", \"smoother\": \"JACOBI_L1\", \"presweeps\": 2, \"selector\": \"PMIS\", \"coarsest_sweeps\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"interp_max_elements\": 4, \"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 2}, \"max_iters\": 2000, \"monitor_residual\": 1, \"gmres_n_restart\": 10, \"convergence\": \"RELATIVE_INI_CORE\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "IDRMSYNC_DILU"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"max_iters\": 4000, \"preconditioner\": {\"scope\": \"precond\", \"max_iters\": 1, \"solver\": \"MULTICOLOR_DILU\"}, \"use_scalar_norm\": 1, \"solver\": \"IDRMSYNC\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"subspace_dim_s\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 1e-08, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "IDR_DILU"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"max_iters\": 4000, \"print_grid_stats\": 1, \"preconditioner\": {\"scope\": \"precond\", \"max_iters\": 1, \"solver\": \"MULTICOLOR_DILU\"}, \"use_scalar_norm\": 1, \"solver\": \"IDR\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"subspace_dim_s\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 1e-08, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"amg\", \"solver\": \"NOSOLVER\"}, \"use_scalar_norm\": 1, \"solver\": \"BLOCK_JACOBI\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PBICGSTAB"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"interpolator\": \"D2\", \"presweeps\": 1, \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"PBICGSTAB\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PBICGSTAB_CLASSICAL_JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 1, \"interpolator\": \"D2\", \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"PBICGSTAB\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PBICGSTAB_NOPREC"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"amg_solver\", \"solver\": \"NOSOLVER\"}, \"use_scalar_norm\": 1, \"solver\": \"PBICGSTAB\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 0.1, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCGF_CLASSICAL_V_JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 1, \"interpolator\": \"D2\", \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"PCGF\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCGF_V"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 1, \"max_iters\": 1, \"interpolator\": \"D2\", \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"PCGF\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_AGGREGATION_JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"algorithm\": \"AGGREGATION\", \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"relaxation_factor\": 0.8, \"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 0, \"interpolator\": \"D2\", \"selector\": \"SIZE_2\", \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"postsweeps\": 3, \"cycle\": \"V\"}, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI\", \"scope\": \"main\", \"tolerance\": 0.0001, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_CLASSICAL_F_JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 1, \"max_iters\": 1, \"interpolator\": \"D2\", \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"F\", \"postsweeps\": 1}, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_CLASSICAL_V_JACOBI"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"aggressive_levels\": 2, \"presweeps\": 1, \"interpolator\": \"D2\", \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_DILU"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"precond\", \"solver\": \"MULTICOLOR_DILU\"}, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 20, \"monitor_residual\": 1, \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_F"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 1, \"interpolator\": \"D2\", \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"F\", \"postsweeps\": 1}, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_NOPREC"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"amg\", \"solver\": \"NOSOLVER\"}, \"use_scalar_norm\": 1, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"monitor_residual\": 1, \"convergence\": \"RELATIVE_INI_CORE\", \"scope\": \"main\", \"tolerance\": 0.1, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "PCG_V"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"print_grid_stats\": 1, \"print_vis_data\": 0, \"solver\": \"AMG\", \"smoother\": {\"scope\": \"jacobi\", \"solver\": \"BLOCK_JACOBI\", \"monitor_residual\": 0, \"print_solve_stats\": 0}, \"print_solve_stats\": 0, \"presweeps\": 1, \"interpolator\": \"D2\", \"max_iters\": 1, \"monitor_residual\": 0, \"store_res_history\": 0, \"scope\": \"amg\", \"max_levels\": 100, \"cycle\": \"V\", \"postsweeps\": 1}, \"solver\": \"PCG\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"tolerance\": 1e-06, \"norm\": \"L2\"}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "V-cheby-aggres-L1-trunc-userLambda"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"interpolator\": \"D2\", \"interp_max_elements\": 4, \"aggressive_levels\": 1, \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"coarse_solver\": \"NOSOLVER\", \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"error_scaling\": 3, \"cycle\": \"V\", \"smoother\": {\"solver\": \"CHEBYSHEV\", \"preconditioner\": {\"solver\": \"JACOBI_L1\", \"max_iters\": 1}, \"max_iters\": 1, \"chebyshev_polynomial_order\": 2, \"chebyshev_lambda_estimate_mode\": 3, \"cheby_max_lambda\": 0.92, \"cheby_min_lambda\": 0.07}, \"tolerance\": 1e-06, \"norm\": \"L2\", \"presweeps\": 0, \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "V-cheby-aggres-L1-trunc"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"interpolator\": \"D2\", \"interp_max_elements\": 4, \"aggressive_levels\": 1, \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"coarse_solver\": \"NOSOLVER\", \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"error_scaling\": 3, \"cycle\": \"V\", \"smoother\": {\"solver\": \"CHEBYSHEV\", \"preconditioner\": {\"solver\": \"JACOBI_L1\", \"max_iters\": 1}, \"max_iters\": 1, \"chebyshev_polynomial_order\": 2, \"chebyshev_lambda_estimate_mode\": 2}, \"tolerance\": 1e-06, \"norm\": \"L2\", \"presweeps\": 0, \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "V-cheby-smoother"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"scope\": \"main\", \"print_grid_stats\": 1, \"solver\": \"AMG\", \"scaling\": \"DIAGONAL_SYMMETRIC\", \"interpolator\": \"D2\", \"aggressive_levels\": 0, \"interp_max_elements\": 4, \"coarse_solver\": \"NOSOLVER\", \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"max_levels\": 100, \"cycle\": \"V\", \"smoother\": {\"solver\": \"CHEBYSHEV\", \"preconditioner\": {\"solver\": \"NOSOLVER\", \"max_iters\": 1}, \"max_iters\": 1, \"chebyshev_polynomial_order\": 4, \"chebyshev_lambda_estimate_mode\": 2}, \"tolerance\": 1e-06, \"norm\": \"L2\", \"presweeps\": 0, \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "V-cheby_poly-smoother"){
        
        cfgString = "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"scope\": \"main\", \"print_grid_stats\": 1, \"solver\": \"AMG\", \"scaling\": \"NONE\", \"interpolator\": \"D2\", \"aggressive_levels\": 0, \"interp_max_elements\": 4, \"max_row_sum\": 0.9, \"print_solve_stats\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"max_levels\": 100, \"cycle\": \"V\", \"smoother\": {\"scope\": \"cheb_smoother\", \"solver\": \"CHEBYSHEV_POLY\", \"chebyshev_polynomial_order\": 2, \"max_iters\": 1, \"preconditioner\": {\"solver\": \"JACOBI_L1\"}, \"tolerance\": 0.0001}, \"tolerance\": 1e-06, \"norm\": \"L2\", \"presweeps\": 0, \"postsweeps\": 3}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "V"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"interpolator\": \"D2\", \"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"V\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else if (configName == "W"){
        
        cfgString = "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 1, \"solver\": \"AMG\", \"print_solve_stats\": 1, \"presweeps\": 1, \"obtain_timings\": 1, \"max_iters\": 2000, \"monitor_residual\": 1, \"interpolator\": \"D2\", \"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"max_levels\": 100, \"cycle\": \"W\", \"tolerance\": 1e-06, \"norm\": \"L2\", \"postsweeps\": 1}}";
        cfg = cfgString.c_str();
    }
    else
    {
        std::cerr << "Error: Unknown configuration name." << std::endl;
        return nullptr;
    }

    return cfg;
}


#ifdef AMGX_SOLVER
std::pair<PetscInt, PetscReal> LinearSolver::solve(petsc::Matrix &A,
                                                  petsc::Vector &x,
                                                  petsc::Vector &b,
                                                  const double tol)
{
    // Specify the desired AMGX configuration name
    const std::string amgxConfigName = CommandLineArgs::read("-amgx", "./configs/AMG_CLASSICAL_PMIS.json");

    int its;
    double rnorm;

    // Initialize AMGX
    AMGX_SAFE_CALL(AMGX_initialize());

    // Create and configure AMGX resources and solver
    AMGX_config_handle config;
    AMGX_resources_handle rsrc;
    AMGX_solver_handle amgx_solver;

    //const char *amgxConfig = getConfig(amgxConfigName);

    std::string config_path_str = amgxConfigName;
    const char* config_path = config_path_str.c_str();

    //AMGX_SAFE_CALL(AMGX_config_create(&config, amgxConfig));
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&config, config_path))
    AMGX_SAFE_CALL(AMGX_resources_create_simple(&rsrc, config));

    // Create and configure AMGX matrix
    AMGX_matrix_handle amgx_A;
    AMGX_SAFE_CALL(AMGX_matrix_create(&amgx_A, rsrc, AMGX_mode_dDDI));

    AMGX_vector_handle amgx_b, amgx_x;
    AMGX_SAFE_CALL(AMGX_vector_create(&amgx_b, rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_vector_create(&amgx_x, rsrc, AMGX_mode_dDDI));

    AMGX_SAFE_CALL(AMGX_solver_create(&amgx_solver, rsrc, AMGX_mode_dDDI, config));

    int n, nnz, *row_ptrs, *col_indices;
    double *values;
    nnz = A.get_nnz();
    n = A.size();

    row_ptrs = new int[n + 1];
    col_indices = new int[nnz];
    values = new double[nnz];

    A.get_CSR(&n, row_ptrs, col_indices, values);

    for(int i = 0; i < n + 1; i++)
    {
      row_ptrs[i] -= 1;
    }

    for(int i = 0; i < nnz; i++)
    {
      col_indices[i] -= 1;
    }


    AMGX_SAFE_CALL(AMGX_matrix_upload_all(amgx_A, n, nnz, 1, 1, row_ptrs, col_indices, values, NULL));

    // Create AMGX vectors for b and x
    double *b_values, *x_values;

    b_values = b.get_array();
    x_values = x.get_array();
    AMGX_SAFE_CALL(AMGX_vector_upload(amgx_b, b.size(), 1, b_values));
    AMGX_SAFE_CALL(AMGX_vector_upload(amgx_x, x.size(), 1, x_values));

    // Set up the AMGX solver
    AMGX_SAFE_CALL(AMGX_solver_setup(amgx_solver, amgx_A));

    // Solve the system using AMGX
    AMGX_SAFE_CALL(AMGX_solver_solve(amgx_solver, amgx_b, amgx_x));

    // Copy the solution back to Petsc
    double *data = new double[x.size()];
    AMGX_SAFE_CALL(AMGX_vector_download(amgx_x, data));
    x.set_data(data);
    delete[] data;
    // Get iteration count and residual norm
    AMGX_SAFE_CALL(AMGX_solver_get_iterations_number(amgx_solver, &its)); // alt: AMGX_solver_get_iterations_number(AMGX_solver_handle obj, int *n)
    AMGX_SAFE_CALL(AMGX_solver_get_iteration_residual(amgx_solver, its, 0, &rnorm)); // AMGX_solver_get_iteration_residual(AMGX_solver_handle obj, int iter,int idx, double *res);

    // Cleanup and finalize AMGX
    AMGX_SAFE_CALL(AMGX_solver_destroy(amgx_solver));
    AMGX_SAFE_CALL(AMGX_vector_destroy(amgx_x));
    AMGX_SAFE_CALL(AMGX_vector_destroy(amgx_b));
    AMGX_SAFE_CALL(AMGX_matrix_destroy(amgx_A));
    AMGX_SAFE_CALL(AMGX_resources_destroy(rsrc));
    AMGX_SAFE_CALL(AMGX_config_destroy(config));
    AMGX_SAFE_CALL(AMGX_finalize());


    delete[] row_ptrs;
    delete[] col_indices;
    delete[] values;

    // Tamanho do array de caracteres para o nome do arquivo
    //const int filenameSize = 100; // Ajuste conforme necessário

    // Array de caracteres para o nome do arquivo
    //char filename[filenameSize];

    // Construindo o nome do arquivo
    //snprintf(filename, filenameSize, "solution_%s.txt", amgxConfigName.c_str());

    // Chamando a função save_array com o array de caracteres
    //x.save_array(filename);

    // Return iteration count and residual norm
    return std::make_pair((PetscInt)its, (PetscReal)rnorm);
}
#endif

#ifdef PETSC_SOLVER


std::pair<PetscInt, PetscReal> LinearSolver::solve (petsc::Matrix & A,
                                                    petsc::Vector & x,
                                                    petsc::Vector & b,
                                                    const double tol)
{
  PetscInt its;
  PetscReal rnorm;
  KSPConvergedReason reason;

  // old code
  //ierr = KSPSetOperators(_ksp, A.mat(), A.mat(), SAME_NONZERO_PATTERN);

  ierr = KSPSetOperators(_ksp, A.mat(), A.mat());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // use NON-zero initial guess
//  ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
//  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSetFromOptions(_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSetUp(_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSolve(_ksp,b.vec(),x.vec());
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetIterationNumber(_ksp,&its);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetResidualNorm(_ksp, &rnorm);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPGetConvergedReason(_ksp, &reason);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  if (reason < 0)
  {
    cout << "PETSC Converged Reason = " << reason << endl;
    print_error("petsc_linear_solver.cpp", "solve", "KSP not converged");    
  }
  
  std::pair<PetscInt,PetscReal> ir(its,rnorm);

  //save matrix in CSR format
  //A.store_CSR();
  //b.save_array("rhs.txt");
  //x.save_array("solution_petsc.txt");

  return ir;
}

#endif

   
void LinearSolver::use_umfpack()
{
  // command line arguments
  //     -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package umfpack

  // -ksp_type preonly
  ierr = KSPSetType(_ksp, KSPPREONLY);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // -pc_type lu
  ierr = PCSetType(_pc, PCLU);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
 
  // -pc_factor_mat_solver_package umfpack
  ierr = PCFactorSetMatSolverPackage(_pc, "umfpack");
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void LinearSolver::use_mumps()
{
  // command line arguments
  //     -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps

  // -ksp_type preonly
  ierr = KSPSetType(_ksp, KSPPREONLY);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // -pc_type lu
  ierr = PCSetType(_pc, PCLU);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = PCFactorSetMatSolverPackage(_pc, "mumps");
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}


void LinearSolver::use_fieldsplit(int N)
{
  /*
   This option should be used with the following command line arguments
    -pc_type fieldsplit
    -fieldsplit_X_pc_type hypre -fieldsplit_X_pc_hypre_type boomeramg \
                                -fieldsplit_X_ksp_type gmres\
                                -fieldsplit_X_ksp_rtol 0.001

    -fieldsplit_Y_pc_type hypre -fieldsplit_Y_pc_hypre_type boomeramg \
                                -fieldsplit_Y_ksp_type gmres \
                                -fieldsplit_Y_ksp_rtol 0.001

    -fieldsplit_Z_pc_type hypre -fieldsplit_Z_pc_hypre_type boomeramg \
                                -fieldsplit_Z_ksp_type gmres \
                                -fieldsplit_Z_ksp_rtol 0.001

    or at least some variation of it.
   */
  std::cout << "Solve using 3 Field Split" << std::endl;

  // create arrays with dofs numbers
  PetscInt xidxs[N/3], yidxs[N/3], zidxs[N/3];
  for (int i = 0; i < N/3; ++i)
  {
    xidxs[i] = i;         // x-displacement dofs
    yidxs[i] = i + N/3;   // y-displacement dofs
    zidxs[i] = i + 2*N/3; // z-displacement dofs
  }

  // create PETSc IndexSets (IS)
  IS X_IS, Y_IS, Z_IS;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/3, xidxs, PETSC_COPY_VALUES, &X_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/3, yidxs, PETSC_COPY_VALUES, &Y_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/3, zidxs, PETSC_COPY_VALUES, &Z_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Solve the system
  //ierr = KSPCreate(PETSC_COMM_WORLD, &_ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // setup solver
  ierr = KSPSetType(_ksp, KSPGMRES);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSetFromOptions(_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // setup preconditioner
  ierr = KSPGetPC(_ksp, &_pc);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // define the fieldsplit for the global matrix
  ierr = PCFieldSplitSetIS(_pc, "X", X_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PCFieldSplitSetIS(_pc, "Y", Y_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PCFieldSplitSetIS(_pc, "Z", Z_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void LinearSolver::use_fieldsplit_2D(int N)
{
  std::cout << "Solve using 2 Field Split\n";

  // Create Index Sets
  PetscInt X_indices[N/2];
  PetscInt Y_indices[N/2];

  for (int i = 0; i < N/2; ++i)
  {
    X_indices[i] = i;
    Y_indices[i] = i + N/2;
  }

  IS X_IS, Y_IS;  
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/2, X_indices, PETSC_COPY_VALUES, &X_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);  
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, N/2, Y_indices, PETSC_COPY_VALUES, &Y_IS);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Solve the system
  //ierr = KSPCreate(PETSC_COMM_WORLD, &_ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  //ierr = KSPSetType(_ksp, KSPGMRES);CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = KSPSetFromOptions(_ksp);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = KSPGetPC(_ksp, &_pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Define the fieldsplit for the global matrix
  ierr = PCFieldSplitSetIS(_pc, "X", X_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PCFieldSplitSetIS(_pc, "Y", Y_IS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

void LinearSolver::use_hypre()
{
  // cmdline: 

  // -ksp_type preonly
  //ierr = KSPSetType(_ksp, KSPGMRES);
  //CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // -pc_type lu
  //ierr = PCSetType(_pc, PCHYPRE);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void LinearSolver::view()
{
  ierr = KSPView(_ksp,PETSC_VIEWER_STDOUT_WORLD);
  CHKERRABORT(PETSC_COMM_WORLD,ierr);
}

} // namespace PETSc
