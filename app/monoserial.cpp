// Finite Element Solver for Monodomain Model
// Bernardo M. Rocha

#include <string>

#include "util/util.hpp"
#include "util/igb.hpp"
#include "util/datatype.hpp"
#include "util/timer.hpp"
#include "util/timing.hpp"
#include "util/time_parameters.hpp"
#include "odes/cells.h"
#include "odes/cellmodel.h"
//#include "odes/cellmodel_factory.h"
#include "odes/ode_solver.h"
#include "odes/luo_rudy.h"
#include "odes/ten_tusscher.h"
#include "fem/fem.h"
#include "fem/vtkio.hpp"
#include "cardiac/stimulus.hpp"
#include "linalg/linalg.h"

using namespace std;

// some definitions
#define MAXSTR 256
#define SURF2VOL   0.14       // surface to volume constant
#define UM2_TO_CM2 0.00000001 // unit conversion
#define CSR                   // sparse matrix format

// prototypes
void
monodomain (char *filename, Real dt, Real ftime, Real tp,
	    char *outfile);


static const int tsize=9;

#ifdef USE_PETSC
static char help[] = "Monodomain Solver.\n\n";
#endif

/*
 static char title[tsize][80] = {
  "                                                                       ",
  "  .oooooo.                            .o8   o8o                        ",
  " d8P'  `Y8b                          \"888   `\"'                      ",
  "888           .oooo.   oooo d8b  .oooo888  oooo   .oooo.   oooo    ooo ",
  "888          `P  )88b  `888\"\"8P d88' `888  `888  `P  )88b   `88b..8P'",
  "888           .oP\"888   888     888   888   888   .oP\"888     Y888'  ",
  "`88b    ooo  d8(  888   888     888   888   888  d8(  888   .o8\"'88b  ",
  " `Y8bood8P'  `Y888\"\"8o d888b    `Y8bod88P\" o888o `Y888\"\"8o o88'   888o",
  ""};
*/

void usage(int argc, char** argv)
{
  cout << endl;
  cout << " Usage: monodomain <model> <msh> <dt> <tend> <pr> <outfile>" 
       << endl << endl;
  cout << "\t<model>    string that identifies the ionic model" << endl;
  cout << "\t<mesh>     prefix of the mesh file" << endl;
  cout << "\t<dt>       time step (ms)" << endl;
  cout << "\t<tend>     total time of simulation (ms)" << endl;
  cout << "\t<pr>       print rate to output file" << endl;
  cout << "\t<outfile>  output file (IGB format)" << endl << endl;
}

int main(int argc, char *argv[])
{
  char filename[MAXSTR];
  char outfile[MAXSTR];
  char ionicname[MAXSTR];

  Real dt, ftime, tp;
  string ifile, ofile, cellname;

  if (argc < 7) {
    usage(argc, argv);
    exit(-1);
  }

  sprintf(filename,"%s",argv[1]);
  dt    = atof(argv[2]);
  ftime = atof(argv[3]);
  tp    = atof(argv[4]);
  sprintf(outfile,"%s",argv[5]);
  sprintf(ionicname,"%s",argv[6]);

  ifile = argv[1];
  ofile = argv[5];
  cellname = argv[6];

  //for(int i=0; i<tsize; i++)
  //  printf("%s\n",title[i]);

#ifdef USE_PETSC
  PetscMPIInt rank; // processor rank 
  PetscMPIInt size; // size of communicator

  PetscInitialize(&argc, &argv, (char *) 0, help);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  // PETSc data types
  Mat A, K, M;
  Vec v0, v1;
  KSP ksp;
  PC pc;
  PetscErrorCode ierr;  
#endif

  // MONODOMAIN CODE HERE
  /*
  Monodomain monodomain;
  monodomain.init(argc,argv);
  monodomain.solve();
  */
    
  int ndof, neq, it;
  Real time, *sv;  
  gzFile gzOutFile;

  TimeParameters tip = TimeParameters (dt, ftime, tp);

  Timer timer;
  timer = Timer();

  Mesh msh; 
  msh.read(ifile);

  Stimuli stimuli;
  stimuli.read(ifile);

  ndof = msh.n_points;

  // setup model and cells
  CellModel *cellmodel = CellModel::create(ionicname);
  cellmodel->solver("ForwardEuler");
  cellmodel->setup(dt, ftime);
  
  Cells *cells = new Cells(ndof,cellmodel);
  neq = cells->get_ode_size();

  SparseMatrix Ki, Mi, Ai, invC;
  Vector vm_0(ndof);
  Vector vm_1(ndof);
  Vector vec_b(ndof);

  // open output file
  IGBheader *igbh;
  igbh = create_igb (ofile.c_str(), &gzOutFile, ndof, tip.get_size(), sizeof(Real));

  // -----------------------------------------------------------------
  // setup FEM matrices
  // -----------------------------------------------------------------

  const Real kappa = (SURF2VOL/dt) * UM2_TO_CM2;
  const Real ikapp = 0.5*(1.0/kappa); // Crank Nicolson

  Timing asm_new;
  asm_new.start();
  
  H1FESpace fespace(&msh);
  cout << " Assembling stiffness" << endl;
  fem_assemble_stiffness (msh, fespace ,Ki);
  cout << " Assembling mass" << endl;
  fem_assemble_mass (msh, fespace, Mi);
  fem_setup_matrices2 (Ki, Mi, Ai, ikapp);
  calc_diag_prec (Ai, invC);
  
  asm_new.stop();

  cout << " Assembling total time: " 
       << asm_new.get_total_time() << endl;

#ifdef USE_PETSC
  // PETSC matrices creation
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n_points, n_points); CHKERRQ(ierr);
  ierr = MatSetType(A, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &K); CHKERRQ(ierr);
  ierr = MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, n_points, n_points); CHKERRQ(ierr);
  ierr = MatSetType(K, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetFromOptions(K); CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &M); CHKERRQ(ierr);
  ierr = MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, n_points, n_points); CHKERRQ(ierr);
  ierr = MatSetType(M, MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetFromOptions(M); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, " Assembling Matrix begin (PETSc)\n");
  //ierr = PetscGetFEMatrices(msh, A, K, M); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Assembling Matrix end (PETSc)\n\n");

  // Assemble the matrix
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
#endif
  
  timer.stop_setup();

  // -----------------------------------------------------------------
  // loop for time integration (with operator splitting)
  // -----------------------------------------------------------------

  tip.init_time_loop(); 

  // initial conditions
  cells->init(); 
  cells->get_var(0,vm_1);

  sv = cells->get_state_vars();

  // save results at time=0.0, start time and stimulus
  gzwrite (gzOutFile, vm_1.ref(), ndof*sizeof(Real));
  stimuli.check(tip.it(), tip.time(), neq, sv, msh);

  int cstep = 0;

  // start loop in time
  while( !tip.finished() )
  {
    tip.increase_time();

    it = tip.it();
    time = tip.time();

    if (tip.time2print()) 
      fprintf(stdout, " Time %6.4f  ", time);
    
    //
    // Advance ODEs in time
    //
    timer.start_ode();

    stimuli.check(it, time, neq, sv, msh);
    cells->advance(time, dt);
    cells->get_var(0,vm_0);
    
    timer.stop_ode();

    // Solve the parabolic PDE
    //    b = M * vm0 (sparse matrix vector multiplication)
    //    A * vm1 = b (solve)
    timer.start_parab();

    //cout << "PCG" << time <<  endl;

    vec_b = Mi * vm_0;   
    pcg_solve (Ai, invC, vec_b, vm_1, it, tip.pr(), ndof, 1.0e-8);
    cells->set_var(0, vm_1);

    timer.stop_parab();
     
    if (tip.time2print()) {
      gzwrite (gzOutFile, vm_1.ref(), ndof*sizeof(Real));

      char buf[256];
      sprintf(buf,"vm_%03d.vtu",cstep);
      cstep++;
      string vtunew(buf);
      write_vtu_file_r(msh, vm_1.ref(), vtunew);

      fprintf (stdout, "\n");
      fflush  (stdout);
    }

  }

  // print timings and cleanup memory
  timer.stop();
  timer.print();
  
  gzclose (gzOutFile);

  delete igbh;
  delete cells;
  delete cellmodel;

  // END OF MONODOMAIN

#ifdef USE_PETSC
  MPI_Barrier(PETSC_COMM_WORLD);
  ierr = PetscFinalize(); CHKERRQ(ierr);
#endif

  return 0;
}
