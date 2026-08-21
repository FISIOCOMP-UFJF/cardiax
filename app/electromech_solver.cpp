#include "util/log.hpp"
#include "util/command_line_args.h"
#include <filesystem>
//#include "cardiac/coupled_electromechanic.hpp"
#include "cardiac/cardiac_cycle.hpp"
#include "../src/util/command_line_args.h"

using namespace std;
namespace fs = std::filesystem;


static char help[] = "coupled monodomain-mechanics solver";

void usage()
{
  cout << endl;
  cout << " Usage: monomech [OPTIONS]" << endl << endl;
  cout << "    -f    <meshbase>   prefix of the mesh file" << endl;
  cout << "    -dt   <dt>         time step (ms)" << endl;
  cout << "    -t    <time>       with -circ 1: MAXIMUM simulation time; the" << endl;
  cout << "                       run is -beats * 1/hr and -t only caps it." << endl;
  cout << "                       Without -circ: total time of simulation." << endl;
  cout << "    -c    <cell>       name of the cell model (FHN,NP,LR1,TT2)" << endl; 
  cout << "    -m    <odesolver>  ODE solver (ExplicitEuler, ImplicitEuler)" << endl;
  cout << "    -cond <condtype>   conductivity type (0 1 ... 5)" << endl;
  cout << "    -ep   <EP model>   Monodomain or Bidomain (mono or bido)" << endl;
  cout << "    -circ <0|1>        closed-loop 0D circulation (Regazzoni 2022)" << endl;
  cout << "    -beats <n>         number of beats with -circ 1" << endl;
  cout << "    -hr   <Hz>         heart rate (default 1.0 Hz = 60 bpm)" << endl;
  cout << "    -bpm  <bpm>        heart rate in beats per minute" << endl;
  cout << "    -bcl  <time>       warm-up BCL (default: the heart period)" << endl;
  cout << "    -tact <s>          ventricular activation instant (0)" << endl;
  cout << "    -emd  <s>          electromechanical delay, stimulus->Ta (0.04)" << endl;
  cout << "    -tascale <x>       GLOBAL contractility scale (default 1.0);" << endl;
  cout << "                       multiplies Ta at every node. The per-material" << endl;
  cout << "                       ta_scale in the XML is a regional fraction" << endl;
  cout << "                       (1 = myocardium, 0 = valve plug/scar) and" << endl;
  cout << "                       multiplies on top of this one." << endl;
  cout << "    -circ_vtot <mL>    stressed blood volume (default: from the" << endl;
  cout << "                       0D initial state, 1617.9 mL)" << endl;
  cout << "    -amgx <file>       AMGX config (./configs/CG_DILU.json)" << endl;
  cout << "    -dumpflags <0|1>   print every flag read, with [cli]/[default]" << endl;
  cout << "                       markers, at the end of config() (default 1)" << endl;
  cout << "    -dtc  <dt>         circulation time step in s (default 1e-3)" << endl;
  cout << "    -prc  <n>          save the mesh every n circulation steps" << endl;
  cout << "    -logc <n>          detailed log every n circulation steps (25)" << endl;
  cout << "    -c    <model>      cell model (Kerkoff2003 or ToRORdLand)" << endl;
  cout << "    -tunit <s|ms>      EP time unit (default s)" << endl;
  cout << "                       ms: use -t 800 -dt 0.02 for ToRORdLand" << endl;
  cout << "                       s : use -t 0.8 -dt 0.001 for Kerkoff2003" << endl;
  cout << "    -stimamp <A>       LAT stimulus amplitude, ToRORdLand (-53)" << endl;
  cout << "    -stimdur <d>       LAT stimulus duration in solver units (0.002)" << endl;
  cout << "    -kpa2p <f>         model pressure units per kPa" << endl;
  cout << "                       1000 if material is in Pa (default)," << endl;
  cout << "                       1 if in kPa. Also sets the mmHg" << endl;
  cout << "                       conversion used by the circulation." << endl;
  cout << endl;
  cout << " Nem toda flag lida esta listada acima. Para o inventario completo" << endl;
  cout << " de uma rodada, com valor resolvido e marca [cli]/[default], veja o" << endl;
  cout << " bloco \"flags desta rodada\" impresso no fim da configuracao." << endl;
  cout << endl;
  exit(0);
}

int main(int argc, const char* argv[])
{
  string basename, meshname, ep_model;
  int condtype, use_circ, num_beats;
  double dt_circ;

  if (argc <= 1) usage();

  CommandLineArgs::init(argc, argv);
  if (CommandLineArgs::wants_help()) usage();
  basename = CommandLineArgs::read("-f","emptymesh");
  ep_model = CommandLineArgs::read("-ep","mono");
  condtype = CommandLineArgs::read("-cond",0);
  use_circ = CommandLineArgs::read("-circ",0);
  num_beats = CommandLineArgs::read("-beats",1);
  dt_circ  = CommandLineArgs::read("-dtc",1.0e-3);

  // check and clean output directory
 if ( fs::exists("output") ) 
  {
    fs::remove_all("output");
    assert(!fs::exists("output"));
  }
  
  fs::create_directory("output");
  fs::create_directories("output/vm_text/");

  // start solution
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  // start PDE solver
  {
    CardiacElectromechanic model(ep_model);
    model.config(basename);
    cout << "Config OK" << endl;
    model.ref().set_conductivity(condtype);
    if (use_circ)
    {
      cout << "Iniciando solve com circulacao 0D fechada" << endl;
      model.solve_circulation(num_beats, dt_circ);
    }
    else
    {
      cout << "Iniciando solve"<<endl;
      model.solve();
    }
  }
  // end of PDE solver

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}