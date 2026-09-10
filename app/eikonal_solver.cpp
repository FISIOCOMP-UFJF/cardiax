#include <iostream>
#include <string>
#include "util/toml.hpp"
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "cardiac/eikonal.hpp"
#include "cardiac/coupled_purkinje.hpp"
#include <filesystem>

static char help[] = "Eikonal Cardiac Solver.\n\n";

using namespace std;
namespace fs = std::filesystem;

void usage()
{
  cout << endl;
  cout << " Usage: ";
  cout << "Eikonal [OPTIONS]" << endl << endl;
  cout << "    -f     <config.toml>          TOML input file" << endl;
  cout << endl;
  exit(0);
}

void print_banner()
{
  cout << "========================================" << endl;
  cout << "     E I K O N A L   S O L V E R        " << endl;
  cout << "========================================" << endl;
}

// Map the conductivity_type string from the TOML onto the Eikonal enum.
// Falls back to M_TRANSVERSE (with a warning) on an unknown name.
static int conductivity_from_string(const std::string & s)
{
  if (s == "S_ISOTROPIC")   return Eikonal::S_ISOTROPIC;
  if (s == "S_TRANSVERSE")  return Eikonal::S_TRANSVERSE;
  if (s == "S_ORTHOTROPIC") return Eikonal::S_ORTHOTROPIC;
  if (s == "M_ISOTROPIC")   return Eikonal::M_ISOTROPIC;
  if (s == "M_TRANSVERSE")  return Eikonal::M_TRANSVERSE;
  if (s == "M_ORTHOTROPIC") return Eikonal::M_ORTHOTROPIC;

  cout << " *** WARNING: unknown conductivity_type '" << s
       << "', using M_TRANSVERSE." << endl;
  return Eikonal::M_TRANSVERSE;
}

int main(int argc, const char *argv[])
{
  string mshname, meshbase, cellmodel, odesolver;

  // Parse command line options: -f is the TOML input file.
  CommandLineArgs::init(argc, argv);
  mshname = CommandLineArgs::read("-f", "config.toml");

  // ---- read the TOML config -------------------------------------------
  toml::table cfg;
  try
  {
    cfg = toml::parse_file(mshname);
  }
  catch (const toml::parse_error & err)
  {
    cerr << " *** ERROR parsing '" << mshname << "':\n" << err << endl;
    return 1;
  }

  // [problem]
  meshbase = cfg["problem"]["mesh"].value_or("emptymesh");

  // [numerics]  ->  setup(dt, T, print_rate, print_rate_apd)
  double dt = cfg["numerics"]["timestep"].value_or(0.1);
  double T  = cfg["numerics"]["total_time"].value_or(10.0);
  double pr = cfg["numerics"]["print_rate"].value_or(1.0);
  double pa = cfg["numerics"]["print_rate_apd"].value_or(1.0);

  // [cell]  ->  setup(c, m)
  cellmodel = cfg["cell"]["model"].value_or("Kerkoff2003");
  odesolver = cfg["cell"]["ode_solver"].value_or("ExplicitEuler");

  // [physical]
  string cond_str = cfg["physical"]["conductivity_type"].value_or("M_TRANSVERSE");

  // [stimulus.lat]  ->  apply_lat_stimulus(amplitude, duration, period)
  bool   lat_stim_on = cfg["stimulus"]["lat"]["enabled"].value_or(false);
  double lat_amp     = cfg["stimulus"]["lat"]["amplitude"].value_or(1.0);
  double lat_dur     = cfg["stimulus"]["lat"]["duration"].value_or(0.002);
  double lat_per     = cfg["stimulus"]["lat"]["period"].value_or(0.0);

  // typefile = meshbase + ".typ";

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  Eikonal* eikonal = new Eikonal();

  print_banner();
  msg(("Reading config file: " + mshname).c_str());

  // setup() takes non-const std::string references, so pass named variables.
  eikonal->setup(meshbase, cellmodel, odesolver, dt, T, pr, pa);
  eikonal->set_conductivity(conductivity_from_string(cond_str));
  eikonal->init();
  eikonal->solve(meshbase);

  // // LAT-driven stimulus for ionic cell models (no-op for phenomenological
  // // Kerckhoffs). Only meaningful once the LAT has been computed by solve().
  // if (lat_stim_on)
  //   eikonal->apply_lat_stimulus(lat_amp, lat_dur, lat_per);

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}

// OLD CODE

// #include <iostream>
// #include <string>
// #include "util/util.hpp"
// #include "util/command_line_args.h"
// #include "cardiac/eikonal.hpp"
// #include "cardiac/coupled_purkinje.hpp"
// #include<filesystem>

// static char help[] = "Eikonal Cardiac Solver.\n\n";

// using namespace std;
// namespace fs = std::filesystem;

// void usage()
// {
//   cout << endl;
//   cout << " Usage: ";
//   cout << "Eikonal [OPTIONS]" << endl << endl;
//   cout << "    -f     <meshbase>              prefix of the mesh file" << endl;
//   cout << endl;
//   exit(0);
// }

// void print_banner()
// {
//   cout << "========================================" << endl;
//   cout << "     E I K O N A L   S O L V E R        " << endl;
//   cout << "========================================" << endl;
// }

// int main(int argc, const char *argv[])
// {
//   double dt = 0.1, T = 10.0, tp = 1.0; 
//   string mshname, pkmshname, cellmodel, odesolver, typefile, restfilename;

//   // Parse command line options
//   CommandLineArgs::init(argc, argv);
//   T  = CommandLineArgs::read("-t",10.0);
//   mshname = CommandLineArgs::read("-f", "emptymesh");
//   cellmodel = CommandLineArgs::read("-c","Kerkoff2003"); 
//   odesolver = CommandLineArgs::read("-m","ExplicitEuler");

//   typefile  = mshname + ".typ";

//   // Start PETSc
//   PetscMPIInt rank;
//   PetscMPIInt size;
//   PetscErrorCode ierr;

//   ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
//   ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

//   Eikonal* eikonal = new Eikonal();

//   print_banner();
//   msg("Reading parameters file");

//   eikonal->setup(mshname, cellmodel, odesolver, dt, T, tp, tp);
//   eikonal->init();
//   eikonal->solve(mshname);   

//   msg("Done.");
//   ierr = PetscFinalize(); CHKERRQ(ierr);
  
//   return 0;
// }