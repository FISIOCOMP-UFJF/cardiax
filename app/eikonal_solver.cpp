#include <iostream>
#include <string>
#include "util/toml.hpp"
#include "util/util.hpp"
#include "util/command_line_args.h"
#include "pdes/eikonal.hpp"
#include "mesh/mesh.hpp"
#include <filesystem>

static char help[] = "Eikonal Pure PDE Solver.\n\n";

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

static int conductivity_from_string(const std::string & s)
{
  if (s == "S_ISOTROPIC")   return S_ISOTROPIC;
  if (s == "S_TRANSVERSE")  return S_TRANSVERSE;
  if (s == "S_ORTHOTROPIC") return S_ORTHOTROPIC;
  if (s == "M_ISOTROPIC")   return M_ISOTROPIC;
  if (s == "M_TRANSVERSE")  return M_TRANSVERSE;
  if (s == "M_ORTHOTROPIC") return M_ORTHOTROPIC;

  cout << " *** WARNING: unknown conductivity_type '" << s
       << "', using M_TRANSVERSE." << endl;
  return M_TRANSVERSE;
}

int main(int argc, const char *argv[])
{
  string mshname, meshbase;

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

  // [physical]
  string cond_str = cfg["physical"]["conductivity_type"].value_or("M_TRANSVERSE");

  // Start PETSc
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, (char ***)&argv, (char *) 0, help); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  print_banner();
  msg(("Reading config file: " + mshname).c_str());

  // 1. Instancia e carrega a malha baseada no XML
  Mesh* mesh = new Mesh();
  mesh->read_xml(meshbase);

  // 2. Cria o solver puro
  Eikonal* eikonal = new Eikonal();

  // 3. Aplica propriedades físicas e resolve a propagação
  eikonal->set_conductivity(conductivity_from_string(cond_str));
  eikonal->solve(mesh, meshbase);

  msg("Done.");
  ierr = PetscFinalize(); CHKERRQ(ierr);

  // Limpeza de memória
  delete eikonal;
  delete mesh;

  return 0;
}