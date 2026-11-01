/*
 * Helpers shared by the tests: small generated meshes (written as the XML format
 * the solvers read) and a way to run the solver executables and read back what
 * they wrote.
 *
 * Meshes are generated rather than stored because the problem they need (size,
 * markers, fibres, loads) is part of what each test is about, and a test reads
 * better when the input is visible next to the assertions.
 *
 * Generated box meshes are [0,Lx]x[0,Ly]x[0,Lz] with these boundary markers:
 *   x=0: 1   x=Lx: 2   y=0: 3   y=Ly: 4   z=0: 5   z=Lz: 6
 */
#ifndef CARDIAX_TEST_FIXTURES_HPP
#define CARDIAX_TEST_FIXTURES_HPP

#include <map>
#include <string>
#include <vector>
#include <armadillo>

namespace cardiax_test
{

struct BoxSpec
{
  int nx = 2, ny = 2, nz = 2;
  double Lx = 1.0, Ly = 1.0, Lz = 1.0;
  arma::vec3 fiber = {1.0, 0.0, 0.0};
  arma::vec3 sheet = {0.0, 1.0, 0.0};
  std::string extra_mesh_xml;  //!< inserted inside <mesh> (e.g. <pvloop>)
  std::string prefix_xml;      //!< inserted before <mesh> (e.g. <elasticity>)
};

//! Node index in a generated box mesh
inline int box_node(const BoxSpec & b, int i, int j, int k)
{ return i + j * (b.nx + 1) + k * (b.nx + 1) * (b.ny + 1); }

//! Number of nodes of a generated box mesh
inline int box_num_nodes(const BoxSpec & b)
{ return (b.nx + 1) * (b.ny + 1) * (b.nz + 1); }

//! Node coordinates of a generated box mesh
inline arma::vec3 box_point(const BoxSpec & b, int i, int j, int k)
{
  return arma::vec3({i * b.Lx / b.nx, j * b.Ly / b.ny, k * b.Lz / b.nz});
}

//! Conforming tetrahedral box (6 tets per cell), triangles on the boundary
void write_tet_box(const std::string & path, const BoxSpec & spec);

//! Hexahedral box, quadrilaterals on the boundary
void write_hex_box(const std::string & path, const BoxSpec & spec);

//! Triangle mesh of [0,Lx]x[0,Ly] (nz, Lz ignored); boundary markers 1..4 as above
void write_tri_square(const std::string & path, const BoxSpec & spec);

//! Quadrilateral mesh of [0,Lx]x[0,Ly] (nz, Lz ignored)
void write_quad_square(const std::string & path, const BoxSpec & spec);

//! Line mesh with n elements along x, length L; boundary markers 1 and 2
void write_cable(const std::string & path, int n, double L);

//! <prescribed_displacement> entries fixing the given directions on a set of nodes
std::string fix_nodes_xml(const std::vector<int> & nodes, const std::vector<int> & dirs);

//! Writes a text file
void write_text(const std::string & path, const std::string & contents);

//! Reads a whole text file ("" if it does not exist)
std::string read_text(const std::string & path);

// ---- running the solvers ---------------------------------------------------

struct AppResult
{
  int exit_code = -1;
  std::string log;    //!< stdout + stderr of the run
  std::string dir;    //!< directory where the run happened
};

//! Runs build/app/<app> in <dir> (created if needed), with stdout+stderr captured.
//! `args` go to the shell as they are.
AppResult run_app(const std::string & app, const std::string & args,
                  const std::string & dir = ".");

//! Path of an input in the repository, e.g. repo_file("configs/CG_DILU.json")
std::string repo_file(const std::string & relative);

//! Copies a repository input into the current directory and returns its local name.
//! Solvers write next to their input, so tests never run on the original file.
std::string copy_input(const std::string & relative, const std::string & dir = ".");

//! First number after the LAST occurrence of `key` in `text`. NaN if absent.
double number_after(const std::string & text, const std::string & key);

//! All numbers following each occurrence of `key` in `text`
std::vector<double> numbers_after(const std::string & text, const std::string & key);

//! Numeric table of a whitespace-separated text file. Lines starting with # or
//! that do not begin with a number are skipped.
arma::mat read_table(const std::string & path);

//! True if an HDF5 file contains the dataset
bool h5_has(const std::string & file, const std::string & dataset);

//! Reads an HDF5 dataset as doubles (rows x cols; a 1-D dataset gives a column)
arma::mat h5_read(const std::string & file, const std::string & dataset);

//! Names of all datasets in an HDF5 file, with full paths
std::vector<std::string> h5_datasets(const std::string & file);

//! True if two HDF5 files hold the same datasets, within a relative tolerance.
//! Datasets whose names are in `ignore` are not compared.
bool h5_equal(const std::string & a, const std::string & b, double rtol,
              std::string * why = nullptr,
              const std::vector<std::string> & ignore = {});

//! Regression baseline: compares `value` to the entry `key` of
//! tests/golden/baseline.txt. Setting the environment variable
//! CARDIAX_UPDATE_GOLDEN=1 rewrites the entry instead. The file holds values the
//! code produced when the baseline was taken; they catch unintended changes, and
//! do not prove the value is physically right.
void check_golden(const char * file, int line, const std::string & key,
                  double value, double rtol);

#define CHECK_GOLDEN(key, value, rtol)                                         \
  cardiax_test::check_golden(__FILE__, __LINE__, key, value, rtol)

} // namespace cardiax_test

#endif
