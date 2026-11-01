/*
 * Tests of the AMGX backend of petsc::LinearSolver, on the GPU (label `gpu`).
 *
 * What is covered and why:
 *  - Backend selection: `-linsolver petsc` and `-linsolver amgx` pick what they say,
 *    and with no flag a PETSc option selects PETSc, as CLAUDE.md documents. Selection is
 *    decided once per process, so this is checked through the executables.
 *  - AMGX and the direct PETSc solver reach the same answer on a Poisson problem, and
 *    on one mechanics solve, whose matrix is the largest and worst conditioned path.
 *  - A missing AMGX config falls back to the built-in one.
 *
 * The cases are skipped when the build has no AMGX or the machine has no NVIDIA GPU.
 * Problems are small on purpose: this exercises the wrapper, it is not a benchmark.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <cstdlib>
#include <filesystem>

using namespace cardiax_test;
namespace fs = std::filesystem;

namespace
{

void require_gpu()
{
#ifndef AMGX_SOLVER
  SKIP_TEST("built without AMGX");
#endif
  if (std::system("nvidia-smi -L > /dev/null 2>&1") != 0)
    SKIP_TEST("no NVIDIA GPU available");
}

void write_poisson(const std::string & dir)
{
  BoxSpec b;
  b.nx = 8; b.ny = b.nz = 4;
  fs::create_directories(dir);
  write_tet_box(dir + "/m.xml", b);
  write_text(dir + "/c.toml",
    "[problem]\nmesh = \"m.xml\"\n"
    "[physical]\nsource = \"2\"\n"
    "dirichlet = [ { marker = 1, value = 0.0 }, { marker = 2, value = 0.0 } ]\n"
    "neumann = [ { marker = 3, value = 0.0 }, { marker = 4, value = 0.0 }, { marker = 5, value = 0.0 }, { marker = 6, value = 0.0 } ]\n"
    "[output]\nfilename = \"out\"\n"
    "[validation]\nexact = \"x*(1-x)\"\n");
}

std::string amgx_flags()
{
  return "-linsolver amgx -amgx " + repo_file("configs/CG_DILU.json");
}

} // namespace

TEST_CASE("Backend selection: the flag decides, and PETSc options imply PETSc")
{
  write_poisson("sel_petsc");
  const AppResult p = run_app("poisson", "-f c.toml -linsolver petsc -ksp_type preonly -pc_type lu", "sel_petsc");
  CHECK_EQ(p.exit_code, 0);
  CHECK(p.log.find("backend: PETSc KSP") != std::string::npos);

  write_poisson("sel_opt");
  const AppResult o = run_app("poisson", "-f c.toml -ksp_type preonly -pc_type lu", "sel_opt");
  CHECK_EQ(o.exit_code, 0);
  CHECK(o.log.find("backend: PETSc KSP (opcao PETSc") != std::string::npos);

  require_gpu();
  write_poisson("sel_amgx");
  const AppResult a = run_app("poisson", "-f c.toml " + amgx_flags(), "sel_amgx");
  CHECK_EQ(a.exit_code, 0);
  CHECK(a.log.find("backend: AMGX") != std::string::npos);
}

TEST_CASE("AMGX: Poisson gives the same solution as the direct solver")
{
  require_gpu();
  write_poisson("amgx_p");
  write_poisson("direct_p");
  const AppResult a = run_app("poisson", "-f c.toml " + amgx_flags(), "amgx_p");
  const AppResult d = run_app("poisson", "-f c.toml -linsolver petsc -ksp_type preonly -pc_type lu", "direct_p");
  REQUIRE(a.exit_code == 0);
  REQUIRE(d.exit_code == 0);
  const arma::mat ua = h5_read("amgx_p/out.h5", "/point_data/vm");
  const arma::mat ud = h5_read("direct_p/out.h5", "/point_data/vm");
  CHECK_MSG(arma::abs(ua - ud).max() < 1e-3, "max gap " << arma::abs(ua - ud).max());
}

TEST_CASE("AMGX: a mechanics solve converges to the state of the direct solver")
{
  require_gpu();
  const char * input = "testes/passive/teste_unitario_guccione.xml";
  const std::string n1 = copy_input(input, "amgx_m");
  const std::string n2 = copy_input(input, "direct_m");
  const AppResult a = run_app("nonlinearelas", "-f " + n1 + " -num_threads 2 " + amgx_flags(), "amgx_m");
  const AppResult d = run_app("nonlinearelas",
      "-f " + n2 + " -num_threads 2 -linsolver petsc -ksp_type preonly -pc_type lu", "direct_m");
  REQUIRE(a.exit_code == 0);
  REQUIRE(d.exit_code == 0);
  const std::string base = n1.substr(0, n1.rfind(".xml"));
  const double va = number_after(a.log, "End inner volume:");
  const double vd = number_after(d.log, "End inner volume:");
  CHECK_REL(va, vd, 1e-4);
  CHECK(fs::exists("amgx_m/" + base + "_output_nl.h5"));
}

TEST_CASE("AMGX: a missing config file falls back to the built-in solver")
{
  require_gpu();
  write_poisson("amgx_fb");
  const AppResult a = run_app("poisson", "-f c.toml -linsolver amgx -amgx does_not_exist.json", "amgx_fb");
  CHECK_EQ(a.exit_code, 0);
  CHECK_MSG(number_after(a.log, "L2 error norm:") < 5e-3, number_after(a.log, "L2 error norm:"));
}

TEST_CASE("-linsolver amgx in a build without AMGX falls back with a message")
{
#ifdef AMGX_SOLVER
  SKIP_TEST("this build has AMGX");
#endif
  write_poisson("no_amgx");
  const AppResult a = run_app("poisson", "-f c.toml -linsolver amgx", "no_amgx");
  CHECK(a.log.find("compilado sem AMGX") != std::string::npos);
}
