/*
 * Shared main() of every test executable. Linked through cardiax_test_support.
 *
 * What it does around the test cases:
 *   1. Creates a fresh scratch directory for this executable and chdir()s into
 *      it, because solvers delete ./output and write next to their input files.
 *   2. Starts PETSc and CommandLineArgs with the PETSc direct solver selected, so
 *      that results do not depend on a GPU. A test may pass other options on the
 *      command line (see the `gpu` tests in tests/CMakeLists.txt).
 *   3. Redirects stdout of each case into a log file, shown only if it fails.
 *   4. Installs an atexit guard: the library calls exit(0) on some error paths,
 *      which would otherwise look like a passing run to CTest.
 *
 * Exit code: 0 if all cases passed (known bugs still failing count as passing).
 *
 * Options: --list prints the case names; any other argument that does not start
 * with '-' selects the cases whose name contains it.
 */
#include "test_framework.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>

#include "petscsys.h"
#include "petscoptions.h"
#include "util/command_line_args.h"

namespace fs = std::filesystem;

namespace cardiax_test
{

std::vector<TestCase> & registry()
{
  static std::vector<TestCase> cases;
  return cases;
}

std::string source_dir() { return CARDIAX_SOURCE_DIR; }
std::string app_dir()    { return CARDIAX_APP_DIR; }

static std::vector<std::string> g_failures;

void report_failure(const char * file, int line, const std::string & what)
{
  std::ostringstream os;
  os << fs::path(file).filename().string() << ":" << line << ": " << what;
  g_failures.push_back(os.str());
}

} // namespace cardiax_test

using namespace cardiax_test;

static int  g_real_stdout   = -1;
static bool g_main_finished = false;
static std::string g_current_case;

static void say(const std::string & s)
{
  const ssize_t ignored = write(g_real_stdout, s.data(), s.size());
  (void)ignored;
}

//! Fires when something called exit() from inside a case
static void exit_guard()
{
  if (g_main_finished) return;
  say("\n*** FAIL: the code under test called exit() during the case: " + g_current_case + "\n");
  _exit(1);
}

//! Sends stdout and stderr to `path`, so cout, printf, HDF5 and PETSc messages are
//! all captured
static void capture_stdout(const std::string & path)
{
  fflush(stdout);
  fflush(stderr);
  std::cout.flush();
  std::cerr.flush();
  const int fd = open(path.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
  if (fd >= 0) { dup2(fd, 1); dup2(fd, 2); close(fd); }
}

static void print_log_tail(const std::string & path, std::size_t max_lines)
{
  std::ifstream in(path);
  std::vector<std::string> lines;
  std::string l;
  while (std::getline(in, l)) lines.push_back(l);
  const std::size_t first = lines.size() > max_lines ? lines.size() - max_lines : 0;
  say("    --- last lines of the case output (" + path + ") ---\n");
  for (std::size_t i = first; i < lines.size(); i++) say("    | " + lines[i] + "\n");
}

int main(int argc, char ** argv)
{
  // --- arguments ------------------------------------------------------------
  std::vector<std::string> filters;
  bool list_only = false;

  std::vector<const char *> args;
  args.push_back(argv[0]);
  // The default for tests: PETSc direct solver, no GPU
  static const char * defaults[] = {"-linsolver", "petsc"};
  bool user_picked_backend = false;
  for (int i = 1; i < argc; i++)
    if (!std::strcmp(argv[i], "-linsolver")) user_picked_backend = true;
  if (!user_picked_backend)
    for (const char * d : defaults) args.push_back(d);

  for (int i = 1; i < argc; i++)
  {
    if (!std::strcmp(argv[i], "--list")) { list_only = true; continue; }
    if (argv[i][0] != '-' && (i == 1 || argv[i-1][0] != '-'))
    { filters.push_back(argv[i]); continue; }
    args.push_back(argv[i]);
  }
  args.push_back(nullptr);

  if (list_only)
  {
    for (const TestCase & c : registry()) std::cout << c.name << "\n";
    return 0;
  }

  // --- scratch directory ----------------------------------------------------
  const std::string exe_name = fs::path(argv[0]).filename().string();
  const fs::path scratch = fs::path(CARDIAX_SCRATCH_ROOT) / exe_name;
  fs::remove_all(scratch);
  fs::create_directories(scratch);
  fs::current_path(scratch);

  // --- stdout bookkeeping and exit guard -------------------------------------
  g_real_stdout = dup(1);
  std::atexit(exit_guard);

  // --- PETSc ----------------------------------------------------------------
  int    pargc = static_cast<int>(args.size()) - 1;
  char ** pargv = const_cast<char **>(args.data());
  CommandLineArgs::init(pargc, args.data());

  capture_stdout("petsc_init.log");
  PetscErrorCode ierr = PetscInitialize(&pargc, &pargv, nullptr, "cardiax test");
  if (ierr) { say("PetscInitialize failed\n"); _exit(2); }
  if (!user_picked_backend)
  {
    PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-pc_type", "lu");
  }

  // --- the cases ------------------------------------------------------------
  int n_pass = 0, n_fail = 0, n_skip = 0, n_known = 0;
  say("[" + exe_name + "] scratch: " + scratch.string() + "\n");

  int index = 0;
  for (const TestCase & tc : registry())
  {
    index++;
    if (!filters.empty() &&
        std::none_of(filters.begin(), filters.end(),
                     [&](const std::string & f)
                     { return tc.name.find(f) != std::string::npos; }))
      continue;

    const std::string log = "case_" + std::to_string(index) + ".log";
    capture_stdout(log);
    g_failures.clear();
    g_current_case = tc.name;

    bool skipped = false;
    std::string skip_reason;
    try
    {
      tc.fn();
    }
    catch (const AbortCase &) {}
    catch (const SkipCase & s) { skipped = true; skip_reason = s.why; }
    catch (const std::exception & e)
    {
      g_failures.push_back(std::string("uncaught exception: ") + e.what());
    }
    catch (const char * e)
    {
      g_failures.push_back(std::string("uncaught exception: ") + e);
    }
    catch (...)
    {
      g_failures.push_back("uncaught exception of unknown type");
    }
    fflush(stdout);
    std::cout.flush();

    const bool failed = !g_failures.empty();

    if (skipped)
    {
      n_skip++;
      say("  SKIP       " + tc.name + " (" + skip_reason + ")\n");
    }
    else if (!tc.known_bug.empty())
    {
      if (failed)
      {
        n_known++;
        say("  KNOWN BUG  " + tc.name + "\n             " + tc.known_bug + "\n");
      }
      else
      {
        n_fail++;
        say("  FAIL       " + tc.name + "\n             this case is marked as a known "
            "bug but it passes now: replace TEST_CASE_KNOWN_BUG by TEST_CASE.\n");
      }
    }
    else if (failed)
    {
      n_fail++;
      say("  FAIL       " + tc.name + "\n");
      for (const std::string & f : g_failures) say("             " + f + "\n");
      print_log_tail(log, 25);
    }
    else
    {
      n_pass++;
      say("  ok         " + tc.name + "\n");
    }
  }

  capture_stdout("petsc_finalize.log");
  PetscFinalize();

  std::ostringstream summary;
  summary << "[" << exe_name << "] " << n_pass << " passed, " << n_fail << " failed, "
          << n_known << " known bugs, " << n_skip << " skipped\n";
  say(summary.str());

  g_main_finished = true;

  // Exit straight away. Some solvers keep state in function-local statics that
  // are still alive at this point, and their destructors would run after PETSc
  // was finalized.
  _exit(n_fail > 0 ? 1 : 0);
}
