/*
 * Minimal test framework for Cardiax.
 *
 * Why not GoogleTest/Catch2: it would be one more dependency for install_deps.sh,
 * and the suite needs things they do not give out of the box:
 *   - PETSc must be initialised once per process, around all test cases;
 *   - the library calls exit(0) on some error paths, which must count as a
 *     FAILURE and not as a successful run (see the atexit guard in test_main.cpp);
 *   - the library prints a lot to stdout, so each case's output is captured to a
 *     log file and only shown if the case fails.
 *
 * Usage:
 *   TEST_CASE("description") { CHECK(a == b); CHECK_NEAR(x, 1.0, 1e-12); }
 *   TEST_CASE_KNOWN_BUG("description", "why") { ... }   // see below
 *
 * Assertions:
 *   CHECK(cond)              records a failure and continues
 *   CHECK_EQ(a, b)           same, prints both values
 *   CHECK_NEAR(a, b, tol)    |a - b| <= tol
 *   CHECK_REL(a, b, rtol)    |a - b| <= rtol * max(|a|, |b|, tiny)
 *   CHECK_THROWS(expr)       expr must throw
 *   REQUIRE(cond)            records a failure and aborts the current case
 *   SKIP_TEST(reason)        stops the current case and reports it as skipped
 *
 * TEST_CASE_KNOWN_BUG documents a defect that exists in the code today. It runs
 * like any other case and must keep failing: the run stays green and prints
 * KNOWN BUG. Once the code is fixed the case starts passing, which is reported
 * as an error so that the marker gets removed and the case becomes a regular
 * regression test.
 */
#ifndef CARDIAX_TEST_FRAMEWORK_HPP
#define CARDIAX_TEST_FRAMEWORK_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cardiax_test
{

struct TestCase
{
  std::string name;
  void (*fn)();
  std::string known_bug;  //!< empty for a regular case
};

std::vector<TestCase> & registry();

struct Registrar
{
  Registrar(const char * name, void (*fn)(), const char * known_bug = "")
  {
    registry().push_back({name, fn, known_bug});
  }
};

struct AbortCase {};               //!< thrown by REQUIRE
struct SkipCase { std::string why; };  //!< thrown by SKIP_TEST

//! Records a failure of the running case
void report_failure(const char * file, int line, const std::string & what);

//! Path of the repository (for fixtures in testes/, examples/, configs/)
std::string source_dir();

//! Path of the directory with the solver executables
std::string app_dir();

template <typename T>
std::string show(const T & v)
{
  std::ostringstream os;
  os.precision(16);
  os << v;
  return os.str();
}

inline double rel_gap(double a, double b)
{
  const double scale = std::max({std::fabs(a), std::fabs(b),
                                 std::numeric_limits<double>::min()});
  return std::fabs(a - b) / scale;
}

} // namespace cardiax_test

#define CXT_CONCAT2(a, b) a##b
#define CXT_CONCAT(a, b) CXT_CONCAT2(a, b)

#define TEST_CASE(name)                                                        \
  static void CXT_CONCAT(cxt_case_, __LINE__)();                               \
  static cardiax_test::Registrar CXT_CONCAT(cxt_reg_, __LINE__)(               \
      name, &CXT_CONCAT(cxt_case_, __LINE__));                                 \
  static void CXT_CONCAT(cxt_case_, __LINE__)()

#define TEST_CASE_KNOWN_BUG(name, why)                                         \
  static void CXT_CONCAT(cxt_case_, __LINE__)();                               \
  static cardiax_test::Registrar CXT_CONCAT(cxt_reg_, __LINE__)(               \
      name, &CXT_CONCAT(cxt_case_, __LINE__), why);                            \
  static void CXT_CONCAT(cxt_case_, __LINE__)()

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond))                                                               \
      cardiax_test::report_failure(__FILE__, __LINE__, "CHECK(" #cond ")");    \
  } while (0)

#define CHECK_MSG(cond, msg)                                                   \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::ostringstream cxt_os;                                               \
      cxt_os << "CHECK(" #cond "): " << msg;                                   \
      cardiax_test::report_failure(__FILE__, __LINE__, cxt_os.str());          \
    }                                                                          \
  } while (0)

#define CHECK_EQ(a, b)                                                         \
  do {                                                                         \
    const auto cxt_a = (a);                                                    \
    const auto cxt_b = (b);                                                    \
    if (!(cxt_a == cxt_b))                                                     \
      cardiax_test::report_failure(__FILE__, __LINE__,                         \
          std::string("CHECK_EQ(" #a ", " #b "): ") +                          \
          cardiax_test::show(cxt_a) + " != " + cardiax_test::show(cxt_b));     \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                  \
  do {                                                                         \
    const double cxt_a = (a), cxt_b = (b), cxt_t = (tol);                      \
    if (!(std::fabs(cxt_a - cxt_b) <= cxt_t))                                  \
      cardiax_test::report_failure(__FILE__, __LINE__,                         \
          std::string("CHECK_NEAR(" #a ", " #b "): ") +                        \
          cardiax_test::show(cxt_a) + " vs " + cardiax_test::show(cxt_b) +     \
          " (tol " + cardiax_test::show(cxt_t) + ")");                         \
  } while (0)

#define CHECK_REL(a, b, rtol)                                                  \
  do {                                                                         \
    const double cxt_a = (a), cxt_b = (b), cxt_t = (rtol);                     \
    if (!(cardiax_test::rel_gap(cxt_a, cxt_b) <= cxt_t))                       \
      cardiax_test::report_failure(__FILE__, __LINE__,                         \
          std::string("CHECK_REL(" #a ", " #b "): ") +                         \
          cardiax_test::show(cxt_a) + " vs " + cardiax_test::show(cxt_b) +     \
          " (rel gap " + cardiax_test::show(cardiax_test::rel_gap(cxt_a, cxt_b)) + \
          ", tol " + cardiax_test::show(cxt_t) + ")");                         \
  } while (0)

#define CHECK_THROWS(expr)                                                     \
  do {                                                                         \
    bool cxt_threw = false;                                                    \
    try { expr; } catch (...) { cxt_threw = true; }                            \
    if (!cxt_threw)                                                            \
      cardiax_test::report_failure(__FILE__, __LINE__,                         \
                                   "CHECK_THROWS(" #expr "): did not throw");  \
  } while (0)

#define REQUIRE(cond)                                                          \
  do {                                                                         \
    if (!(cond)) {                                                             \
      cardiax_test::report_failure(__FILE__, __LINE__, "REQUIRE(" #cond ")");  \
      throw cardiax_test::AbortCase{};                                         \
    }                                                                          \
  } while (0)

#define SKIP_TEST(reason) throw cardiax_test::SkipCase{reason}

#endif
