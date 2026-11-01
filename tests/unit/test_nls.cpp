/*
 * Unit tests of src/nls: the Newton solvers on small problems with known roots.
 *
 * What is covered and why:
 *  - Newton on the 2x2 problem of Kelley ("Solving nonlinear equations with Newton's
 *    method"): converges to the root (1, 1), with the quadratic rate that identifies a
 *    correct Newton step.
 *  - Newton on a problem it cannot solve: it must stop and throw after its maximum
 *    number of iterations, and not loop forever or return as if it had converged.
 *  - NewtonLineSearch: reaches the same root; on a problem where the full Newton step
 *    overshoots (arctan), the line search is what makes it converge.
 *
 * Not covered: BFGSSolver. It depends on PardisoSolver and the Pardiso library, which
 * is not part of the environment built by install_deps.sh.
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include <cmath>
#include <vector>

#include "nls/newton.hpp"
#include "nls/newton_ls.hpp"
#include "util/util.hpp"

using namespace cardiax_test;

namespace
{

//! F(x) = (x0^2 + x1^2 - 2, exp(x0 - 1) + x1^2 - 2); root at (1, 1)
class Kelley : public NonlinearProblem
{
public:
  double x[2] = {2.0, 0.5};
  std::vector<double> history;

  void evaluate(petsc::Vector & r) override
  {
    r.set(0, x[0]*x[0] + x[1]*x[1] - 2.0);
    r.set(1, std::exp(x[0]-1.0) + x[1]*x[1] - 2.0);
    r.assemble();
  }
  void jacobian(petsc::Matrix & K) override
  {
    K.set(0, 0, 2*x[0]);          K.set(0, 1, 2*x[1]);
    K.set(1, 0, std::exp(x[0]-1.0)); K.set(1, 1, 2*x[1]);
    K.assemble();
  }
  void update(petsc::Vector & u, double s) override
  {
    double dx[2];
    u.get_data(dx);
    x[0] -= s * dx[0];
    x[1] -= s * dx[1];
  }
  bool converged(petsc::Vector &) override
  {
    const double err = std::hypot(x[0] - 1.0, x[1] - 1.0);
    history.push_back(err);
    return err < 1e-12;
  }
  uint size() override { return 2; }
};

//! F(x) = atan(x). The full Newton step diverges for |x0| > 1.39
class Arctan : public NonlinearProblem
{
public:
  double x;
  explicit Arctan(double x0) : x(x0) {}

  void evaluate(petsc::Vector & r) override { r.set(0, std::atan(x)); r.assemble(); }
  void jacobian(petsc::Matrix & K) override { K.set(0, 0, 1.0 / (1.0 + x*x)); K.assemble(); }
  void update(petsc::Vector & u, double s) override
  {
    double dx;
    u.get_data(&dx);
    x -= s * dx;
  }
  bool converged(petsc::Vector &) override { return std::fabs(x) < 1e-10; }
  uint size() override { return 1; }
};

//! F(x) = x^2 + 1 has no real root
class NoRoot : public NonlinearProblem
{
public:
  double x = 3.0;
  void evaluate(petsc::Vector & r) override { r.set(0, x*x + 1.0); r.assemble(); }
  void jacobian(petsc::Matrix & K) override { K.set(0, 0, 2.0*x); K.assemble(); }
  void update(petsc::Vector & u, double s) override
  {
    double dx;
    u.get_data(&dx);
    x -= s * dx;
  }
  bool converged(petsc::Vector &) override { return false; }
  uint size() override { return 1; }
};

} // namespace

TEST_CASE("Newton: converges to the root of the Kelley problem")
{
  Kelley p;
  Newton solver(&p);
  solver.init();
  const int its = solver.solve();
  CHECK_NEAR(p.x[0], 1.0, 1e-12);
  CHECK_NEAR(p.x[1], 1.0, 1e-12);
  CHECK(its > 1);
  CHECK(its < 15);
}

TEST_CASE("Newton: the error decreases quadratically near the root")
{
  Kelley p;
  Newton solver(&p);
  solver.init();
  solver.solve();

  // e_{k+1} ~ C e_k^2: measure the order where the errors are small but not at roundoff
  double order = 0.0;
  for (std::size_t k = 1; k + 1 < p.history.size(); k++)
  {
    const double e0 = p.history[k-1], e1 = p.history[k], e2 = p.history[k+1];
    if (e1 < 1e-1 && e2 > 1e-13 && e0 > e1)
      order = std::log(e2 / e1) / std::log(e1 / e0);
  }
  CHECK_MSG(order > 1.7, "observed order " << order);
}

TEST_CASE("Newton: gives up with an error when it cannot converge")
{
  NoRoot p;
  Newton solver(&p);
  solver.init();
  CHECK_THROWS(solver.solve());
}

TEST_CASE("Newton: the plain step diverges on arctan from far away")
{
  // This is the premise of the line search test below
  Arctan p(3.0);
  Newton solver(&p);
  solver.init();
  bool failed = false;
  try { solver.solve(); } catch (...) { failed = true; }
  CHECK(failed || !std::isfinite(p.x));
}

TEST_CASE("NewtonLineSearch: reaches the root of the Kelley problem")
{
  Kelley p;
  NewtonLineSearch solver(&p);
  solver.init();
  solver.solve();
  CHECK_NEAR(p.x[0], 1.0, 1e-10);
  CHECK_NEAR(p.x[1], 1.0, 1e-10);
}

TEST_CASE("NewtonLineSearch: converges on arctan, where the plain step diverges")
{
  Arctan p(3.0);
  NewtonLineSearch solver(&p);
  solver.init();
  solver.solve();
  CHECK_NEAR(p.x, 0.0, 1e-10);
}
