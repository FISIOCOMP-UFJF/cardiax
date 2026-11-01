/*
 * Unit tests of src/linalg: PETSc wrappers, the backend selection and the 4th order
 * tensors of the mechanics.
 *
 * What is covered and why:
 *  - petsc::Vector: creation, set/add, assembly, dot/norm, scale, get/set of raw data,
 *    copy. Every solver stores its solution here.
 *  - petsc::Matrix: assembly, matvec, the symmetric flag, zero_rows_cols (how every
 *    Dirichlet boundary condition is applied, including the right-hand side
 *    correction), and get_CSR (the path of the matrix into AMGX).
 *  - petsc::LinearSolver: direct and iterative solves on a known system, repeated
 *    solves, and solve_constant_matrix (used by the EP solvers, where the matrix is
 *    factorised once).
 *  - Tensor4: index layout, the tensor product and unit tensor used in the elasticity
 *    tangents, and the Voigt conversion in get_matrix().
 *
 * The AMGX backend has its own test (test_amgx.cpp, label gpu).
 */
#include "test_framework.hpp"
#include "fixtures.hpp"

#include "linalg/petsc_vector.hpp"
#include "linalg/petsc_matrix.hpp"
#include "linalg/petsc_linear_solver.hpp"
#include "linalg/tensor.hpp"

using namespace cardiax_test;

namespace
{

// Symmetric positive definite tridiagonal 1D Laplacian, with its exact solution
void laplacian_1d(int n, petsc::Matrix & A)
{
  A.create(n, n, 3);
  for (int i = 0; i < n; i++)
  {
    A.set(i, i, 2.0);
    if (i > 0)     A.set(i, i - 1, -1.0);
    if (i < n - 1) A.set(i, i + 1, -1.0);
  }
  A.assemble();
}

arma::vec to_arma(petsc::Vector & v)
{
  arma::vec a(v.size());
  v.get_data(a.memptr());
  return a;
}

} // namespace

// ---- Vector -----------------------------------------------------------------

TEST_CASE("Vector: size, set, add and assemble")
{
  petsc::Vector v;
  v.create(5);
  CHECK_EQ(v.size(), 5);

  v.set(1.0);
  v.set(2, 3.0);
  v.assemble();
  v.add(2, 0.5);
  v.add(4, 2.0);
  v.assemble();

  const arma::vec a = to_arma(v);
  CHECK_NEAR(a(0), 1.0, 0.0);
  CHECK_NEAR(a(2), 3.5, 0.0);
  CHECK_NEAR(a(4), 3.0, 0.0);
}

TEST_CASE("Vector: indexed set and add of several entries")
{
  petsc::Vector v;
  v.create(4);
  v = 0.0;
  const int idx[2] = {1, 3};
  const double val[2] = {2.0, 4.0};
  v.set(2, idx, val);
  v.assemble();
  v.add(2, idx, val);
  v.assemble();
  const arma::vec a = to_arma(v);
  CHECK_NEAR(a(0), 0.0, 0.0);
  CHECK_NEAR(a(1), 4.0, 0.0);
  CHECK_NEAR(a(3), 8.0, 0.0);
}

TEST_CASE("Vector: dot, l2norm, scale and add_vec")
{
  petsc::Vector x, y;
  x.create(3);
  y.create(3);
  x.set(2.0);
  y.set(3.0);
  CHECK_NEAR(x.dot(y), 18.0, 1e-14);
  CHECK_NEAR(x.l2norm(), std::sqrt(12.0), 1e-14);

  x.scale(0.5);
  CHECK_NEAR(x.l2norm(), std::sqrt(3.0), 1e-14);

  x.add_vec(y);
  CHECK_NEAR(to_arma(x)(1), 4.0, 1e-14);
}

TEST_CASE("Vector: raw data in and out, whole and by range")
{
  petsc::Vector v;
  v.create(4);
  const double in[4] = {1.0, 2.0, 3.0, 4.0};
  v.set_data(in);

  double out[4] = {0, 0, 0, 0};
  v.get_data(out);
  for (int i = 0; i < 4; i++) CHECK_NEAR(out[i], in[i], 0.0);

  double * a = v.get_array();
  a[0] = 10.0;
  v.restore_array(a);
  CHECK_NEAR(to_arma(v)(0), 10.0, 0.0);
}

TEST_CASE("Vector: copy_values copies the first n values")
{
  petsc::Vector a, b;
  a.create(3);
  b.create(3);
  a.set(7.0);
  b.set(1.0);
  a.copy_values(3, b);
  CHECK_NEAR(to_arma(a)(2), 1.0, 0.0);
}

TEST_CASE("Vector: place_array makes the vector use external memory")
{
  petsc::Vector v;
  v.create(3);
  double ext[3] = {5.0, 6.0, 7.0};
  v.place_array(ext);
  CHECK_NEAR(v.l2norm(), std::sqrt(25.0 + 36.0 + 49.0), 1e-14);
  v.reset_array();
}

// ---- Matrix -----------------------------------------------------------------

TEST_CASE("Matrix: assembly and matvec")
{
  petsc::Matrix A;
  laplacian_1d(4, A);
  CHECK_EQ(A.size(), 4u);
  CHECK_EQ(A.get_nnz(), 10u);
  CHECK(A.is_sym());
  CHECK(!A.is_null());

  petsc::Vector x, y;
  x.create(4);
  y.create(4);
  const double xv[4] = {1, 2, 3, 4};
  x.set_data(xv);
  A.mult(x, y);
  const arma::vec r = to_arma(y);
  CHECK_NEAR(r(0), 0.0, 1e-14);
  CHECK_NEAR(r(1), 0.0, 1e-14);
  CHECK_NEAR(r(2), 0.0, 1e-14);
  CHECK_NEAR(r(3), 5.0, 1e-14);
}

TEST_CASE("Matrix: add accumulates, including block add")
{
  petsc::Matrix A;
  A.create(2, 2, 2);
  A = 0.0;
  const int rows[2] = {0, 1};
  const double block[4] = {1.0, 2.0, 3.0, 4.0};
  A.add(2, 2, rows, rows, block);
  A.assemble();
  A.add(0, 0, 10.0);
  A.assemble();

  petsc::Vector x, y;
  x.create(2);
  y.create(2);
  x.set(1.0);
  A.mult(x, y);
  const arma::vec r = to_arma(y);
  CHECK_NEAR(r(0), 13.0, 1e-14);
  CHECK_NEAR(r(1), 7.0, 1e-14);
  CHECK(!A.is_sym());
}

TEST_CASE("Matrix: operator= zeroes the entries and keeps the sparsity")
{
  petsc::Matrix A;
  laplacian_1d(3, A);
  const uint nnz = A.get_nnz();
  A = 0.0;
  A.assemble();
  CHECK_EQ(A.get_nnz(), nnz);

  petsc::Vector x, y;
  x.create(3);
  y.create(3);
  x.set(1.0);
  A.mult(x, y);
  CHECK_NEAR(y.l2norm(), 0.0, 0.0);
}

TEST_CASE("Matrix: zero_rows_cols keeps symmetry and moves the column into the RHS")
{
  // Dirichlet u_0 = 2 on a 3x3 Laplacian. The eliminated column must be moved to
  // the right-hand side, or the rest of the solution would be wrong.
  petsc::Matrix A;
  laplacian_1d(3, A);
  petsc::Vector x, b;
  x.create(3);
  b.create(3);
  x = 0.0;
  x.set(0, 2.0);
  x.assemble();
  b = 0.0;

  int row = 0;
  A.zero_rows_cols(1, &row, 1.0, x, b);
  CHECK(A.is_sym());

  const arma::vec rhs = to_arma(b);
  CHECK_NEAR(rhs(0), 2.0, 1e-14);
  CHECK_NEAR(rhs(1), 2.0, 1e-14);   // -A(1,0)*u_0
  CHECK_NEAR(rhs(2), 0.0, 1e-14);

  petsc::LinearSolver ls;
  ls.init();
  petsc::Vector sol;
  sol.create(3);
  sol = 0.0;
  ls.solve(A, sol, b);
  const arma::vec u = to_arma(sol);
  // 2u1 - u2 = 2, -u1 + 2u2 = 0
  CHECK_NEAR(u(0), 2.0, 1e-12);
  CHECK_NEAR(u(1), 4.0 / 3.0, 1e-12);
  CHECK_NEAR(u(2), 2.0 / 3.0, 1e-12);
}

TEST_CASE("Matrix: get_CSR describes the same matrix")
{
  petsc::Matrix A;
  laplacian_1d(3, A);
  int n = 0, ia[4], ja[7];
  double v[7];
  A.get_CSR(&n, ia, ja, v);
  CHECK_EQ(n, 3);

  // Whatever the index base, the row pointers span the 7 nonzeros
  CHECK_EQ(ia[3] - ia[0], 7);
  const int base = ia[0];
  arma::mat dense(3, 3, arma::fill::zeros);
  for (int i = 0; i < 3; i++)
    for (int k = ia[i] - base; k < ia[i + 1] - base; k++)
      dense(i, ja[k] - base) = v[k];

  const arma::mat expected = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
  CHECK_NEAR(arma::abs(dense - expected).max(), 0.0, 0.0);
}

// ---- LinearSolver -----------------------------------------------------------

TEST_CASE("LinearSolver: the direct solve recovers a manufactured solution")
{
  const int n = 50;
  petsc::Matrix A;
  laplacian_1d(n, A);
  arma::vec exact = arma::randu<arma::vec>(n);
  petsc::Vector x, b, xe;
  x.create(n);
  b.create(n);
  xe.create(n);
  xe.set_data(exact.memptr());
  A.mult(xe, b);
  x = 0.0;

  petsc::LinearSolver ls;
  ls.init();
  CHECK(!ls.is_null());
  ls.solve(A, x, b);
  CHECK_NEAR(arma::abs(to_arma(x) - exact).max(), 0.0, 1e-11);
}

TEST_CASE("LinearSolver: CG with Jacobi reaches the same solution through the options")
{
  const int n = 30;
  petsc::Matrix A;
  laplacian_1d(n, A);
  arma::vec exact(n);
  for (int i = 0; i < n; i++) exact(i) = std::sin(0.3 * i);
  petsc::Vector x, b, xe;
  x.create(n);
  b.create(n);
  xe.create(n);
  xe.set_data(exact.memptr());
  A.mult(xe, b);
  x = 0.0;

  PetscOptionsSetValue(nullptr, "-ksp_type", "cg");
  PetscOptionsSetValue(nullptr, "-pc_type", "jacobi");
  PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-12");
  {
    petsc::LinearSolver ls;
    ls.init();
    const auto res = ls.solve(A, x, b, 1e-12);
    CHECK(res.first > 1);    // iterated: the options were honoured
  }
  PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
  PetscOptionsSetValue(nullptr, "-pc_type", "lu");
  PetscOptionsClearValue(nullptr, "-ksp_rtol");

  CHECK_NEAR(arma::abs(to_arma(x) - exact).max(), 0.0, 1e-8);
}

TEST_CASE("LinearSolver: solve_constant_matrix can be called repeatedly with new RHS")
{
  const int n = 20;
  petsc::Matrix A;
  laplacian_1d(n, A);
  petsc::LinearSolver ls;
  ls.init();

  for (double scale : {1.0, -2.0, 3.5})
  {
    arma::vec exact(n);
    for (int i = 0; i < n; i++) exact(i) = scale * (i + 1);
    petsc::Vector x, b, xe;
    x.create(n);
    b.create(n);
    xe.create(n);
    xe.set_data(exact.memptr());
    A.mult(xe, b);
    x = 0.0;
    ls.solve_constant_matrix(A, x, b);
    CHECK_NEAR(arma::abs(to_arma(x) - exact).max(), 0.0, 1e-10);
  }
}

TEST_CASE("LinearSolver: the same solver follows a matrix that changes between solves")
{
  // The Newton loop does this: same LinearSolver, new tangent values
  const int n = 10;
  petsc::Matrix A;
  laplacian_1d(n, A);
  petsc::LinearSolver ls;
  ls.init();

  petsc::Vector x, b;
  x.create(n);
  b.create(n);
  b.set(1.0);
  x = 0.0;
  ls.solve(A, x, b);
  const arma::vec first = to_arma(x);

  for (int i = 0; i < n; i++) A.set(i, i, 4.0);
  A.assemble();
  x = 0.0;
  ls.solve(A, x, b);
  const arma::vec second = to_arma(x);

  // Check against the new matrix, not against a stale factorisation
  petsc::Vector Ax;
  Ax.create(n);
  A.mult(x, Ax);
  CHECK_NEAR(arma::abs(to_arma(Ax) - to_arma(b)).max(), 0.0, 1e-11);
  CHECK(arma::abs(first - second).max() > 1e-3);
}

// ---- Tensor4 ----------------------------------------------------------------

TEST_CASE("Tensor4: zero and component access")
{
  Tensor4 T;
  T.zero();
  T(0, 1, 2, 0) = 3.0;
  CHECK_NEAR(T(0, 1, 2, 0), 3.0, 0.0);
  CHECK_NEAR(T(1, 0, 2, 0), 0.0, 0.0);
  CHECK_NEAR(T(0, 1, 0, 2), 0.0, 0.0);
}

TEST_CASE("Tensor4: tensor_product is A_ij B_kl")
{
  arma::mat33 A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  arma::mat33 B = {{2, 0, 1}, {0, 3, 0}, {1, 0, 4}};
  Tensor4 T = tensor_product(A, B);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          CHECK_NEAR(T(i, j, k, l), A(i, j) * B(k, l), 0.0);
}

TEST_CASE("Tensor4: unit_tensor applied to a matrix returns the matrix")
{
  Tensor4 I = unit_tensor();
  arma::mat33 A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      double r = 0.0;
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          r += I(i, j, k, l) * A(k, l);
      CHECK_NEAR(r, A(i, j), 0.0);
    }
}

TEST_CASE("Tensor4: 3D get_matrix gives the Voigt matrix of the isotropic tensor")
{
  // C_ijkl = lam d_ij d_kl + mu (d_ik d_jl + d_il d_jk)
  const double lam = 2.0, mu = 3.0;
  arma::mat33 d = arma::eye(3, 3);
  Tensor4 C;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          C(i, j, k, l) = lam * d(i, j) * d(k, l) + mu * (d(i, k) * d(j, l) + d(i, l) * d(j, k));

  arma::mat D(6, 6, arma::fill::zeros);
  C.get_matrix(D);

  CHECK_NEAR(D(0, 0), lam + 2 * mu, 1e-14);
  CHECK_NEAR(D(1, 1), lam + 2 * mu, 1e-14);
  CHECK_NEAR(D(2, 2), lam + 2 * mu, 1e-14);
  CHECK_NEAR(D(0, 1), lam, 1e-14);
  CHECK_NEAR(D(1, 2), lam, 1e-14);
  CHECK_NEAR(D(3, 3), mu, 1e-14);
  CHECK_NEAR(D(4, 4), mu, 1e-14);
  CHECK_NEAR(D(5, 5), mu, 1e-14);
  CHECK_NEAR(D(0, 3), 0.0, 1e-14);
  CHECK_NEAR(arma::abs(D - D.t()).max(), 0.0, 1e-14);
}

TEST_CASE_KNOWN_BUG("Tensor4: 2D get_matrix gives the Voigt matrix of the isotropic tensor",
                    "src/linalg/tensor.cpp: D(2,0) = D(2,0) and D(2,1) = D(2,1) assign the "
                    "entries to themselves; the third row keeps whatever D held before")
{
  const double lam = 2.0, mu = 3.0;
  arma::mat33 d = arma::eye(3, 3);
  Tensor4 C;
  // A tensor with nonzero (xy, xx) coupling, as an anisotropic material has
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          C(i, j, k, l) = lam * d(i, j) * d(k, l) + mu * (d(i, k) * d(j, l) + d(i, l) * d(j, k));
  C(0, 1, 0, 0) = C(1, 0, 0, 0) = 1.5;

  arma::mat D(3, 3);
  D.fill(99.0);    // what is left in the output must not matter
  C.get_matrix(D);
  CHECK_NEAR(D(0, 0), lam + 2 * mu, 1e-14);
  CHECK_NEAR(D(0, 1), lam, 1e-14);
  CHECK_NEAR(D(2, 2), mu, 1e-14);
  CHECK_NEAR(D(2, 0), 1.5, 1e-14);
  CHECK_NEAR(D(2, 1), 0.0, 1e-14);
}
